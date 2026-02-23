#!/usr/bin/env Rscript

# =============================================================================
# RNA-seq DESeq2 Analysis Script
# Called by Snakemake rule: deseq2
# Input: merged gene count matrix from featureCounts + metadata + comparisons
# Output: QC plots, DE tables, volcano/MA plots, heatmaps, GO enrichment, PCA
# =============================================================================

start_time <- Sys.time()
cat("=== DESeq2 Analysis Script ===\n")
cat("Start time:", format(start_time), "\n\n")

suppressPackageStartupMessages({
    library(DESeq2)
    library(ggplot2)
    library(ggrepel)
    library(pheatmap)
    library(RColorBrewer)
    library(dplyr)
    library(tidyr)
    library(reshape2)
    library(scales)
    library(gridExtra)
    library(clusterProfiler)
    library(enrichplot)
    library(AnnotationDbi)
    library(data.table)
})
cat("All libraries loaded.\n")

# =============================================================================
# Argument Parsing
# =============================================================================
args <- commandArgs(trailingOnly = TRUE)

parse_args <- function(args) {
    args_list <- list()
    i <- 1
    while (i <= length(args)) {
        if (startsWith(args[i], "--")) {
            key <- sub("^--", "", args[i])
            if (i + 1 <= length(args) && !startsWith(args[i + 1], "--")) {
                args_list[[key]] <- args[i + 1]
                i <- i + 2
            } else {
                args_list[[key]] <- TRUE
                i <- i + 1
            }
        } else {
            i <- i + 1
        }
    }
    return(args_list)
}

opts <- parse_args(args)

counts_file      <- opts$counts
metadata_file    <- opts$metadata
comparisons_file <- opts$comparisons
outdir           <- opts$outdir
species          <- opts$species
lfc_threshold    <- as.numeric(opts$lfc_threshold)
padj_threshold   <- as.numeric(opts$padj_threshold)
n_threads        <- as.integer(opts$threads)

cat("Counts file:      ", counts_file, "\n")
cat("Metadata file:    ", metadata_file, "\n")
cat("Comparisons file: ", comparisons_file, "\n")
cat("Output directory: ", outdir, "\n")
cat("Species:          ", species, "\n")
cat("LFC threshold:    ", lfc_threshold, "\n")
cat("padj threshold:   ", padj_threshold, "\n")
cat("Threads:          ", n_threads, "\n\n")

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# Species-Specific OrgDb + ID Mapping
# =============================================================================
OrgDb_species <- NULL
mapping_table <- NULL

if (species == "mouse") {
    if (requireNamespace("org.Mm.eg.db", quietly = TRUE)) {
        library(org.Mm.eg.db)
        OrgDb_species <- org.Mm.eg.db
        cat("Loaded org.Mm.eg.db\n")
    } else { message("Warning: org.Mm.eg.db not installed.") }
} else if (species == "human") {
    if (requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
        library(org.Hs.eg.db)
        OrgDb_species <- org.Hs.eg.db
        cat("Loaded org.Hs.eg.db\n")
    } else { message("Warning: org.Hs.eg.db not installed.") }
} else {
    cat("Species '", species, "' - no OrgDb available. GO enrichment will be skipped.\n")
}

if (!is.null(OrgDb_species)) {
    tryCatch({
        ensembl_keys <- keys(OrgDb_species, keytype = "ENSEMBL")
        if (length(ensembl_keys) > 0) {
            mapping_raw <- AnnotationDbi::select(OrgDb_species, keys = ensembl_keys,
                                                  columns = c("ENTREZID", "SYMBOL"), keytype = "ENSEMBL")
            mapping_dt <- data.table(mapping_raw)
            mapping_table <- mapping_dt[, lapply(.SD, function(x) first(na.omit(x))), by = ENSEMBL]
            colnames(mapping_table) <- c("ensembl_id", "entrez_id", "symbol")
            cat("ID mapping table created:", nrow(mapping_table), "entries\n")
        }
    }, error = function(e) { message("Warning: Could not create ID mapping table: ", e$message) })
}

get_gene_symbol <- function(ensembl_id) {
    if (is.null(mapping_table)) return(NA_character_)
    idx <- match(ensembl_id, mapping_table$ensembl_id)
    if (is.na(idx)) return(NA_character_)
    return(mapping_table$symbol[idx])
}

# Vectorized version for efficiency
get_gene_symbols <- function(ensembl_ids) {
    if (is.null(mapping_table)) return(rep(NA_character_, length(ensembl_ids)))
    idx <- match(ensembl_ids, mapping_table$ensembl_id)
    symbols <- mapping_table$symbol[idx]
    return(as.character(symbols))
}

# =============================================================================
# Load Count Matrix
# =============================================================================
cat("\n--- Loading Count Matrix ---\n")
counts_raw <- read.delim(counts_file, header = TRUE, row.names = 1, check.names = FALSE)

gene_lengths <- counts_raw$Length
names(gene_lengths) <- rownames(counts_raw)
count_matrix <- as.matrix(counts_raw[, !colnames(counts_raw) %in% "Length", drop = FALSE])

# Strip Ensembl version numbers
rownames(count_matrix) <- gsub("\\..*", "", rownames(count_matrix))
names(gene_lengths) <- gsub("\\..*", "", names(gene_lengths))

cat("Count matrix dimensions:", nrow(count_matrix), "genes x", ncol(count_matrix), "samples\n")
cat("Samples:", paste(colnames(count_matrix), collapse = ", "), "\n")

# =============================================================================
# Load Metadata + Match Samples
# =============================================================================
cat("\n--- Loading Metadata ---\n")
metadata <- read.delim(metadata_file, header = TRUE, stringsAsFactors = FALSE)
cat("Metadata samples:", paste(metadata$sample_name, collapse = ", "), "\n")

required_cols <- c("sample_name", "replicate", "condition")
missing_cols <- setdiff(required_cols, colnames(metadata))
if (length(missing_cols) > 0) stop("FATAL: Missing metadata columns: ", paste(missing_cols, collapse = ", "))

common_samples <- intersect(colnames(count_matrix), metadata$sample_name)
if (length(common_samples) == 0) stop("FATAL: No matching sample names between count matrix and metadata!")
cat("Matched samples:", length(common_samples), "/", ncol(count_matrix), "\n")

count_matrix <- count_matrix[, common_samples, drop = FALSE]
metadata <- metadata[metadata$sample_name %in% common_samples, , drop = FALSE]
metadata <- metadata[match(colnames(count_matrix), metadata$sample_name), ]
metadata$condition <- factor(metadata$condition)
metadata$replicate <- factor(metadata$replicate)

if (!identical(colnames(count_matrix), metadata$sample_name)) stop("FATAL: Sample order mismatch!")

# =============================================================================
# RPKM / TPM Calculation
# =============================================================================
cat("\n--- Calculating RPKM and TPM ---\n")
geneRPKM <- NULL
geneTPM <- NULL
tryCatch({
    valid_genes <- gene_lengths[rownames(count_matrix)] > 0
    gl <- gene_lengths[rownames(count_matrix)][valid_genes]
    cm <- count_matrix[valid_genes, , drop = FALSE]

    rpk <- cm / (gl / 1000)
    lib_millions <- colSums(cm) / 1e6
    lib_millions[lib_millions == 0] <- NA

    geneRPKM <- t(t(rpk) / lib_millions)
    colnames(geneRPKM) <- paste0(colnames(geneRPKM), "_RPKM")

    rpk_sums <- colSums(rpk, na.rm = TRUE) / 1e6
    rpk_sums[rpk_sums == 0] <- NA
    geneTPM <- t(t(rpk) / rpk_sums)
    colnames(geneTPM) <- paste0(colnames(geneTPM), "_TPM")

    outtable <- data.frame(ensembl_id = rownames(geneRPKM), geneRPKM, geneTPM, check.names = FALSE)
    if (!is.null(mapping_table)) {
        outtable <- merge(outtable, as.data.frame(mapping_table), by = "ensembl_id", all.x = TRUE)
    }
    write.table(outtable, file = file.path(outdir, "geneRPKM_TPM_table.tsv"),
                row.names = FALSE, quote = FALSE, sep = "\t")
    cat("RPKM/TPM table saved.\n")
}, error = function(e) { message("Warning: RPKM/TPM calculation failed: ", e$message) })

# Export raw counts with gene IDs
cnt_out <- data.frame(gene_id = rownames(count_matrix), count_matrix, check.names = FALSE)
write.table(cnt_out, file = file.path(outdir, "geneCOUNT_table.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("Raw counts table saved.\n")

# =============================================================================
# QC Plots
# =============================================================================
cat("\n--- Generating QC Plots ---\n")

# --- Histograms ---
tryCatch({
    row_sums <- rowSums(count_matrix, na.rm = TRUE)
    num_zeros <- sum(count_matrix == 0, na.rm = TRUE)
    perc_zeros <- (num_zeros / length(count_matrix)) * 100

    p_raw <- ggplot(data.frame(counts = log1p(row_sums)), aes(x = counts)) +
        geom_histogram(binwidth = 0.1, fill = "blue", color = "black") +
        labs(title = "Histogram of log1p(Raw Counts)",
             x = "Log1p-transformed Total Counts per Gene", y = "Frequency") +
        theme_minimal() +
        annotate("text", x = Inf, y = Inf, label = sprintf("%% Zeros: %.1f%%", perc_zeros),
                 hjust = 1.1, vjust = 1.5, size = 4, color = "red")
    ggsave(file.path(outdir, "histogram_raw_counts.png"), p_raw, width = 10, height = 6, dpi = 300)

    # Filtered: genes expressed in >= 20% of samples
    prop_expr <- rowSums(count_matrix > 0) / ncol(count_matrix)
    keep <- prop_expr >= 0.2
    if (sum(keep) > 0) {
        filt_mat <- count_matrix[keep, , drop = FALSE]
        filt_sums <- rowSums(filt_mat, na.rm = TRUE)
        filt_perc <- sum(filt_mat == 0) / length(filt_mat) * 100

        p_filt <- ggplot(data.frame(counts = log1p(filt_sums)), aes(x = counts)) +
            geom_histogram(binwidth = 0.1, fill = "blue", color = "black") +
            labs(title = "Histogram of log1p(Filtered Counts)",
                 subtitle = "Genes expressed > 0 in at least 20% of samples",
                 x = "Log1p-transformed Total Counts per Gene", y = "Frequency") +
            theme_minimal() +
            annotate("text", x = Inf, y = Inf, label = sprintf("%% Zeros: %.1f%%", filt_perc),
                     hjust = 1.1, vjust = 1.5, size = 4, color = "red") +
            annotate("text", x = Inf, y = Inf,
                     label = sprintf("Genes retained: %d / %d", sum(keep), nrow(count_matrix)),
                     hjust = 1.1, vjust = 3, size = 4, color = "blue")
        ggsave(file.path(outdir, "histogram_filtered_counts.png"), p_filt, width = 10, height = 6, dpi = 300)
    }
    cat("Histograms saved.\n")
}, error = function(e) { message("ERROR generating histograms: ", e$message) })

# --- RPKM/TPM Boxplots and Library Size ---
if (!is.null(geneRPKM) && !is.null(geneTPM)) {
    tryCatch({
        pseudo <- 0.1

        plot_rpkm <- reshape::melt(as.data.frame(geneRPKM), variable.name = "Sample", value.name = "RPKM")
        plot_rpkm$Sample <- gsub("_RPKM$", "", plot_rpkm$Sample)
        plot_rpkm$log2val <- log2(plot_rpkm$RPKM + pseudo)

        plot_tpm <- reshape::melt(as.data.frame(geneTPM), variable.name = "Sample", value.name = "TPM")
        plot_tpm$Sample <- gsub("_TPM$", "", plot_tpm$Sample)
        plot_tpm$log2val <- log2(plot_tpm$TPM + pseudo)

        p1 <- ggplot(plot_rpkm, aes(x = Sample, y = log2val)) +
            geom_boxplot(outlier.shape = NA) + theme_bw() +
            theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
            labs(title = "log2(RPKM + 0.1)", y = "log2(RPKM + 0.1)") +
            coord_cartesian(ylim = quantile(plot_rpkm$log2val, c(0.01, 0.99), na.rm = TRUE))

        p2 <- ggplot(plot_tpm, aes(x = Sample, y = log2val)) +
            geom_boxplot(outlier.shape = NA) + theme_bw() +
            theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
            labs(title = "log2(TPM + 0.1)", y = "log2(TPM + 0.1)") +
            coord_cartesian(ylim = quantile(plot_tpm$log2val, c(0.01, 0.99), na.rm = TRUE))

        p_combined <- grid.arrange(p1, p2, ncol = 2)
        ggsave(file.path(outdir, "rpkm_tpm_distribution_boxplot.png"), p_combined,
               width = 20, height = 7, dpi = 300)
        cat("RPKM/TPM boxplots saved.\n")

        # Library size barplot
        lib_rpkm <- data.frame(Sample = gsub("_RPKM$", "", names(colSums(geneRPKM, na.rm = TRUE))),
                                Size = colSums(geneRPKM, na.rm = TRUE), Type = "RPKM")
        lib_tpm  <- data.frame(Sample = gsub("_TPM$", "", names(colSums(geneTPM, na.rm = TRUE))),
                                Size = colSums(geneTPM, na.rm = TRUE), Type = "TPM")
        lib_df <- rbind(lib_rpkm, lib_tpm)

        p3 <- ggplot(lib_df, aes(x = Sample, y = Size, fill = Type)) +
            geom_bar(stat = "identity", position = position_dodge()) +
            theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
            labs(title = "Library sizes across samples", x = "Sample", y = "Total RPKM/TPM") +
            scale_y_continuous(labels = comma)
        ggsave(file.path(outdir, "library_sizes_plot.png"), p3, width = 15, height = 6, dpi = 300)
        cat("Library size plot saved.\n")

        # Top 20 variable genes heatmap
        tpm_matrix <- as.matrix(geneTPM)
        row_labels <- rownames(tpm_matrix)
        symbols <- get_gene_symbols(gsub("\\..*", "", row_labels))
        rownames(tpm_matrix) <- ifelse(is.na(symbols) | symbols == "", row_labels,
                                        paste0(symbols, " (", row_labels, ")"))

        gene_var <- apply(tpm_matrix, 1, var, na.rm = TRUE)
        top_20 <- order(gene_var, decreasing = TRUE, na.last = NA)[1:min(20, sum(!is.na(gene_var)))]
        if (length(top_20) > 0) {
            log2_mat <- log2(tpm_matrix[top_20, , drop = FALSE] + pseudo)
            colnames(log2_mat) <- gsub("_TPM$", "", colnames(log2_mat))
            png(file.path(outdir, "top_20_variable_genes_heatmap.png"), width = 1400, height = 1000, res = 150)
            pheatmap(log2_mat, scale = "row", main = "Top 20 Variable Genes (log2(TPM+0.1))",
                     fontsize_row = 7, fontsize_col = 9,
                     color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100))
            dev.off()
            cat("Top variable genes heatmap saved.\n")
        }
    }, error = function(e) { message("ERROR generating RPKM/TPM plots: ", e$message) })
}

# =============================================================================
# Helper Functions for DE Plots
# =============================================================================

create_ma_plot <- function(res, title, filepath) {
    tryCatch({
        res_df <- as.data.frame(res)
        res_df$significant <- !is.na(res_df$padj) &
            res_df$padj < padj_threshold &
            abs(res_df$log2FoldChange) > lfc_threshold

        n_up   <- sum(res_df$significant & res_df$log2FoldChange > 0, na.rm = TRUE)
        n_down <- sum(res_df$significant & res_df$log2FoldChange < 0, na.rm = TRUE)

        p <- ggplot(res_df, aes(x = log10(baseMean + 1), y = log2FoldChange, color = significant)) +
            geom_point(alpha = 0.4, size = 0.5) +
            scale_color_manual(values = c("grey40", "red")) +
            geom_hline(yintercept = c(-lfc_threshold, lfc_threshold), linetype = "dashed", color = "blue") +
            geom_hline(yintercept = 0, color = "black", linewidth = 0.3) +
            labs(title = title,
                 subtitle = sprintf("Up: %d | Down: %d (padj < %s, |LFC| > %s)",
                                    n_up, n_down, padj_threshold, lfc_threshold),
                 x = "log10(Mean Expression + 1)", y = "log2 Fold Change") +
            theme_minimal() + theme(legend.position = "none")
        ggsave(filepath, p, width = 10, height = 7, dpi = 300)
        cat("  MA plot saved.\n")
    }, error = function(e) { message("  Error creating MA plot: ", e$message) })
}

create_volcano_plot <- function(res, title, filepath) {
    tryCatch({
        res_df <- as.data.frame(res)
        res_df <- res_df[!is.na(res_df$padj), ]
        res_df$neg_log10_padj <- -log10(res_df$padj)
        res_df$category <- "NS"
        res_df$category[res_df$padj < padj_threshold & res_df$log2FoldChange > lfc_threshold] <- "Up"
        res_df$category[res_df$padj < padj_threshold & res_df$log2FoldChange < -lfc_threshold] <- "Down"
        res_df$category <- factor(res_df$category, levels = c("NS", "Up", "Down"))

        n_up   <- sum(res_df$category == "Up")
        n_down <- sum(res_df$category == "Down")

        # Label top significant genes
        res_df$ensembl_id <- rownames(res_df)
        res_df$symbol <- get_gene_symbols(res_df$ensembl_id)
        res_df$label <- ""
        sig_rows <- res_df[res_df$category != "NS", ]
        if (nrow(sig_rows) > 0) {
            top_n <- min(20, nrow(sig_rows))
            top_genes <- sig_rows[order(sig_rows$padj), ][1:top_n, ]
            res_df$label[rownames(res_df) %in% rownames(top_genes)] <-
                ifelse(!is.na(top_genes$symbol) & top_genes$symbol != "",
                       top_genes$symbol, top_genes$ensembl_id)
        }

        p <- ggplot(res_df, aes(x = log2FoldChange, y = neg_log10_padj, color = category)) +
            geom_point(alpha = 0.4, size = 0.8) +
            scale_color_manual(values = c("NS" = "grey60", "Up" = "red", "Down" = "blue")) +
            geom_vline(xintercept = c(-lfc_threshold, lfc_threshold), linetype = "dashed", alpha = 0.5) +
            geom_hline(yintercept = -log10(padj_threshold), linetype = "dashed", alpha = 0.5) +
            geom_text_repel(aes(label = label), size = 2.5, max.overlaps = 20, show.legend = FALSE) +
            labs(title = title,
                 subtitle = sprintf("Up: %d | Down: %d (padj < %s, |LFC| > %s)",
                                    n_up, n_down, padj_threshold, lfc_threshold),
                 x = "log2 Fold Change", y = "-log10(adjusted p-value)") +
            theme_minimal() + theme(legend.position = "bottom")
        ggsave(filepath, p, width = 10, height = 8, dpi = 300)
        cat("  Volcano plot saved.\n")
    }, error = function(e) { message("  Error creating volcano plot: ", e$message) })
}

create_heatmap <- function(mat, col_data, title, filepath) {
    tryCatch({
        # Add gene symbols as row labels
        symbols <- get_gene_symbols(rownames(mat))
        row_labels <- ifelse(is.na(symbols) | symbols == "", rownames(mat),
                             paste0(symbols, " (", rownames(mat), ")"))
        rownames(mat) <- row_labels

        annotation_col <- data.frame(Condition = col_data$condition,
                                      row.names = col_data$sample_name)

        png(filepath, width = 1200, height = max(400, nrow(mat) * 18 + 200), res = 150)
        pheatmap(mat, scale = "row", annotation_col = annotation_col,
                 main = title, fontsize_row = 6, fontsize_col = 8,
                 color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100),
                 clustering_distance_rows = "euclidean", clustering_method = "complete")
        dev.off()
        cat("  Heatmap saved.\n")
    }, error = function(e) { message("  Error creating heatmap: ", e$message) })
}

perform_GO_enrichment <- function(res, comp_name, comp_dir) {
    if (is.null(OrgDb_species)) {
        message("  Skipping GO enrichment: no OrgDb loaded.")
        return(invisible(NULL))
    }
    tryCatch({
        res_df <- as.data.frame(res)
        sig_genes <- rownames(res_df)[!is.na(res_df$padj) &
            res_df$padj < padj_threshold &
            abs(res_df$log2FoldChange) > lfc_threshold]

        if (length(sig_genes) < 5) {
            message("  Skipping GO enrichment: fewer than 5 significant genes.")
            return(invisible(NULL))
        }

        universe <- rownames(res_df)

        for (ont in c("BP", "MF", "CC")) {
            cat(sprintf("  Running GO %s enrichment...\n", ont))
            ego <- enrichGO(gene = sig_genes,
                            universe = universe,
                            OrgDb = OrgDb_species,
                            keyType = "ENSEMBL",
                            ont = ont,
                            pAdjustMethod = "BH",
                            pvalueCutoff = 0.05,
                            qvalueCutoff = 0.2,
                            readable = TRUE)

            if (!is.null(ego) && nrow(ego) > 0) {
                write.csv(as.data.frame(ego),
                          file = file.path(comp_dir, paste0("GO_", ont, "_enrichment.csv")),
                          row.names = FALSE)
                tryCatch({
                    p <- dotplot(ego, showCategory = min(20, nrow(ego)),
                                 title = paste0("GO ", ont, ": ", comp_name))
                    ggsave(file.path(comp_dir, paste0("GO_", ont, "_dotplot.png")),
                           p, width = 10, height = 8, dpi = 300)
                    cat(sprintf("  GO %s dotplot saved.\n", ont))
                }, error = function(e) {
                    message(sprintf("  Warning: Could not create GO %s dotplot: %s", ont, e$message))
                })
            } else {
                message(sprintf("  No significant GO %s terms found.", ont))
            }
        }
    }, error = function(e) { message("  ERROR in GO enrichment: ", e$message) })
}

# =============================================================================
# DESeq2 Analysis
# =============================================================================
cat("\n--- Running DESeq2 ---\n")

count_matrix_int <- round(count_matrix)
storage.mode(count_matrix_int) <- "integer"

dds <- DESeqDataSetFromMatrix(countData = count_matrix_int, colData = metadata, design = ~ condition)
cat("DESeqDataSet created.\n")

cat("Running DESeq()...\n")
dds <- DESeq(dds)
cat("DESeq() completed.\n")
cat("Available results names:", paste(resultsNames(dds), collapse = ", "), "\n")

# =============================================================================
# Load or Generate Comparisons
# =============================================================================
all_conditions <- levels(dds$condition)
comparisons <- list()

if (!is.null(comparisons_file) && file.exists(comparisons_file) && file.info(comparisons_file)$size > 0) {
    cat("Reading comparisons from:", comparisons_file, "\n")
    comp_raw <- tryCatch({
        read.delim(comparisons_file, header = FALSE, stringsAsFactors = FALSE, sep = "\t")
    }, error = function(e) {
        # Fallback: try space-separated
        read.delim(comparisons_file, header = FALSE, stringsAsFactors = FALSE, sep = " ")
    })
    if (ncol(comp_raw) == 1) {
        # If tab didn't split, try space
        comp_raw <- read.delim(comparisons_file, header = FALSE, stringsAsFactors = FALSE, sep = " ")
    }
    for (i in 1:nrow(comp_raw)) {
        comparisons[[i]] <- c(trimws(comp_raw[i, 1]), trimws(comp_raw[i, 2]))
    }
} else {
    cat("No comparisons file found. Generating all pairwise comparisons.\n")
    for (i in 1:(length(all_conditions) - 1)) {
        for (j in (i + 1):length(all_conditions)) {
            comparisons[[length(comparisons) + 1]] <- c(all_conditions[i], all_conditions[j])
        }
    }
}
cat("Comparisons to perform:", length(comparisons), "\n")
for (i in seq_along(comparisons)) {
    cat(sprintf("  %d. %s vs %s\n", i, comparisons[[i]][1], comparisons[[i]][2]))
}

# =============================================================================
# Run Comparisons
# =============================================================================
cat("\n--- Performing Differential Expression Comparisons ---\n")

comparison_results <- list()

for (comp in comparisons) {
    condition1 <- comp[1]
    condition2 <- comp[2]

    # Determine comparison name
    # Column 1 = Treatment (Numerator)
    # Column 2 = Reference (Denominator)
    
    treatment <- condition1
    reference <- condition2
    comparison_name <- paste(treatment, "vs", reference, sep = "_")
    comparison_dir <- file.path(outdir, comparison_name)
    dir.create(comparison_dir, showWarnings = FALSE, recursive = TRUE)
    cat(sprintf("\n===== %s =====\n", comparison_name))

    # Get results: positive LFC = upregulated in treatment
    res <- NULL
    tryCatch({
        res <- results(dds, contrast = c("condition", treatment, reference), alpha = padj_threshold)
    }, error = function(e) {
        message(sprintf("  ERROR getting results for '%s': %s", comparison_name, e$message))
    })

    if (is.null(res)) {
        message(sprintf("  Skipping comparison '%s' due to results error.", comparison_name))
        next
    }

    res <- res[order(res$padj), ]
    comparison_results[[comparison_name]] <- res

    # --- Export all genes ---
    res_df <- as.data.frame(res)
    res_df$ensembl_id <- rownames(res_df)
    res_df$symbol <- get_gene_symbols(res_df$ensembl_id)
    res_df <- res_df[, c("ensembl_id", "symbol", setdiff(colnames(res_df), c("ensembl_id", "symbol")))]
    write.csv(res_df, file = file.path(comparison_dir, "all_genes.csv"), row.names = FALSE)

    # --- Export significant genes ---
    sig_df <- subset(res_df, !is.na(padj) & padj < padj_threshold & abs(log2FoldChange) > lfc_threshold)
    write.csv(sig_df, file = file.path(comparison_dir, "significant_genes.csv"), row.names = FALSE)

    n_sig  <- nrow(sig_df)
    n_up   <- sum(sig_df$log2FoldChange > 0, na.rm = TRUE)
    n_down <- sum(sig_df$log2FoldChange < 0, na.rm = TRUE)
    cat(sprintf("  Total genes tested: %d\n", nrow(res_df)))
    cat(sprintf("  Significant (padj < %s, |LFC| > %s): %d (Up: %d, Down: %d)\n",
                padj_threshold, lfc_threshold, n_sig, n_up, n_down))

    # --- MA Plot ---
    create_ma_plot(res, sprintf("MA Plot: %s", comparison_name),
                   file.path(comparison_dir, "ma_plot.png"))

    # --- Volcano Plot ---
    create_volcano_plot(res, sprintf("Volcano Plot: %s", comparison_name),
                        file.path(comparison_dir, "volcano_plot.png"))

    # --- Heatmap of Top DEGs ---
    if (n_sig > 0) {
        num_heat <- min(50, n_sig)
        top_genes <- sig_df$ensembl_id[1:num_heat]

        # Subset to just the two conditions and drop unused factor levels
        dds_subset <- dds[, dds$condition %in% c(treatment, reference)]
        dds_subset$condition <- droplevels(dds_subset$condition)

        if (ncol(dds_subset) > 1 && nlevels(dds_subset$condition) > 1) {
            tryCatch({
                vst_counts <- assay(vst(dds_subset, blind = FALSE))
                heatmap_mat <- vst_counts[rownames(vst_counts) %in% top_genes, , drop = FALSE]
                if (nrow(heatmap_mat) > 0) {
                    create_heatmap(heatmap_mat, colData(dds_subset),
                                   sprintf("Top %d DEGs: %s (VST)", nrow(heatmap_mat), comparison_name),
                                   file.path(comparison_dir, "top_degs_heatmap.png"))
                }
            }, error = function(e) { message("  Error generating heatmap: ", e$message) })
        }
    } else {
        cat("  No significant genes found for heatmap.\n")
    }

    # --- GO Enrichment ---
    perform_GO_enrichment(res, comparison_name, comparison_dir)

    cat(sprintf("===== Finished: %s =====\n", comparison_name))
}

# =============================================================================
# PCA Plot (all samples)
# =============================================================================
cat("\n--- Generating PCA Plot ---\n")
tryCatch({
    vst_data <- vst(dds, blind = FALSE)
    pcaData <- plotPCA(vst_data, intgroup = c("condition", "sample_name"), returnData = TRUE)
    percentVar <- round(100 * attr(pcaData, "percentVar"))

    pca_plot <- ggplot(pcaData, aes(x = PC1, y = PC2, color = condition, label = sample_name)) +
        geom_point(size = 3) +
        ggrepel::geom_text_repel(size = 3, box.padding = 0.5, point.padding = 0.5,
                                  max.overlaps = Inf, show.legend = FALSE,
                                  segment.color = "grey50", min.segment.length = 0) +
        xlab(paste0("PC1: ", percentVar[1], "% variance")) +
        ylab(paste0("PC2: ", percentVar[2], "% variance")) +
        theme_minimal(base_size = 12) +
        theme(legend.position = "right",
              panel.grid.major = element_line(colour = "grey90"),
              panel.grid.minor = element_blank()) +
        ggtitle("PCA of RNA-seq Samples (VST transformed)") +
        scale_color_brewer(palette = "Set1")
    ggsave(file.path(outdir, "pca_plot.png"), pca_plot, width = 10, height = 7, dpi = 300)
    cat("PCA plot saved.\n")
}, error = function(e) { message("ERROR creating PCA plot: ", e$message) })

# =============================================================================
# Summary Table
# =============================================================================
cat("\n--- Creating Summary Table ---\n")
summary_df <- data.frame(
    Comparison = character(),
    Total_Genes = integer(),
    Significant_Genes = integer(),
    Up_regulated = integer(),
    Down_regulated = integer(),
    stringsAsFactors = FALSE
)

if (length(comparison_results) > 0) {
    for (comp_name in names(comparison_results)) {
        res <- comparison_results[[comp_name]]
        sig   <- sum(!is.na(res$padj) & res$padj < padj_threshold & abs(res$log2FoldChange) > lfc_threshold)
        up    <- sum(!is.na(res$padj) & res$padj < padj_threshold & res$log2FoldChange > lfc_threshold)
        down  <- sum(!is.na(res$padj) & res$padj < padj_threshold & res$log2FoldChange < -lfc_threshold)
        summary_df <- rbind(summary_df, data.frame(
            Comparison = comp_name,
            Total_Genes = nrow(res),
            Significant_Genes = sig,
            Up_regulated = up,
            Down_regulated = down
        ))
    }
}

write.csv(summary_df, file = file.path(outdir, "deseq2_summary.csv"), row.names = FALSE)
cat("Summary table saved.\n")
print(summary_df)

# =============================================================================
# Wrap Up
# =============================================================================
end_time <- Sys.time()
cat("\n=== DESeq2 Analysis Complete ===\n")
cat("End time:", format(end_time), "\n")
cat("Total execution time:", format(difftime(end_time, start_time)), "\n")

# Explicitly exit
quit(save = "no", status = 0)