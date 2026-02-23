import pandas as pd
import re, os, sys

configfile: "config/config.yaml"

# --- Setup & Configuration ---
EXP_ID = config.get("id")
if not EXP_ID: sys.exit("CRITICAL ERROR: No Experiment ID provided. Run with ./run_pipeline.sh -e EXP_ID")

SAMPLE_SHEET = f"config/samples_{EXP_ID}.tsv"
OUT_DIR = f"results_{EXP_ID}"
TMP_DIR = f"{OUT_DIR}/tmp"
BENCH_DIR = f"{OUT_DIR}/benchmarks"

# --- DESeq2 Configuration ---
RUN_DESEQ2 = str(config.get("run_deseq2", "False")).lower() in ["true", "yes", "1"]
COMPARISONS_FILE = f"config/comparisons_{EXP_ID}.tsv"
METADATA_FILE = f"config/metadata_{EXP_ID}.tsv"

# Only load/validate the sample sheet if we are NOT in setup mode
if EXP_ID == "SETUP":
    samples = pd.DataFrame()
    SAMPLES_PE = []
    SAMPLES_SE = []
else:
    if not os.path.exists(SAMPLE_SHEET): sys.exit(f"Error: Sample sheet {SAMPLE_SHEET} not found!")
    samples = pd.read_csv(SAMPLE_SHEET, sep="\t").set_index("sample_name", drop=False)

    SAMPLES_PE = [s for s in samples.index if (not pd.isna(samples.loc[s, "fq2"]) and samples.loc[s, "fq2"] != "NA")]
    SAMPLES_SE = [s for s in samples.index if s not in SAMPLES_PE]


# --- Helper Functions ---

# DYNAMIC RUNTIME CALCULATOR
# usage: runtime=lambda wc, attempt: get_runtime(base_minutes, attempt)
def get_runtime(base, attempt):
    # Logic: base time * (1, 3, 6). 
    # e.g. base 60 -> 60m, 180m, 360m
    multipliers = [1, 3, 6]
    return base * multipliers[min(attempt, 3) - 1]

def get_raw_fastqs(wc):
    if wc.sample in SAMPLES_PE:
        return [samples.loc[wc.sample, "fq1"], samples.loc[wc.sample, "fq2"]]
    return [samples.loc[wc.sample, "fq1"]]

def get_deseq2_targets():
    """Returns DESeq2 output files if enabled, or empty list."""
    if not RUN_DESEQ2:
        return []
    return [f"{OUT_DIR}/deseq2/deseq2_summary.csv", f"{OUT_DIR}/deseq2/pca_plot.png"]

def get_bed_ref(wc):
    g = samples.loc[wc.sample, "genome"]
    return f"references/{g}.bed12"

rule all:
    input:
        expand(f"{OUT_DIR}/bigwig/{{sample}}_fwd.bw", sample=samples.index),
        expand(f"{OUT_DIR}/bigwig/{{sample}}_rev.bw", sample=samples.index),
        f"{OUT_DIR}/counts/gene_counts_matrix.tsv",
        f"{OUT_DIR}/qc/multiqc_report.html",
        expand(f"{OUT_DIR}/qc/fastqc/{{sample}}_1_fastqc.html", sample=SAMPLES_PE),
        expand(f"{OUT_DIR}/qc/fastqc/{{sample}}_2_fastqc.html", sample=SAMPLES_PE),
        expand(f"{OUT_DIR}/qc/fastqc/{{sample}}_fastqc.html", sample=SAMPLES_SE),
        *get_deseq2_targets()
        
#========================================================================================================================================================
# 0. FastQC: Raw data quality control
#========================================================================================================================================================
rule fastqc_pe:
    input: get_raw_fastqs
    output:
        h1 = f"{OUT_DIR}/qc/fastqc/{{sample}}_1_fastqc.html",
        z1 = f"{OUT_DIR}/qc/fastqc/{{sample}}_1_fastqc.zip",
        h2 = f"{OUT_DIR}/qc/fastqc/{{sample}}_2_fastqc.html",
        z2 = f"{OUT_DIR}/qc/fastqc/{{sample}}_2_fastqc.zip"
    params: 
        outdir = f"{OUT_DIR}/qc/fastqc"
    wildcard_constraints: sample="|".join(re.escape(s) for s in SAMPLES_PE) if SAMPLES_PE else "NO_PE"
    threads: 4
    resources: mem_mb=8000, runtime=lambda wc, attempt: get_runtime(30, attempt)
    container: "docker://quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0"
    shell:
        """
        mkdir -p {params.outdir}
        
        # 1. Create temporary symlinks with predictable names so FastQC outputs match Snakemake expectations
        ln -sf $(readlink -f {input[0]}) {params.outdir}/{wildcards.sample}_1.fastq.gz
        ln -sf $(readlink -f {input[1]}) {params.outdir}/{wildcards.sample}_2.fastq.gz
        
        # 2. Run FastQC on the symlinks
        fastqc -t {threads} --outdir {params.outdir} \
            {params.outdir}/{wildcards.sample}_1.fastq.gz \
            {params.outdir}/{wildcards.sample}_2.fastq.gz
            
        # 3. Cleanup links (Outputs remain)
        rm {params.outdir}/{wildcards.sample}_1.fastq.gz {params.outdir}/{wildcards.sample}_2.fastq.gz
        """

rule fastqc_se:
    input: get_raw_fastqs
    output:
        h1 = f"{OUT_DIR}/qc/fastqc/{{sample}}_fastqc.html",
        z1 = f"{OUT_DIR}/qc/fastqc/{{sample}}_fastqc.zip"
    params: 
        outdir = f"{OUT_DIR}/qc/fastqc"
    wildcard_constraints: sample="|".join(re.escape(s) for s in SAMPLES_SE) if SAMPLES_SE else "NO_SE"
    threads: 2
    resources: mem_mb=8000, runtime=lambda wc, attempt: get_runtime(30, attempt)
    container: "docker://quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0"
    shell:
        """
        mkdir -p {params.outdir}
        
        # 1. Create temporary symlink
        ln -sf $(readlink -f {input[0]}) {params.outdir}/{wildcards.sample}.fastq.gz
        
        # 2. Run FastQC
        fastqc -t {threads} --outdir {params.outdir} \
            {params.outdir}/{wildcards.sample}.fastq.gz
            
        # 3. Cleanup
        rm {params.outdir}/{wildcards.sample}.fastq.gz
        """        
        
#========================================================================================================================================================
# 1. Trimming: Cleans up adapters and low-quality bases from reads.
#========================================================================================================================================================
rule trim_pe:
    input: get_raw_fastqs
    output:
        r1=f"{OUT_DIR}/trimmed/{{sample}}_val_1.fq.gz",
        r2=f"{OUT_DIR}/trimmed/{{sample}}_val_2.fq.gz"
    log: f"{OUT_DIR}/logs/trim/{{sample}}.log"
    benchmark: f"{BENCH_DIR}/trim_pe/{{sample}}.tsv"
    wildcard_constraints: sample="|".join(re.escape(s) for s in SAMPLES_PE) if SAMPLES_PE else "NO_PE"
    threads: 8
    resources: mem_mb=4000, runtime=lambda wc, attempt: get_runtime(40, attempt)
    container: "docker://quay.io/biocontainers/trim-galore:0.6.10--hdfd78af_0"
    shell: "trim_galore --paired --cores {threads} --gzip -o " + OUT_DIR + "/trimmed --basename {wildcards.sample} {input} > {log} 2>&1"

rule trim_se:
    input: get_raw_fastqs
    output: r1=f"{OUT_DIR}/trimmed/{{sample}}_trimmed.fq.gz"
    log: f"{OUT_DIR}/logs/trim/{{sample}}.log"
    benchmark: f"{BENCH_DIR}/trim_se/{{sample}}.tsv"
    wildcard_constraints: sample="|".join(re.escape(s) for s in SAMPLES_SE) if SAMPLES_SE else "NO_SE"
    threads: 8
    resources: mem_mb=4000, runtime=lambda wc, attempt: get_runtime(40, attempt)
    container: "docker://quay.io/biocontainers/trim-galore:0.6.10--hdfd78af_0"
    shell: "trim_galore --cores {threads} --gzip -o " + OUT_DIR + "/trimmed --basename {wildcards.sample} {input} > {log} 2>&1"

#========================================================================================================================================================
# 2. Mapping: Aligns reads to the reference genome using STAR (splice-aware, RNA-seq parameters).
#========================================================================================================================================================
rule star_pe:
    input:
        r1=f"{OUT_DIR}/trimmed/{{sample}}_val_1.fq.gz",
        r2=f"{OUT_DIR}/trimmed/{{sample}}_val_2.fq.gz"
    output:
        bam=f"{OUT_DIR}/mapped/{{sample}}_Aligned.sortedByCoord.out.bam",
        log_final=f"{OUT_DIR}/mapped/{{sample}}_Log.final.out"
    log: f"{OUT_DIR}/logs/star/{{sample}}.log"
    benchmark: f"{BENCH_DIR}/star_pe/{{sample}}.tsv"
    wildcard_constraints: sample="|".join(re.escape(s) for s in SAMPLES_PE) if SAMPLES_PE else "NO_PE"
    params:
        idx=lambda wc: config["refs"][samples.loc[wc.sample, "genome"]]["star_index"],
        gtf=lambda wc: config["refs"][samples.loc[wc.sample, "genome"]]["gtf"],
        m=config["star"]["multimap_max"],
        p=f"{OUT_DIR}/mapped/{{sample}}_",
        t=f"{TMP_DIR}/{{sample}}_STARtmp"
    threads: 16
    resources: mem_mb=40000, disk_mb=40000, runtime=lambda wc, attempt: get_runtime(60, attempt)
    container: "docker://quay.io/biocontainers/star:2.7.4a--0"
    shell:
        r"""
        mkdir -p {TMP_DIR}
        rm -rf {params.t}
        STAR --runMode alignReads \
            --genomeDir {params.idx} \
            --sjdbGTFfile {params.gtf} \
            --readFilesIn {input.r1} {input.r2} \
            --readFilesCommand zcat \
            --runThreadN {threads} \
            --outSAMtype BAM SortedByCoordinate \
            --outFilterType BySJout \
            --outFilterMultimapNmax {params.m} \
            --outFilterMismatchNoverLmax 0.05 \
            --outFilterIntronMotifs RemoveNoncanonical \
            --outSAMmultNmax 1 \
            --outMultimapperOrder Random \
            --outSAMattributes NH HI NM MD AS XS \
            --outSAMunmapped Within \
            --outBAMsortingThreadN {threads} \
            --alignSJoverhangMin 8 \
            --alignSJDBoverhangMin 1 \
            --alignIntronMin 20 \
            --alignIntronMax 1000000 \
            --alignMatesGapMax 1000000 \
            --sjdbScore 1 \
            --outFileNamePrefix {params.p} \
            --outTmpDir {params.t} > {log} 2>&1
        mv {params.p}Aligned.sortedByCoord.out.bam {output.bam}
        rm -rf {params.t}
        """

rule star_se:
    input: r1=f"{OUT_DIR}/trimmed/{{sample}}_trimmed.fq.gz"
    output:
        bam=f"{OUT_DIR}/mapped/{{sample}}_Aligned.sortedByCoord.out.bam",
        log_final=f"{OUT_DIR}/mapped/{{sample}}_Log.final.out"
    log: f"{OUT_DIR}/logs/star/{{sample}}.log"
    benchmark: f"{BENCH_DIR}/star_se/{{sample}}.tsv"
    wildcard_constraints: sample="|".join(re.escape(s) for s in SAMPLES_SE) if SAMPLES_SE else "NO_SE"
    params:
        idx=lambda wc: config["refs"][samples.loc[wc.sample, "genome"]]["star_index"],
        gtf=lambda wc: config["refs"][samples.loc[wc.sample, "genome"]]["gtf"],
        m=config["star"]["multimap_max"],
        p=f"{OUT_DIR}/mapped/{{sample}}_",
        t=f"{TMP_DIR}/{{sample}}_STARtmp"
    threads: 16
    resources: mem_mb=40000, disk_mb=40000, runtime=lambda wc, attempt: get_runtime(60, attempt)
    container: "docker://quay.io/biocontainers/star:2.7.4a--0"
    shell:
        r"""
        mkdir -p {TMP_DIR}
        rm -rf {params.t}
        STAR --runMode alignReads \
            --genomeDir {params.idx} \
            --sjdbGTFfile {params.gtf} \
            --readFilesIn {input.r1} \
            --readFilesCommand zcat \
            --runThreadN {threads} \
            --outSAMtype BAM SortedByCoordinate \
            --outFilterType BySJout \
            --outFilterMultimapNmax {params.m} \
            --outFilterMismatchNoverLmax 0.05 \
            --outFilterIntronMotifs RemoveNoncanonical \
            --outSAMmultNmax 1 \
            --outMultimapperOrder Random \
            --outSAMattributes NH HI NM MD AS XS \
            --outSAMunmapped Within \
            --outBAMsortingThreadN {threads} \
            --alignSJoverhangMin 8 \
            --alignSJDBoverhangMin 1 \
            --alignIntronMin 20 \
            --alignIntronMax 1000000 \
            --alignMatesGapMax 1000000 \
            --sjdbScore 1 \
            --outFileNamePrefix {params.p} \
            --outTmpDir {params.t} > {log} 2>&1
        mv {params.p}Aligned.sortedByCoord.out.bam {output.bam}
        rm -rf {params.t}
        """

rule samtools_index:
    input: f"{OUT_DIR}/mapped/{{sample}}_Aligned.sortedByCoord.out.bam"
    output: f"{OUT_DIR}/mapped/{{sample}}_Aligned.sortedByCoord.out.bam.bai"
    resources: mem_mb=8000, runtime=lambda wc, attempt: get_runtime(15, attempt)
    container: "docker://quay.io/biocontainers/samtools:1.19--h50ea8bc_0"
    shell: "samtools index {input}"

#========================================================================================================================================================
# 3. Strand-specific BigWigs: Forward and Reverse strand BigWigs using deepTools bamCoverage.
#    Uses --filterRNAstrand for strand separation. CPM normalized, binSize 10.
#========================================================================================================================================================
rule bigwig:
    input:
        bam=f"{OUT_DIR}/mapped/{{sample}}_Aligned.sortedByCoord.out.bam",
        bai=f"{OUT_DIR}/mapped/{{sample}}_Aligned.sortedByCoord.out.bam.bai"
    output:
        fwd=f"{OUT_DIR}/bigwig/{{sample}}_fwd.bw",
        rev=f"{OUT_DIR}/bigwig/{{sample}}_rev.bw"
    log: f"{OUT_DIR}/logs/bigwig/{{sample}}.log"
    benchmark: f"{BENCH_DIR}/bigwig/{{sample}}.tsv"
    threads: 4
    resources: mem_mb=6000, runtime=lambda wc, attempt: get_runtime(60, attempt)
    container: "docker://quay.io/biocontainers/deeptools:3.5.6--pyhdfd78af_0"
    shell:
        """
        mkdir -p {TMP_DIR}

        export TMPDIR={TMP_DIR}
        export TEMP={TMP_DIR}
        export TMP={TMP_DIR}
        export MPLCONFIGDIR={TMP_DIR}

        mkdir -p {OUT_DIR}/bigwig

        # Forward strand (sense)
        bamCoverage --bam {input.bam} --outFileName {output.fwd} \
            --outFileFormat bigwig \
            --binSize 10 \
            --normalizeUsing CPM \
            --filterRNAstrand forward \
            --numberOfProcessors {threads} > {log} 2>&1

        # Reverse strand (antisense)
        bamCoverage --bam {input.bam} --outFileName {output.rev} \
            --outFileFormat bigwig \
            --binSize 10 \
            --normalizeUsing CPM \
            --filterRNAstrand reverse \
            --numberOfProcessors {threads} >> {log} 2>&1
        """

#========================================================================================================================================================
# 4. RSeQC quality control 
#========================================================================================================================================================

# --- RSeQC Rules ---
# 1. Convert GTF to BED12 (Required for RSeQC)
#    Replaces bedops gtf2bed to avoid formatting errors with Gencode GTFs.
rule gtf_to_bed:
    input: lambda wc: config["refs"][wc.genome]["gtf"]
    output: "references/{genome}.bed12" 
    log: f"{OUT_DIR}/logs/rseqc/gtf2bed_{{genome}}.log" 
    resources: mem_mb=8000, runtime=lambda wc, attempt: get_runtime(30, attempt)
    run:
        import sys
        import csv
        import gzip
        
        # Function to open file (gz or plain)
        def open_file(fname):
            if fname.endswith(".gz"):
                return gzip.open(fname, "rt")
            return open(fname, "r")

        sys.stderr = open(log[0], "w")
        print(f"Converting GTF: {input[0]} to BED12: {output[0]}", file=sys.stderr)

        transcripts = {}

        # 1. Parse GTF
        # We group exons by transcript_id to build the BED12 structure
        with open_file(input[0]) as f:
            reader = csv.reader(f, delimiter="\t", quoting=csv.QUOTE_NONE)
            for row in reader:
                if len(row) < 9 or row[0].startswith("#"):
                    continue
                
                feature = row[2]
                # We only care about exons and CDS to build the model
                if feature not in ["exon", "CDS"]:
                    continue

                # Parse attributes to find transcript_id
                # Attributes look like: gene_id "ENSG00.."; transcript_id "ENST00..";
                attr_str = row[8]
                try:
                    tid = attr_str.split('transcript_id "')[1].split('"')[0]
                except IndexError:
                    continue # Skip lines without transcript_id

                if tid not in transcripts:
                    transcripts[tid] = {
                        "chrom": row[0],
                        "strand": row[6],
                        "exons": [], # list of (start, end)
                        "cds": []    # list of (start, end)
                    }

                start = int(row[3]) - 1 # 0-based conversion
                end = int(row[4])

                if feature == "exon":
                    transcripts[tid]["exons"].append((start, end))
                elif feature == "CDS":
                    transcripts[tid]["cds"].append((start, end))

        # 2. Write BED12
        with open(output[0], "w") as out:
            for tid, data in transcripts.items():
                chrom = data["chrom"]
                strand = data["strand"]
                
                # Sort exons by start
                sorted_exons = sorted(data["exons"], key=lambda x: x[0])
                if not sorted_exons:
                    continue

                tx_start = sorted_exons[0][0]
                tx_end = sorted_exons[-1][1]

                # Determine CDS (thick) start/end
                # If no CDS (non-coding), thickStart = thickEnd = tx_start
                if data["cds"]:
                    sorted_cds = sorted(data["cds"], key=lambda x: x[0])
                    thick_start = sorted_cds[0][0]
                    thick_end = sorted_cds[-1][1]
                else:
                    thick_start = tx_start
                    thick_end = tx_start

                # Block calculations
                block_count = len(sorted_exons)
                block_sizes = []
                block_starts = []

                for (estart, eend) in sorted_exons:
                    block_sizes.append(str(eend - estart))
                    block_starts.append(str(estart - tx_start))

                # BED12 columns:
                # 1. chrom, 2. start, 3. end, 4. name, 5. score, 6. strand, 
                # 7. thickStart, 8. thickEnd, 9. itemRgb, 10. blockCount, 
                # 11. blockSizes, 12. blockStarts
                bed_row = [
                    chrom, str(tx_start), str(tx_end), tid, "0", strand,
                    str(thick_start), str(thick_end), "0", str(block_count),
                    ",".join(block_sizes), ",".join(block_starts)
                ]
                out.write("\t".join(bed_row) + "\n")
        
        print(f"Conversion complete. Processed {len(transcripts)} transcripts.", file=sys.stderr)



# 2. Check duplicates

rule read_duplication:
    input:
        bam = f"{OUT_DIR}/mapped/{{sample}}_Aligned.sortedByCoord.out.bam",
        bai = f"{OUT_DIR}/mapped/{{sample}}_Aligned.sortedByCoord.out.bam.bai"
    output:
        pos  = f"{OUT_DIR}/qc/dedup/{{sample}}.pos.DupRate.xls",
        seq  = f"{OUT_DIR}/qc/dedup/{{sample}}.seq.DupRate.xls",
        plot = f"{OUT_DIR}/qc/dedup/{{sample}}.DupRate_plot.r"
    log: f"{OUT_DIR}/logs/dedup/{{sample}}.log"
    benchmark: f"{BENCH_DIR}/dedup/{{sample}}.tsv"
    resources: mem_mb=8000, runtime=lambda wc, attempt: get_runtime(50, attempt)
    threads: 2
    container: "docker://quay.io/biocontainers/rseqc:5.0.4--pyhdfd78af_1"
    shell:
        """
        read_duplication.py \
            -i {input.bam} \
            -o {OUT_DIR}/qc/dedup/{wildcards.sample} \
            > {log} 2>&1
        """

# 3. Infer Strandedness using RSeQC
rule infer_strandedness:
    input:
        bam = f"{OUT_DIR}/mapped/{{sample}}_Aligned.sortedByCoord.out.bam",
        bai = f"{OUT_DIR}/mapped/{{sample}}_Aligned.sortedByCoord.out.bam.bai",
        bed = get_bed_ref
    output:
        txt = f"{OUT_DIR}/qc/rseqc/{{sample}}_infer_experiment.txt",
        flag = f"{OUT_DIR}/qc/rseqc/{{sample}}_strandedness.flag"
    log: f"{OUT_DIR}/logs/rseqc/infer_{{sample}}.log"
    resources: mem_mb=4000, runtime=lambda wc, attempt: get_runtime(15, attempt)
    container: "docker://quay.io/biocontainers/rseqc:5.0.4--pyhdfd78af_1"
    shell:
        """
        # 1. Run RSeQC inference (sample 200k reads)
        infer_experiment.py -r {input.bed} -i {input.bam} -s 200000 > {output.txt} 2> {log}

        # 2. Parse output to determine featureCounts flag (0, 1, or 2)
        # We use a small python one-liner inside the shell to handle the parsing logic
        
        python3 -c '
import sys
fwd = 0.0
rev = 0.0
with open("{output.txt}") as f:
    for line in f:
        # Pattern for Forward signal (1++,1-- or ++,--)
        if "1++,1--" in line or "++,--" in line:
            parts = line.strip().split(":")
            if len(parts) > 1: fwd = float(parts[-1])
            
        # Pattern for Reverse signal (1+-,1-+ or +-,-+)
        if "1+-,1-+" in line or "+-,-+" in line:
            parts = line.strip().split(":")
            if len(parts) > 1: rev = float(parts[-1])

# Logic: If > 0.8 signal, valid strand. Else unstranded (0).
# 2 = Reverse (Standard Illumina TruSeq)
# 1 = Forward
# 0 = Unstranded
res = "0"
if fwd > 0.8: res = "1"
elif rev > 0.8: res = "2"

with open("{output.flag}", "w") as out:
    out.write(res)
'
        """


# 4. Read Distribution (Critical for checking gDNA contamination/mRNA enrichment)
rule rseqc_read_dist:
    input:
        bam = f"{OUT_DIR}/mapped/{{sample}}_Aligned.sortedByCoord.out.bam",
        bai = f"{OUT_DIR}/mapped/{{sample}}_Aligned.sortedByCoord.out.bam.bai",
        bed = get_bed_ref
    output: f"{OUT_DIR}/qc/rseqc/{{sample}}_read_distribution.txt"
    log: f"{OUT_DIR}/logs/rseqc/read_dist_{{sample}}.log"
    resources: mem_mb=8000, runtime=lambda wc, attempt: get_runtime(40, attempt)
    container: "docker://quay.io/biocontainers/rseqc:5.0.4--pyhdfd78af_1"
    shell:
        "read_distribution.py -i {input.bam} -r {input.bed} > {output} 2> {log}"

# 5. Junction Annotation (Checks splicing efficiency (Optional))
#rule rseqc_junction_annotation:
#    input:
#        bam = f"{OUT_DIR}/mapped/{{sample}}_Aligned.sortedByCoord.out.bam",
#        bai = f"{OUT_DIR}/mapped/{{sample}}_Aligned.sortedByCoord.out.bam.bai",
#        bed = get_bed_ref
#    output: f"{OUT_DIR}/qc/rseqc/{{sample}}_junction_annotation.txt"
#    log: f"{OUT_DIR}/logs/rseqc/junction_{{sample}}.log"
#    params: prefix=f"{OUT_DIR}/qc/rseqc/{{sample}}_junction"
#    resources: mem_mb=8000, runtime=lambda wc, attempt: get_runtime(40, attempt)
#    container: "docker://quay.io/biocontainers/rseqc:5.0.4--pyhdfd78af_1"
#    shell:
#        """
#        # junction_annotation produces multiple files/plots, we dump the stats to stdout/txt
#        junction_annotation.py -i {input.bam} -r {input.bed} -o {params.prefix} > {output} 2> {log}
#        """
        
#========================================================================================================================================================
# 5. featureCounts: Counts reads per gene using the GTF annotation (part of the Subread package).
#========================================================================================================================================================
rule featurecounts:
    input:
        bam=f"{OUT_DIR}/mapped/{{sample}}_Aligned.sortedByCoord.out.bam",
        bai=f"{OUT_DIR}/mapped/{{sample}}_Aligned.sortedByCoord.out.bam.bai",
        strand_flag=f"{OUT_DIR}/qc/rseqc/{{sample}}_strandedness.flag"
    output:
        counts=f"{OUT_DIR}/counts/{{sample}}_counts.txt",
        summary=f"{OUT_DIR}/counts/{{sample}}_counts.txt.summary"
    log: f"{OUT_DIR}/logs/featurecounts/{{sample}}.log"
    benchmark: f"{BENCH_DIR}/featurecounts/{{sample}}.tsv"
    params:
        gtf=lambda wc: config["refs"][samples.loc[wc.sample, "genome"]]["gtf"],
        pe_flags=lambda wc: "-p --countReadPairs -B" if wc.sample in SAMPLES_PE else ""
    threads: 4
    resources: mem_mb=8000, runtime=lambda wc, attempt: get_runtime(30, attempt)
    container: "docker://quay.io/biocontainers/subread:2.1.1--h577a1d6_0"
    shell:
        r"""
        # Read the auto-detected strand flag
        STRAND=$(cat {input.strand_flag})
        
        echo "Running featureCounts with auto-detected strandedness: $STRAND" > {log}
        
        featureCounts \
            -a {params.gtf} \
            -o {output.counts} \
            -T {threads} \
            -s $STRAND \
            -t exon \
            -g gene_id \
            --minOverlap 1 \
            --fracOverlap 0 \
            --primary \
            {params.pe_flags} \
            {input.bam} >> {log} 2>&1
        """
#========================================================================================================================================================
# 6. Merge Counts: Combines per-sample featureCounts output into a single gene count matrix.
#========================================================================================================================================================
rule merge_counts:
    input: expand(f"{OUT_DIR}/counts/{{sample}}_counts.txt", sample=samples.index)
    output: f"{OUT_DIR}/counts/gene_counts_matrix.tsv"
    log: f"{OUT_DIR}/logs/merge_counts.log"
    resources: mem_mb=4000, runtime=lambda wc, attempt: get_runtime(5, attempt)
    run:
        import pandas as pd

        merged = None
        for f in input:
            sample_name = os.path.basename(f).replace("_counts.txt", "")
            df = pd.read_csv(f, sep="\t", comment="#")
            count_col = df.columns[-1]
            df = df.rename(columns={count_col: sample_name})

            if merged is None:
                merged = df[["Geneid", "Length", sample_name]]
            else:
                merged = merged.merge(df[["Geneid", sample_name]], on="Geneid", how="outer")

        merged = merged.fillna(0)
        for col in merged.columns:
            if col not in ["Geneid", "Length"]:
                merged[col] = merged[col].astype(int)

        merged.to_csv(str(output), sep="\t", index=False)
        with open(str(log), "w") as logf:
            logf.write(f"Merged {len(input)} samples. Matrix shape: {merged.shape}\n")

#========================================================================================================================================================
# 7. DESeq2: Differential expression analysis with QC plots, volcano/MA plots, heatmaps, and GO enrichment.
#    Only runs if run_deseq2=True is passed via config. Uses custom-built container with all R packages.
#========================================================================================================================================================
rule deseq2:
    input:
        counts=f"{OUT_DIR}/counts/gene_counts_matrix.tsv",
        metadata=METADATA_FILE,
        script="run_deseq2.R"
    output:
        summary=f"{OUT_DIR}/deseq2/deseq2_summary.csv",
        pca=f"{OUT_DIR}/deseq2/pca_plot.png"
    log: f"{OUT_DIR}/logs/deseq2.log"
    benchmark: f"{BENCH_DIR}/deseq2.tsv"
    params:
        species=lambda wc: samples.iloc[0]["genome"],
        outdir=f"{OUT_DIR}/deseq2",
        lfc_threshold=config.get("deseq2", {}).get("lfc_threshold", 1.0),
        padj_threshold=config.get("deseq2", {}).get("padj_threshold", 0.05),
        comparisons_file=COMPARISONS_FILE
    threads: 8
    resources: mem_mb=16000, runtime=lambda wc, attempt: get_runtime(60, attempt)
    container: "docker://gavinkarsters/rnaseq-r_env:v1"
    shell:
        r"""
        export LC_ALL=en_US.UTF-8
        export LANG=en_US.UTF-8

        Rscript {input.script} \
            --counts {input.counts} \
            --metadata {input.metadata} \
            --comparisons {params.comparisons_file} \
            --outdir {params.outdir} \
            --species {params.species} \
            --lfc_threshold {params.lfc_threshold} \
            --padj_threshold {params.padj_threshold} \
            --threads {threads} > {log} 2>&1
        """
#========================================================================================================================================================
# 8. MultiQC: Compiles QC metrics (STAR logs, Trim Galore logs, featureCounts summaries, benchmarks).
#========================================================================================================================================================
rule count_pipeline_stats:
    input:
        b_pe=expand(f"{BENCH_DIR}/{{r}}/{{s}}.tsv", r=["trim_pe","star_pe"], s=SAMPLES_PE),
        b_se=expand(f"{BENCH_DIR}/{{r}}/{{s}}.tsv", r=["trim_se","star_se"], s=SAMPLES_SE),
        b_all=expand(f"{BENCH_DIR}/{{r}}/{{s}}.tsv", r=["bigwig","featurecounts"], s=samples.index)
    output:
        time=f"{OUT_DIR}/qc/execution_time_mqc.tsv",
        mem=f"{OUT_DIR}/qc/memory_usage_mqc.tsv"
    resources: mem_mb=1000, runtime=lambda wc, attempt: get_runtime(10, attempt)
    run:
        import os, pandas as pd

        R_MAP = {
            "trim_pe": "1. Trim", "trim_se": "1. Trim",
            "star_pe": "2. Align", "star_se": "2. Align",
            "bigwig": "3. BigWig", "featurecounts": "4. Counting"
        }
        met = {"Time": {}, "Mem": {}}

        all_bench = [(f, False) for f in input.b_pe + input.b_se + input.b_all]

        for f, _ in all_bench:
            try:
                df = pd.read_csv(f, sep="\t")
                s_name = os.path.basename(f).replace(".tsv", "")
                step = R_MAP.get(f.split("/")[-2], f.split("/")[-2])
                met["Time"].setdefault(s_name, {})[step] = round(df.iloc[0]["s"] / 60, 2)
                met["Mem"].setdefault(s_name, {})[step] = round(df.iloc[0]["max_rss"], 1)
            except:
                pass

        for k, f_out, ylab in [("Time", output.time, "Minutes"), ("Mem", output.mem, "Max RSS (MB)")]:
            if not met[k]:
                open(f_out, 'w').close()
                continue
            with open(f_out, 'w') as fh:
                fh.write(f"# plot_type: 'bargraph'\n# section_name: 'Pipeline {k}'\n# ylab: '{ylab}'\n")
                pd.DataFrame.from_dict(met[k], orient='index').sort_index(axis=1).to_csv(fh, sep="\t")
                
rule multiqc:
    input:
        mqc=[f"{OUT_DIR}/qc/{x}_mqc.tsv" for x in ["execution_time", "memory_usage"]],
        star_logs=expand(f"{OUT_DIR}/mapped/{{s}}_Log.final.out", s=samples.index),
        trim_logs=expand(f"{OUT_DIR}/logs/trim/{{s}}.log", s=samples.index),
        fc_summaries=expand(f"{OUT_DIR}/counts/{{s}}_counts.txt.summary", s=samples.index),
        rseqc_dist=expand(f"{OUT_DIR}/qc/rseqc/{{s}}_read_distribution.txt", s=samples.index),
        #rseqc_junc=expand(f"{OUT_DIR}/qc/rseqc/{{s}}_junction_annotation.txt", s=samples.index),
        rseqc_inf=expand(f"{OUT_DIR}/qc/rseqc/{{s}}_infer_experiment.txt", s=samples.index),
        dedup_pos=expand(f"{OUT_DIR}/qc/dedup/{{s}}.pos.DupRate.xls", s=samples.index),
        dedup_seq=expand(f"{OUT_DIR}/qc/dedup/{{s}}.seq.DupRate.xls", s=samples.index),
        fastqc_zips = expand(f"{OUT_DIR}/qc/fastqc/{{sample}}_1_fastqc.zip", sample=SAMPLES_PE) + \
                      expand(f"{OUT_DIR}/qc/fastqc/{{sample}}_2_fastqc.zip", sample=SAMPLES_PE) + \
                      expand(f"{OUT_DIR}/qc/fastqc/{{sample}}_fastqc.zip", sample=SAMPLES_SE)
    output: f"{OUT_DIR}/qc/multiqc_report.html"
    resources: mem_mb=5000, runtime=lambda wc, attempt: get_runtime(10, attempt)
    container: "docker://quay.io/biocontainers/multiqc:1.21--pyhdfd78af_0"
    shell:
        """
        multiqc {OUT_DIR} -o {OUT_DIR}/qc -n multiqc_report.html -f
        """

       
        
#========================================================================================================================================================
# Genome Setup (Optional / One-time use)
#========================================================================================================================================================
# Usage:
#   snakemake --jobs 1 \
#     --cluster "sbatch --partition=cpu --time={resources.runtime} --mem={resources.mem_mb}M --cpus-per-task={threads} --gres=tmpspace:170G" \
#     --use-singularity \
#     references/star_index_human --config id=SETUP

rule download_genome:
    output:
        fasta="references/{genome}.fa",
        gtf="references/{genome}.gtf"
    params:
        url_fasta=lambda wc: "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_45/GRCh38.primary_assembly.genome.fa.gz" if wc.genome == "human" else "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/GRCm38.primary_assembly.genome.fa.gz",
        url_gtf=lambda wc: "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_45/gencode.v45.primary_assembly.annotation.gtf.gz" if wc.genome == "human" else "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/gencode.vM25.primary_assembly.annotation.gtf.gz"
    resources: mem_mb=4000, runtime=40
    shell:
        """
        mkdir -p references
        wget -O {output.fasta}.gz {params.url_fasta}
        gunzip {output.fasta}.gz
        wget -O {output.gtf}.gz {params.url_gtf}
        gunzip {output.gtf}.gz
        """

rule build_star_index:
    input:
        fasta="references/{genome}.fa",
        gtf="references/{genome}.gtf"
    output: directory("references/star_index_{genome}")
    threads: 16
    resources: mem_mb=64000, runtime=100
    container: "docker://quay.io/biocontainers/star:2.7.4a--0"
    params:
        sjdbOverhang=100,
        genomeSAindexNbases=14,
        genomeChrBinNbits=18
    shell:
        """
        mkdir -p {output}
        STAR --runMode genomeGenerate \
             --genomeDir {output} \
             --genomeFastaFiles {input.fasta} \
             --sjdbGTFfile {input.gtf} \
             --sjdbOverhang {params.sjdbOverhang} \
             --genomeSAindexNbases {params.genomeSAindexNbases} \
             --genomeChrBinNbits {params.genomeChrBinNbits} \
             --runThreadN {threads}
        """