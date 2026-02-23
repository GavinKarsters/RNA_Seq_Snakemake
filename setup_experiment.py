import os
import glob
import pandas as pd
import re
import argparse
import sys

# ==============================================================================
#                         USER CONFIGURATION
# ==============================================================================
# EDIT THESE SECTIONS PER EXPERIMENT
# ==============================================================================

# 1. WHERE ARE THE FASTQS?
RAW_DIR = "/hpc/local/Rocky8/umc_kaaij/gkarsters/Rna-seq/Janneke_boudewijn_tom20H2B_d_ala/Run_2/"

# 2. CONDITION MAPPING
# Map sample ID numbers (the 1, 2, 3 you added to the start of filenames) to conditions.
CONDITION_MAPPING = [
    # H2B Group
    # L-Ala is the WT/Reference
    {"samples": [1, 2, 3], "condition": "H2B_L_Ala", "replicates": [1, 2, 3]},
    # D-Ala is the Treatment
    {"samples": [4, 5, 6], "condition": "H2B_D_Ala", "replicates": [1, 2, 3]},

    # TOM20 Group
    # L-Ala is the WT/Reference
    {"samples": [7, 8, 9],   "condition": "TOM20_L_Ala", "replicates": [1, 2, 3]},
    # D-Ala is the Treatment
    {"samples": [10, 11, 12], "condition": "TOM20_D_Ala", "replicates": [1, 2, 3]},
]

# 3. DEFINE DESEQ2 COMPARISONS
# The R script will calculate: Treatment / Reference (Fold Change)
COMPARISONS = [
     # Treatment (D) vs Reference (L)
     ["H2B_D_Ala", "H2B_L_Ala"],
     
     # Treatment (D) vs Reference (L)
     ["TOM20_D_Ala", "TOM20_L_Ala"],
]
# ==============================================================================
#                   LOGIC (DO NOT EDIT)
# ==============================================================================

def get_args():
    parser = argparse.ArgumentParser(description="Generate samples TSV + metadata for RNA-seq pipeline")
    parser.add_argument("-e", "--experiment", required=True, help="Experiment ID (e.g. KAA12315)")
    parser.add_argument("-g", "--genome", default="human", help="Genome build (default: human)")
    parser.add_argument("--no-metadata", action="store_true", help="Skip metadata/comparisons file generation")
    return parser.parse_args()


def main():
    args = get_args()
    print(f"--- Setting up RNA-seq Experiment: {args.experiment} ---")
    print(f"--- Scanning Directory: {RAW_DIR} ---")

    if not os.path.exists(RAW_DIR):
        print(f"CRITICAL ERROR: Directory {RAW_DIR} does not exist.")
        sys.exit(1)

    files_map = {}

    for f in sorted(glob.glob(os.path.join(RAW_DIR, "*.gz"))):
        base = os.path.basename(f)

        # --- Try Pattern 1: Paired-end with read number ---
        # Matches: "100-SampleName_S1_R1_001.fastq.gz", "Sample_R2.fastq.gz", etc.
        pe_regex = re.compile(r"^(\d+)-(.+?)(?:_B23.+|_S\d+.+)?(?:_R|_)([12])(?:_001)?\.f(ast)?q\.gz$")
        match = pe_regex.match(base)

        if match:
            num = int(match.group(1))
            raw_name = match.group(2)
            name_clean = raw_name.split("_B23")[0].split("_S")[0]
            read_num = match.group(3)

            if num not in files_map:
                files_map[num] = {"name": f"{num}-{name_clean}", "R1": None, "R2": None}
            files_map[num][f"R{read_num}"] = f
            continue

        # --- Try Pattern 2: Single-end (no read number) ---
        # Matches: "1-2904F1_zfp869_dTAG_WT_rep1.fastq.gz"
        se_regex = re.compile(r"^(\d+)-(.+?)\.f(ast)?q\.gz$")
        match = se_regex.match(base)

        if match:
            num = int(match.group(1))
            raw_name = match.group(2)
            name_clean = raw_name.split("_B23")[0].split("_S")[0]

            if num not in files_map:
                files_map[num] = {"name": f"{num}-{name_clean}", "R1": None, "R2": None}
            files_map[num]["R1"] = f
            continue

    if not files_map:
        print("Error: No matching FastQ files found.")
        print("Expected patterns:")
        print("  Paired-end: <NUM>-<NAME>_R1_001.fastq.gz / <NUM>-<NAME>_R2_001.fastq.gz")
        print("  Single-end: <NUM>-<NAME>.fastq.gz")
        sys.exit(1)

    # Build condition lookup
    condition_lookup = {}
    for rule in CONDITION_MAPPING:
        for i, sample_num in enumerate(rule["samples"]):
            rep = rule["replicates"][i] if i < len(rule["replicates"]) else i + 1
            condition_lookup[sample_num] = {
                "condition": rule["condition"],
                "replicate": rep
            }

    data = []
    metadata_rows = []

    for num, info in sorted(files_map.items()):
        if not info["R1"]:
            continue

        cond_info = condition_lookup.get(num, {"condition": "FILL_IN", "replicate": "FILL_IN"})

        data.append({
            "sample_name": info["name"],
            "genome": args.genome,
            "fq1": info["R1"],
            "fq2": info["R2"] if info["R2"] else "NA"
        })

        metadata_rows.append({
            "sample_name": info["name"],
            "replicate": cond_info["replicate"],
            "condition": cond_info["condition"]
        })

    # --- Output sample sheet ---
    os.makedirs("config", exist_ok=True)

    sample_file = f"config/samples_{args.experiment}.tsv"
    df = pd.DataFrame(data).sort_values("sample_name")
    df.to_csv(sample_file, sep="\t", index=False)
    print(f"--- SUCCESS: Created {sample_file} ({len(df)} samples) ---")

    # Print preview
    print("\nSample sheet preview:")
    print(df.to_string(index=False))

    # --- Output metadata file (for DESeq2) ---
    if not args.no_metadata:
        metadata_file = f"config/metadata_{args.experiment}.tsv"
        meta_df = pd.DataFrame(metadata_rows).sort_values("sample_name")
        meta_df.to_csv(metadata_file, sep="\t", index=False)
        print(f"\n--- SUCCESS: Created {metadata_file} ---")

        print("\nMetadata preview:")
        print(meta_df.to_string(index=False))

        if any(meta_df["condition"] == "FILL_IN"):
            print("\n!! WARNING: Some samples have condition='FILL_IN'. Edit the metadata file before running DESeq2.")

        # --- Output comparisons file ---
        if COMPARISONS:
            comp_file = f"config/comparisons_{args.experiment}.tsv"
            with open(comp_file, "w") as fh:
                for pair in COMPARISONS:
                    fh.write(f"{pair[0]}\t{pair[1]}\n")
            print(f"\n--- SUCCESS: Created {comp_file} ({len(COMPARISONS)} comparisons) ---")
        else:
            print("\n--- INFO: No comparisons defined. All pairwise comparisons will be performed by DESeq2. ---")


if __name__ == "__main__":
    main()