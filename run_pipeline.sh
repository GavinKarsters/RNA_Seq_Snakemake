#!/bin/bash

mkdir -p containers

# Point Apptainer/Singularity to containers folder
export APPTAINER_CACHEDIR="$PWD/containers"
export SINGULARITY_CACHEDIR="$PWD/containers"

export APPTAINER_TMPDIR="$PWD/containers/tmp"
export SINGULARITY_TMPDIR="$PWD/containers/tmp"

mkdir -p "$PWD/containers/tmp"

# Default values
DRY_RUN=false
TOUCH_TARGET=""
UNLOCK=false
EXP_ID=""
RUN_DESEQ2=false

# Function to detect rules from Snakefile
get_rules() {
    echo "  Available Rules (copy these names for -t):"
    echo "  ------------------------------------------"
    grep "^rule " Snakefile | sed 's/rule //; s/://' | awk '{print "    - " $1}'
    echo ""
}

# Function to show help
show_help() {
    echo "################################################################"
    echo "              RNA-seq Pipeline Wrapper"
    echo "################################################################"
    echo ""
    echo "Usage: ./run_pipeline.sh -e EXP_ID [options]"
    echo ""
    echo "Required:"
    echo "  -e STRING   Experiment ID (e.g., KAA12315). Matches config/samples_KAA12315.tsv"
    echo ""
    echo "Options:"
    echo "  -d          Enable DESeq2 differential expression analysis"
    echo "              Requires: config/metadata_EXP_ID.tsv"
    echo "              Optional: config/comparisons_EXP_ID.tsv (all pairwise if absent)"
    echo "  -n          Dry-run (print what would happen, don't execute)"
    echo "  -t STRING   Touch mode. Marks specific rules as 'done' to avoid rerunning."
    echo "  -u          Unlock directory (if Snakemake crashed previously)"
    echo "  -h          Show this help message"
    echo ""

    echo "================================================================"
    echo "                     WORKFLOW INSTRUCTIONS"
    echo "================================================================"
    echo "STEP 1: SETUP"
    echo "   Edit 'setup_experiment.py' to define raw data paths & conditions."
    echo "   Then generate the sample sheet + metadata:"
    echo "     $ python setup_experiment.py -e KAA12315"
    echo ""
    echo "STEP 2: CHECK"
    echo "   Verify the generated files:"
    echo "     $ column -t -s \$'\t' config/samples_KAA12315.tsv | less -S"
    echo "     $ column -t -s \$'\t' config/metadata_KAA12315.tsv | less -S"
    echo ""
    echo "STEP 3: (Optional) COMPARISONS"
    echo "   Create: config/comparisons_KAA12315.tsv"
    echo "     Tab-separated pairs: reference<TAB>treatment (one per line)"
    echo "     If absent, all pairwise comparisons are performed."
    echo ""
    echo "STEP 4: RUN"
    echo "   Without DESeq2:  ./run_pipeline.sh -e KAA12315"
    echo "   With DESeq2:     ./run_pipeline.sh -e KAA12315 -d"
    echo ""

    echo "================================================================"
    echo "                     UTILITIES & VISUALIZATION"
    echo "================================================================"
    echo "1. Generate Rule Graph (Flowchart):"
    echo "     $ snakemake --rulegraph --config id=KAA12315 | dot -Tpng > pipeline_flow.png"
    echo ""
    echo "2. Generate HTML Report:"
    echo "     $ snakemake --report report.html --config id=KAA12315"
    echo ""

    echo "================================================================"
    echo "                     ADVANCED: PARTIAL RERUNS"
    echo "================================================================"
    echo "Scenario: You changed Mapping params, but don't want to re-Trim."
    echo ""
    echo "   1. Touch the trimming rules (mark them as 'fresh'):"
    echo "      $ ./run_pipeline.sh -e KAA12315 -t \"trim_pe trim_se\""
    echo ""
    echo "   2. Run the pipeline normally:"
    echo "      $ ./run_pipeline.sh -e KAA12315"
    echo ""
    
    get_rules
    exit 1
}

# Parse command line arguments
while getopts "e:t:dnuh" opt; do
    case ${opt} in
        e) EXP_ID=$OPTARG ;;
        d) RUN_DESEQ2=true ;;
        n) DRY_RUN=true ;;
        t) TOUCH_TARGET=$OPTARG ;;
        u) UNLOCK=true ;;
        h) show_help ;;
        *) show_help ;;
    esac
done

# Check if ID was provided
if [ -z "$EXP_ID" ]; then
    echo "Error: Experiment ID (-e) is required."
    echo "Try ./run_pipeline.sh -h for help."
    exit 1
fi

# Define Expected Files
SAMPLE_SHEET="config/samples_${EXP_ID}.tsv"
OUT_DIR="results_${EXP_ID}"
LOG_DIR="${OUT_DIR}/logs/slurm_logs"

# Validation
if [ ! -f "$SAMPLE_SHEET" ]; then
    echo "Error: Sample sheet not found at $SAMPLE_SHEET"
    exit 1
fi

# DESeq2 validation
if [ "$RUN_DESEQ2" = true ]; then
    METADATA_FILE="config/metadata_${EXP_ID}.tsv"
    if [ ! -f "$METADATA_FILE" ]; then
        echo "Error: DESeq2 enabled (-d) but metadata file not found at $METADATA_FILE"
        echo "Create it with columns: sample_name  replicate  condition"
        exit 1
    fi
    echo "DESeq2 enabled. Metadata: $METADATA_FILE"
    
    COMPARISONS_FILE="config/comparisons_${EXP_ID}.tsv"
    if [ -f "$COMPARISONS_FILE" ]; then
        echo "Comparisons file found: $COMPARISONS_FILE"
    else
        echo "No comparisons file found. All pairwise comparisons will be performed."
    fi
fi

# Activate Environment
source $(conda info --base)/etc/profile.d/conda.sh
conda activate snakemake-c

# Count jobs
N_SAMPLES=$(tail -n +2 $SAMPLE_SHEET | wc -l)

# Ensure Log dir exists
mkdir -p $LOG_DIR

# Base Command
CMD=(
    snakemake
    -s Snakefile
    --config id=$EXP_ID run_deseq2=$RUN_DESEQ2
    --use-singularity
    --singularity-args "--cleanenv --bind /hpc"
    --rerun-incomplete
    --printshellcmds
    --restart-times 3
    --latency-wait 60
)

# --- HANDLE MODES ---

# 1. Unlock
if [ "$UNLOCK" = true ]; then
    echo "Unlocking directory..."
    "${CMD[@]}" --unlock
    exit 0
fi

# 2. Touch Mode (Specific Rules)
if [ ! -z "$TOUCH_TARGET" ]; then
    echo "----------------------------------------------------------------"
    echo "WARNING: Touch Mode selected for rule(s): '$TOUCH_TARGET'"
    echo "This will update the timestamp of these rules outputs to NOW."
    echo "Running on 1 core locally..."
    echo "----------------------------------------------------------------"
    
    "${CMD[@]}" --touch --cores 1 -R $TOUCH_TARGET
    
    echo ""
    echo "Touch complete. Now run the pipeline normally to process the rest."
    exit 0
fi

# 3. Dry Run or Cluster Execution
echo "Experiment: $EXP_ID | Samples: $N_SAMPLES | DESeq2: $RUN_DESEQ2"

CLUSTER_CMD="sbatch \
    --partition=cpu \
    --time={resources.runtime} \
    --mem={resources.mem_mb}M \
    --cpus-per-task={threads} \
    --output=$LOG_DIR/slurm-%j.out \
    --error=$LOG_DIR/slurm-%j.err \
    --mail-type=FAIL \
    --mail-user=g.j.karsters@umcutrecht.nl \
    --gres=tmpspace:170G"

if [ "$DRY_RUN" = true ]; then
    echo "Performing Dry-Run..."
    "${CMD[@]}" -npr
else
    echo "Submitting to SLURM..."
    "${CMD[@]}" --jobs 50 --cluster "$CLUSTER_CMD" --cluster-status "$PWD/slurm_status.sh"
fi