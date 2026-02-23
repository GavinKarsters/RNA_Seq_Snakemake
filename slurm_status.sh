#!/bin/bash
# Snakemake passes the full sbatch output, e.g. "Submitted batch job 46568348"
# We need to extract just the numeric job ID
jobid=$(echo "$1" | grep -oP '\d+$')

if [ -z "$jobid" ]; then
    echo "failed"
    exit 0
fi

state=$(sacct -j "$jobid" --format=State --noheader --parsable2 | head -1 | awk '{print $1}')

if [ -z "$state" ]; then
    state=$(squeue -j "$jobid" -h -o "%T" 2>/dev/null)
fi

case "$state" in
    COMPLETED)
        echo "success"
        ;;
    RUNNING|PENDING|CONFIGURING|COMPLETING|REQUEUED|SUSPENDED)
        echo "running"
        ;;
    *)
        echo "failed"
        ;;
esac