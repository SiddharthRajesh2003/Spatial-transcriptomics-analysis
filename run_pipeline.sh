#!/bin/bash

#SBATCH -J nf-spatial
#SBATCH -p general
#SBATCH -o run_pipeline_%j.txt
#SBATCH -e run_pipeline_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=sidrajes@iu.edu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=12
#SBATCH --time=48:00:00
#SBATCH --mem=10GB
#SBATCH -A r00750

echo "Starting Nextflow pipeline"
echo "Job ID: $SLURM_JOB_ID"
echo "Start time: $(date)"

base=/N/project/Krolab/Siddharth/Pipelines/spatial
cd $base

# Create output directory variable (fixed case)
OUTPUT_DIR=${base}/results                   # Fixed: Changed from lowercase 'output'
mkdir -p ${OUTPUT_DIR}/reports               # Create directory for reports

echo "Base directory: $base"
echo "Output directory: $OUTPUT_DIR"
export NXF_WORK="${base}/work"

# Load modules
module load apptainer
module load conda
module load java/21.0.8

conda activate ../../conda_env/envs/nf_env

echo "Starting pipeline execution..."

nextflow run main.nf \
    -resume \
    -profile slurm \
    -with-timeline ${OUTPUT_DIR}/reports/execution_timeline_${SLURM_JOB_ID}.html \
    -with-report ${OUTPUT_DIR}/reports/execution_report_${SLURM_JOB_ID}.html \
    -with-trace ${OUTPUT_DIR}/reports/execution_trace_${SLURM_JOB_ID}.txt \
    -with-dag ${OUTPUT_DIR}/reports/pipeline_dag_${SLURM_JOB_ID}.html

# Capture exit status
EXIT_STATUS=$?

echo ""
echo "Pipeline finished at $(date)"
echo "Exit status: $EXIT_STATUS"

if [ $EXIT_STATUS -eq 0 ]; then
    echo "SUCCESS: Pipeline completed successfully!"
    echo "Results available at: $OUTPUT_DIR"
    echo "Reports available at: ${OUTPUT_DIR}/reports/"
else
    echo "FAILED: Pipeline failed with exit status $EXIT_STATUS"
    echo "Check logs for details"
fi

exit $EXIT_STATUS