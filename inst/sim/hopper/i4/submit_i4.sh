#!/bin/bash
#SBATCH --job-name=probmed-i4
#SBATCH --partition=general
#SBATCH --array=1-1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=8G
#SBATCH --time=12:00:00
#SBATCH --account=2016507
#SBATCH --output=/users/dtofighi/probmed_i4_20261015/logs/i4_%A_%a.out
source /etc/profile.d/modules.sh
module load r/4.4.0-ytj2
command -v Rscript >/dev/null || { echo "FATAL: Rscript not on PATH after module load"; exit 127; }
D=$HOME/probmed_i4_20261015
export I4_PMED_SRC=$D/code
export I4_CORES=${SLURM_CPUS_PER_TASK:-8}
Rscript $D/runner.R $SLURM_ARRAY_TASK_ID $D/parts/i4_${SLURM_ARRAY_TASK_ID}.csv
