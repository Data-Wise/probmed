#!/bin/bash
#SBATCH --job-name=probmed-i6
#SBATCH --partition=general
#SBATCH --array=1,101
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=8G
#SBATCH --time=12:00:00
#SBATCH --account=2016507
#SBATCH --output=/users/dtofighi/probmed_i6_20261019/logs/i6_%A_%a.out
# Pilot: tasks 1 (Gaussian s1 N=50) and 101 (binary s3 N=50). Override the array on the sbatch line for the full grid.
source /etc/profile.d/modules.sh
module load r/4.4.0-ytj2
command -v Rscript >/dev/null || { echo "FATAL: Rscript not on PATH after module load"; exit 127; }
D=$HOME/probmed_i6_20261019
export I6_PMED_SRC=$D/code
export I6_CORES=${SLURM_CPUS_PER_TASK:-8}
Rscript $D/runner.R $SLURM_ARRAY_TASK_ID $D/parts/i6_${SLURM_ARRAY_TASK_ID}.csv
