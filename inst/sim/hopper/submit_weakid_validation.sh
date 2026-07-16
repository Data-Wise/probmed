#!/bin/bash
#SBATCH --job-name=weakidval
#SBATCH --partition=general
#SBATCH --array=1-192
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=12:00:00
#SBATCH --output=/users/dtofighi/weakid_val/logs/wv_%A_%a.out
# 24 cells x 8 chunks = 192 tasks; 250 reps/chunk => nsim = 2000 per cell.
# R_LIBS must put a probmed >= PR #23 FIRST -- the weak-ID fields are required.
# The runner's stale-package guard aborts the task if they are missing.
module load r/4.4.0-ytj2 2>/dev/null
export R_LIBS=$HOME/Rlib/4.4-gauge:$HOME/Rlib/4.4-a15
mkdir -p $HOME/weakid_val/logs $HOME/weakid_val/parts
Rscript $HOME/weakid_val/run_weakid_validation.R
