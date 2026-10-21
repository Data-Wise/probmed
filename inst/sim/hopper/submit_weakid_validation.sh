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
# SLURM runs the batch script in a NON-login shell, where `module` is not defined
# unless the lmod init is sourced. The old form was `module load ... 2>/dev/null`,
# which silently swallowed "module: command not found", left Rscript off PATH, and
# killed all 192 tasks with exit 127 in 0s (job 4277033). It only ever worked
# because sbatch exports the submitter's environment (--export=ALL is the default),
# so a submit from a login shell that had already module-loaded R inherited the
# PATH. Submitting over plain ssh does not. Source the init explicitly and let the
# load fail loudly.
source /etc/profile.d/modules.sh
module load r/4.4.0-ytj2
command -v Rscript >/dev/null || { echo "FATAL: Rscript not on PATH after module load"; exit 127; }
# 4.4-weakid holds probmed >= 0.3.0 (the weak-ID fields) and MUST come first:
# 4.4-gauge still holds the 0.2.0.9000 that produced the coverage grid, which
# predates those fields -- it is kept on the path only for medfit/S7, and is left
# intact so that grid stays reproducible.
export R_LIBS=$HOME/Rlib/4.4-weakid:$HOME/Rlib/4.4-gauge:$HOME/Rlib/4.4-a15
mkdir -p $HOME/weakid_val/logs $HOME/weakid_val/parts
Rscript $HOME/weakid_val/run_weakid_validation.R
