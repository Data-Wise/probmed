#!/bin/bash
#SBATCH --job-name=gaugeboot
#SBATCH --partition=general
#SBATCH --array=1-64
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=12:00:00
#SBATCH --output=/users/dtofighi/gauge_boot/logs/gb_%A_%a.out
# NOTE: this originally read `module load r/4.4.0-ytj2 2>/dev/null`. That works
# ONLY when sbatch is invoked from a shell that has already module-loaded R (the
# submitter's PATH is exported by default), which is how this grid was in fact run.
# Submitted from a plain ssh command it fails silently -- `module` is undefined in
# SLURM's non-login batch shell -- and every task dies with exit 127 in 0s. Fixed
# here after that failure bit the weakid array (job 4277033, 192/192 FAILED).
source /etc/profile.d/modules.sh
module load r/4.4.0-ytj2
command -v Rscript >/dev/null || { echo "FATAL: Rscript not on PATH after module load"; exit 127; }
export R_LIBS=$HOME/Rlib/4.4-gauge:$HOME/Rlib/4.4-a15   # gauge lib FIRST (feature-branch code)
mkdir -p $HOME/gauge_boot/logs
Rscript $HOME/gauge_boot/run_gauge_boot_grid.R
