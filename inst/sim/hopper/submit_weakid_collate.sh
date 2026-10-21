#!/bin/bash
#SBATCH --job-name=wvcollate
#SBATCH --partition=general
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=/users/dtofighi/weakid_val/logs/collate_%j.out
# See submit_weakid_validation.sh: SLURM's batch shell is non-login, so the lmod
# init must be sourced or `module` is undefined and Rscript never reaches PATH.
source /etc/profile.d/modules.sh
module load r/4.4.0-ytj2
command -v Rscript >/dev/null || { echo "FATAL: Rscript not on PATH after module load"; exit 127; }
# See submit_weakid_validation.sh: 4.4-weakid (probmed >= 0.3.0) must precede
# 4.4-gauge, which still holds the pre-weak-ID 0.2.0.9000.
export R_LIBS=$HOME/Rlib/4.4-weakid:$HOME/Rlib/4.4-gauge:$HOME/Rlib/4.4-a15
Rscript $HOME/weakid_val/collate_weakid_validation.R
