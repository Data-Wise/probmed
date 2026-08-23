#!/bin/bash
#SBATCH --job-name=gbcollate
#SBATCH --partition=general
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:20:00
#SBATCH --output=/users/dtofighi/gauge_boot/logs/collate_%j.out
# See submit_gauge_boot.sh: `module` is undefined in SLURM's non-login batch
# shell unless the lmod init is sourced first.
source /etc/profile.d/modules.sh
module load r/4.4.0-ytj2
command -v Rscript >/dev/null || { echo "FATAL: Rscript not on PATH after module load"; exit 127; }
export GAUGE_BOOT_DIR=${GAUGE_BOOT_DIR:-$HOME/gauge_boot}
export R_LIBS=${GAUGE_R_LIBS:-$HOME/Rlib/4.4-gauge:$HOME/Rlib/4.4-a15}
Rscript $GAUGE_BOOT_DIR/collate_gauge_boot.R
