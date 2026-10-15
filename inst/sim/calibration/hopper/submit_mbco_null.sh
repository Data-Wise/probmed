#!/bin/bash
#SBATCH --job-name=mbconull
#SBATCH --partition=general
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=06:00:00
# Composite-null MBCO harness (Phase 8 T2). Submit from $NULL_DIR (package
# inst/sim tree), one task per lattice point:
#   cd $NULL_DIR && NULL_R_LIBS=... sbatch --account=2016507 --array=1-1 \
#     --output=$NULL_DIR/logs/%x_%A_%a.out inst/sim/calibration/hopper/submit_mbco_null.sh
# `module` is undefined in SLURM's non-login shell; source it explicitly.
source /etc/profile.d/modules.sh
module load r/4.4.0-ytj2
command -v Rscript >/dev/null || { echo "FATAL: Rscript not on PATH after module load"; exit 127; }
export R_LIBS=${NULL_R_LIBS:?set NULL_R_LIBS}
cd ${NULL_DIR:?set NULL_DIR}
mkdir -p logs parts
Rscript inst/sim/calibration/calib_mbco_null.R $SLURM_ARRAY_TASK_ID parts 2000 50
