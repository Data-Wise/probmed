#!/bin/bash
#SBATCH --job-name=pmedcal
#SBATCH --partition=general
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=12:00:00
# Calibration array worker. Submit from $CALIB_DIR (a copy of the package's
# inst/sim/calibration tree under inst/sim/calibration/), e.g.
#   cd $CALIB_DIR && CALIB_MODE=main CALIB_R_LIBS=$HOME/Rlib/4.4-repair:$HOME/Rlib/4.4-a15 \
#     sbatch --account=2016507 --array=1-1 --output=$CALIB_DIR/logs/%x_%A_%a.out \
#     inst/sim/calibration/hopper/submit_calib.sh          # pilot: ONE task first
# then --array=1-720 for main (36 cells x 20 chunks), --array=1-15 for stability.
# `module` is undefined in SLURM's non-login shell; source it explicitly.
source /etc/profile.d/modules.sh
module load r/4.4.0-ytj2
command -v Rscript >/dev/null || { echo "FATAL: Rscript not on PATH after module load"; exit 127; }
export R_LIBS=${CALIB_R_LIBS:?set CALIB_R_LIBS}
export CALIB_BASE_SEED=${CALIB_BASE_SEED:-20261014}
echo "base seed $CALIB_BASE_SEED, task $SLURM_ARRAY_TASK_ID"
cd ${CALIB_DIR:?set CALIB_DIR}
mkdir -p logs parts
Rscript inst/sim/calibration/calib_worker.R ${CALIB_MODE:-main} $SLURM_ARRAY_TASK_ID parts 20
