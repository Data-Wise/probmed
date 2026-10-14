#!/bin/bash
#SBATCH --job-name=pmedcalcol
#SBATCH --partition=general
#SBATCH --ntasks=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00
source /etc/profile.d/modules.sh
module load r/4.4.0-ytj2
export R_LIBS=${CALIB_R_LIBS:?set CALIB_R_LIBS}
cd ${CALIB_DIR:?set CALIB_DIR}
Rscript inst/sim/calibration/calib_collate.R parts results
