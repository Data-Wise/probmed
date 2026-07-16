#!/bin/bash
#SBATCH --job-name=gaugeboot
#SBATCH --partition=general
#SBATCH --array=1-64
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=12:00:00
#SBATCH --output=/users/dtofighi/gauge_boot/logs/gb_%A_%a.out
module load r/4.4.0-ytj2 2>/dev/null
export R_LIBS=$HOME/Rlib/4.4-gauge:$HOME/Rlib/4.4-a15   # gauge lib FIRST (feature-branch code)
mkdir -p $HOME/gauge_boot/logs
Rscript $HOME/gauge_boot/run_gauge_boot_grid.R
