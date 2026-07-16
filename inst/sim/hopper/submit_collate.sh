#!/bin/bash
#SBATCH --job-name=gbcollate
#SBATCH --partition=general
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:20:00
#SBATCH --output=/users/dtofighi/gauge_boot/logs/collate_%j.out
module load r/4.4.0-ytj2 2>/dev/null
export R_LIBS=$HOME/Rlib/4.4-gauge:$HOME/Rlib/4.4-a15
Rscript $HOME/gauge_boot/collate_gauge_boot.R
