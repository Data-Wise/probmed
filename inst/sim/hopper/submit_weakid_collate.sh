#!/bin/bash
#SBATCH --job-name=wvcollate
#SBATCH --partition=general
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=/users/dtofighi/weakid_val/logs/collate_%j.out
module load r/4.4.0-ytj2 2>/dev/null
export R_LIBS=$HOME/Rlib/4.4-gauge:$HOME/Rlib/4.4-a15
Rscript $HOME/weakid_val/collate_weakid_validation.R
