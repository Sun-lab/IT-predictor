#!/bin/bash
#SBATCH --job-name=step10
#SBATCH --ntasks=1
#SBATCH --time=5:00
#SBATCH --mem=1G

ml R
R CMD BATCH step10_annotate_var.R
