#!/bin/bash
#SBATCH -J nD
#SBATCH -o ./02_nD_data_%A_%a.out
#SBATCH -N 1 -n 15

singularity exec -B /mnt:/mnt \
/mnt/FCS_local3/backup/Brinkman\ group/current/Alice/singularities/rstudio.im Rscript /mnt/FCS_local3/backup/Brinkman\ group/current/Alice/flowMagic_data/src/02_nD_data.R
