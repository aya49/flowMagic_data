#!/bin/bash
#SBATCH -J 2D
#SBATCH -o ./02_2D_data_%A_%a.out
#SBATCH --nodelist=BIOINFO25L

singularity exec -B /mnt:/mnt \
-H /mnt/FCS_local3/backup/Brinkman\ group/current/Alice \
/mnt/FCS_local3/backup/Brinkman\ group/current/Alice/singularities/rstudio.im Rscript /mnt/FCS_local3/backup/Brinkman\ group/current/Alice/flowMagic_data/src/02_2D_data.R
