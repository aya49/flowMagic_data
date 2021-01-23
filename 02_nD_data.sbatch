#!/bin/bash 
#SBATCH -J 02_nD_data
#SBATCH -o 02_nD_data.out

singularity exec -B /mnt:/mnt \
/mnt/FCS_local2/Brinkman\ group/Alice/singularities/rstudio_202101.im Rscript \
/mnt/FCS_local2/Brinkman\ group/Alice/flowMagic_data/src/02_nD_data.R
