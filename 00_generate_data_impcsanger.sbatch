#!/bin/bash 
#SBATCH -J 00_generate_data_impcsanger
#SBATCH -o 00_generate_data_impcsanger.out

singularity exec -B /mnt:/mnt \
/mnt/FCS_local2/Brinkman\ group/Alice/singularities/rstudio_202101.im Rscript \
/mnt/FCS_local2/Brinkman\ group/Alice/flowMagic_data/00_generate_data_impcsanger.R
