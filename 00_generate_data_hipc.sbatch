#!/bin/bash 
#SBATCH -J 00_generate_data_hipc
#SBATCH -o 00_generate_data_hipc_%A_%a.out

singularity exec -B /mnt:/mnt \
/mnt/FCS_local2/Brinkman\ group/Alice/singularities/rstudio_202101.im Rscript \
/mnt/FCS_local2/Brinkman\ group/Alice/flowMagic_data/src/00_generate_data_hipc.R
