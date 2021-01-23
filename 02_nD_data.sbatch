#!/bin/bash 
#SBATCH -o 02_nD_data.out

/mnt/FCS_local2/Brinkman\ group/Alice/singularities
singularity exec -B /mnt:/mnt \
/mnt/FCS_local2/Brinkman\ group/Alice/singularities/rstudio_202101.im Rscript \
/mnt/FCS_local2/Brinkman\ group/Alice/flowMagic_data/02_nD_data.R
