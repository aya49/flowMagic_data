#!/bin/bash 
#SBATCH -J 04_flowLearn
#SBATCH -o 04_flowLearn.out

singularity exec -B /mnt:/mnt \
/mnt/FCS_local2/Brinkman\ group/Alice/singularities/rstudio_202101.im Rscript \
/mnt/FCS_local2/Brinkman\ group/Alice/flowMagic_data/04_flowLearn.R
