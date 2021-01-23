#!/bin/bash 
#SBATCH -o 03_good_samples.out

/mnt/FCS_local2/Brinkman\ group/Alice/singularities
singularity exec -B /mnt:/mnt \
/mnt/FCS_local2/Brinkman\ group/Alice/singularities/rstudio_202101.im Rscript \
/mnt/FCS_local2/Brinkman\ group/Alice/flowMagic_data/03_good_samples.R
