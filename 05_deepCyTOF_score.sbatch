#!/bin/bash 
#SBATCH -J 05_deepCyTOF_score
#SBATCH -o 05_deepCyTOF_score.out

singularity exec -B /mnt:/mnt \
/mnt/FCS_local2/Brinkman\ group/Alice/singularities/rstudio_202101.im Rscript \
/mnt/FCS_local2/Brinkman\ group/Alice/flowMagic_data/src/05_deepCyTOF_score.R
