#!/bin/bash 
#SBATCH -J 05_deepCyTOF_score
#SBATCH -o 05_deepCyTOF_score.out
#SBATCH --nodelist=BIOINFO17L,BIOINFO22L,BIOINFO23L,BIOINFO25L,BIOINFO26L

singularity exec -B /:/ \
/home/ayue/projects/singularities/rstudio_202101.im Rscript \
/home/ayue/projects/flowMagic_data/src/05_deepCyTOF_score.R
