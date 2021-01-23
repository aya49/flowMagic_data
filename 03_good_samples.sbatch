#!/bin/bash 
#SBATCH -J 03_good_samples
#SBATCH -o 03_good_samples.out
#SBATCH --nodelist=BIOINFO17L,BIOINFO22L,BIOINFO23L,BIOINFO25L,BIOINFO26L

singularity exec -B /:/ \
/home/ayue/projects/singularities/rstudio_202101.im Rscript \
/home/ayue/projects/flowMagic_data/src/03_good_samples.R
