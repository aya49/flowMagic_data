#!/bin/bash 
#SBATCH -J 02_nD_data
#SBATCH -o 02_nD_data.out
#SBATCH --nodelist=BIOINFO17L,BIOINFO22L,BIOINFO23L,BIOINFO25L,BIOINFO26L

singularity exec -B /:/ \
/home/ayue/projects/singularities/rstudio_202101.im Rscript \
/home/ayue/projects/flowMagic_data/src/02_nD_data.R
