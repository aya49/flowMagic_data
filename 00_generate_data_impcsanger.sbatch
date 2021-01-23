#!/bin/bash 
#SBATCH -J 00_generate_data_impcsanger
#SBATCH -o 00_generate_data_impcsanger.out
#SBATCH --nodelist=BIOINFO17L,BIOINFO22L,BIOINFO23L,BIOINFO25L,BIOINFO26L

singularity exec -B /:/ \
/home/ayue/projects/singularities/rstudio.im Rscript \
/home/ayue/projects/flowMagic_data/src/00_generate_data_impcsanger.R
