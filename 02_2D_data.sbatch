#!/bin/bash 
#SBATCH -J 02_2D_data
#SBATCH -o 02_2D_data.out
#SBATCH --nodelist=BIOINFO17L,BIOINFO22L,BIOINFO23L,BIOINFO25L,BIOINFO26L

singularity exec -B /:/ \
/home/ayue/projects/singularities/rstudio.im Rscript \
/home/ayue/projects/flowMagic_data/src/02_2D_data.R
