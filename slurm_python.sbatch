#!/bin/bash 
#SBATCH -J 04_deepCyTOF
#SBATCH -o 04_deepCyTOF.out
#SBATCH --nodelist=BIOINFO17L,BIOINFO22L,BIOINFO23L,BIOINFO25L,BIOINFO26L

singularity exec -B /:/ \
/home/ayue/projects/singularities/autokeras.im python \
/home/ayue/projects/flowMagic_data/src/04_deepCyTOF.py
