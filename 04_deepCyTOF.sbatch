#!/bin/bash 
#SBATCH -J 04_deepCyTOF
#SBATCH -o 04_deepCyTOF.out
#SBATCH --nodelist=BIOINFO17L,BIOINFO22L,BIOINFO23L,BIOINFO25L,BIOINFO26L

singularity exec -B /mnt:/mnt \
/mnt/FCS_local2/Brinkman\ group/Alice/singularities/autokeras.im python \
/mnt/FCS_local2/Brinkman\ group/Alice/flowMagic_data/src/04_deepCyTOF.py
