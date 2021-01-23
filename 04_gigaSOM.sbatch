#!/bin/bash 
#SBATCH -J 05_gigaSOM_score
#SBATCH -o 05_gigaSOM_score.out
#SBATCH --nodelist=BIOINFO17L,BIOINFO22L,BIOINFO23L,BIOINFO25L,BIOINFO26L

singularity exec -B /mnt:/mnt \
/mnt/FCS_local2/Brinkman\ group/Alice/singularities/julia.im julia \
/mnt/FCS_local2/Brinkman\ group/Alice/flowMagic_data/src/05_gigaSOM_score.jl
