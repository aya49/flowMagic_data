#!/bin/bash 
#SBATCH -o ./Outdata/COVID_data_clustering_%A_%a.out

singularity exec -B /mnt/f/Brinkman\ group/:/project/ \
/mnt/f/Docker/BrinkmanLabSingularity/brinkman_lab_singularity_190730.im Rscript \
/project/COVID/data/code/data_clustering.R ${SLURM_ARRAY_TASK_ID}
