#!/bin/bash 

### Written on 2020-06-28
### By Daniel Yokosawa

# Calls script to read JSON file that contains the path to users with 100 users or more.
# The Rscript reads the filepaths_count.json and saves the files_to_analyze.csv 
Rscript /mnt/f/Brinkman\ group/COVID/data/code/2.1.list_files_to_analyze.R

# Identifies the number of lines in the files_to_analyze.csv, use this to create the number of sbatch arrays.
nvals=$(cat /mnt/f/Brinkman\ group/COVID/data/code/support_files/files_to_analyze.csv | wc -l) 
nvals=$(($nvals-1))

# Create the jobs in slurm
sbatch --array=1-$nvals 2.2.combine_gates.sbatch
