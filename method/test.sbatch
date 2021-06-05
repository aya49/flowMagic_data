#!/bin/bash
#SBATCH -J test
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=4
#SBATCH --time=0-00:05:00
#SBATCH --output=log_test.out
#SBATCH --partition=debug
#SBATCH --qos=overcap

source ~/miniconda3/etc/profile.d/conda.sh
conda activate pytorch_seg
module load LIB/CUDA/10.2

hostname
echo $CUDA_AVAILABLE_DEVICES

nvcc --version
nvidia-smi

cd /home/aya43/SETR

srun python ../flowMagic_test/20210423_test.py
