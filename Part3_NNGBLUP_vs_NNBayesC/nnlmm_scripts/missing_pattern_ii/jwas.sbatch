#!/bin/bash -l
#SBATCH --job-name=accuracy
#SBATCH --cpus-per-task=1
#SBATCH --mem=2500
#SBATCH --time=9-99:00:00
#SBATCH --partition=high
#SBATCH --mail-type=NONE
#SBATCH --mail-user=tjzhao@ucdavis.edu



module load julia


srun julia /group/qtlchenggrp/tianjing/res_andrew_1000snp/random_missing/jwas.jl $1 $2 $3