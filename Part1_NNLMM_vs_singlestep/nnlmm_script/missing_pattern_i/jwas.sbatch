#!/bin/bash -l
#SBATCH --job-name=accuracy
#SBATCH --cpus-per-task=10
#SBATCH --mem=140G
#SBATCH --time=9-99:00:00
#SBATCH --partition=bmm
#SBATCH --mail-type=NONE
#SBATCH --mail-user=tjzhao@ucdavis.edu



module load julia


srun julia --threads 20 /group/qtlchenggrp/tianjing/omics_fix_all/missing_pattern_legarra/jwas.jl $1 $2 $3