#!/bin/bash
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --mem=4G
#SBATCH --ntasks-per-node=1
#SBATCH --time=24:00:00
#SBATCH --partition=debug

# Cargar entorno de Anaconda
eval "$(conda shell.bash hook)"
conda activate fenicsproject
export OMP_NUM_THREADS=1
NAME=${SLURM_JOB_NAME}_job_${SLURM_JOB_ID}

python3 meshing.py --config data/gmsh_1.json --case_dir results/$NAME
python3 main.py --config data/gmsh_1.json --case_dir $NAME