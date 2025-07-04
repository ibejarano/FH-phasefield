#!/bin/bash
#SBATCH --job-name=fenics_job
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=12:00:00
#SBATCH --partition=debug

# Cargar entorno de Anaconda
eval "$(conda shell.bash hook)"
conda activate fenicsproject
export OMP_NUM_THREADS=1
NAME=slurm_job_${SLURM_JOB_ID}

python3 meshing.py --config data/gmsh_1.json --case_dir results/$NAME
python3 main.py --config data/gmsh_1.json --case_dir $NAME