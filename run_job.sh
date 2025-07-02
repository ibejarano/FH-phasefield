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
H=10

# Comando a ejecutar — reemplazá esto por lo que quieras
python3 meshing.py gmsh_1 name=$NAME l_max=1.5 h=1e-3 H=$H
python3 main.py gmsh_1 name=$NAME H=$H t_max=0.2 l_max=1.5