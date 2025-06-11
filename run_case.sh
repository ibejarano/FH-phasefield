#!/bin/bash


python meshing.py gmsh_1 name=test_frac_length l_max=0.1 h=1e-3 H=2
mpirun -n $1 python main.py gmsh_1 name=test_frac_length H=2 t_max=0.02