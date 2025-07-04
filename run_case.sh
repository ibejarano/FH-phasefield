#!/bin/bash


python meshing.py --config data/gmsh_1.json --case_dir results/test_comd
python main.py --config data/gmsh_1.json --case_dir test_comd