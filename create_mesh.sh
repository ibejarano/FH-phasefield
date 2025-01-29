#!/bin/bash
gmsh -2 -format msh2 meshes/$1.geo
dolfin-convert meshes/$1.msh meshes/$1.xml
