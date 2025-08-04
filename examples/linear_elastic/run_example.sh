gmsh -2 plate_crack.geo -format msh2 -o plate_crack.msh
dolfin-convert plate_crack.msh mesh.xml
python ex_lefm_full.py