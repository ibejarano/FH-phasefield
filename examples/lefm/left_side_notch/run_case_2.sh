gmsh -2 caso_2.geo -format msh2 -o caso_2.msh
dolfin-convert caso_2.msh caso_2.xml
python caso_2.py