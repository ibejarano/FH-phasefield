gmsh -2 caso_1.geo -format msh2 -o caso_1.msh
dolfin-convert caso_1.msh caso_1.xml
python caso_1.py