gmsh -2 caso_4.geo -format msh2 -o caso_4.msh
dolfin-convert caso_4.msh caso_4.xml
python caso_4.py