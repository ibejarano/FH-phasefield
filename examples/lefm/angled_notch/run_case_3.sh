gmsh -2 caso_3_beta.geo -format msh2 -o caso_3_beta.msh
dolfin-convert caso_3_beta.msh caso_3_beta.xml
python caso_3.py