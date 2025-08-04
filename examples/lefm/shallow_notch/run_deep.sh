gmsh -2 deep.geo -format msh2 -o deep.msh
dolfin-convert deep.msh deep.xml
python deep.py