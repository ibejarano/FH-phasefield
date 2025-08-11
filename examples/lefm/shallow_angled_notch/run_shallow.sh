gmsh -2 shallow.geo -format msh2 -o shallow.msh
dolfin-convert shallow.msh shallow.xml
python shallow.py