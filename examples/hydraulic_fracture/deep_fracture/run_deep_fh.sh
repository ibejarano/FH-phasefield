gmsh -2 deep_fh.geo -format msh2 -o deep_fh.msh
dolfin-convert deep_fh.msh deep_fh.xml
python kgd_2d.py