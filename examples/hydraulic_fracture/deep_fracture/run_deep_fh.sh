gmsh -2 deep_fh.geo -format msh2 -o ./output/deep_fh.msh
dolfin-convert ./output/deep_fh.msh ./output/deep_fh.xml
python kgd_2d.py