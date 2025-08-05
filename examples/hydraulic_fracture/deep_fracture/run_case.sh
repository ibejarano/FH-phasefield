eval "$(conda shell.bash hook)"
conda activate fenicsproject
export OMP_NUM_THREADS=1
NAME=output

mkdir $NAME
cp deep_fh.geo $NAME
cp kgd_2d.py $NAME
gmsh -2 deep_fh.geo -format msh2 -o ./$NAME/deep_fh.msh
dolfin-convert ./$NAME/deep_fh.msh ./$NAME/deep_fh.xml
python kgd_2d.py $NAME