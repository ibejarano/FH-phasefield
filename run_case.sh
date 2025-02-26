#!/bin/bash
START_TIME=$SECONDS
MESH="curving_2"
OUTDIR="results/$1"

if [ -d "$OUTDIR" ];
then
    echo "$1 case already exists."
    echo -e "\nDelete case? ( [N] / y )"
    read ans
    if [ $ans == "y" ];
        then
            rm -rf $OUTDIR
        else
            echo "Script finished."
            exit
    fi
fi

mkdir -p $OUTDIR
LOG_FILE=$OUTDIR/output.log
cp "meshes/$MESH.geo" $OUTDIR
echo "Iniciando simulacion: " `date` | tee $LOG_FILE
gmsh -2 -format msh2 meshes/$MESH.geo -o $OUTDIR/$MESH.msh | tee -a $LOG_FILE
dolfin-convert $OUTDIR/$MESH.msh $OUTDIR/$MESH.xml | tee -a $LOG_FILE
python curving_hydraulic_fracture.py $1 $MESH | tee -a $LOG_FILE
ELAPSED_TIME=$(($SECONDS/60 - $START_TIME/60))
echo "Duracion: $ELAPSED_TIME mins"
echo "Simulacion Terminada: " `date` | tee -a $LOG_FILE