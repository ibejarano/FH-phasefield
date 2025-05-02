#!/bin/bash
RESULTS_DIR="results/$1"
MESH_NAME="curving_H_min"

# Check if the results directory already exists
if [ -d "$RESULTS_DIR" ]; then
  # Ask the user for confirmation
  read -p "Directory '$RESULTS_DIR' already exists. Continue and potentially overwrite? [y/N]: " confirm
  # Convert confirmation to lowercase
  confirm_lower=$(echo "$confirm" | tr '[:upper:]' '[:lower:]')
  # If confirmation is not 'y', exit the script
  if [[ "$confirm_lower" != "y" ]]; then
    echo "Aborting."
    exit 1
  fi
fi

# Run the main python script, passing the case name argument
rm -rf $RESULTS_DIR
mkdir $RESULTS_DIR
gmsh -2 -format msh2 meshes/$MESH_NAME.geo -o $RESULTS_DIR/$MESH_NAME.msh
dolfin-convert $RESULTS_DIR/$MESH_NAME.msh $RESULTS_DIR/$MESH_NAME.xml
python main.py $1