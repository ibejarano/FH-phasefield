#!/bin/bash


casos=("std_dt1e-4" "std_dt2e-4" "std_dt5e-5")
base_dir="./results"

#for caso in "${casos[@]}"; do
for dir in "$base_dir"/*/; do
    caso=$(basename "$dir")
    echo "Enviando: $caso"
    pvpython extract_frpath.py "$caso" --force-offscreen-rendering
done