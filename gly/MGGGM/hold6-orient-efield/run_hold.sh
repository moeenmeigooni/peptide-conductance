#!/bin/bash

n_reps=16
PYTHON_DIR="/Scr/meigoon2/miniconda3/bin"

for (( i=0; i<$n_reps; i++ )); do

echo "Running rep ${i}"
$PYTHON_DIR/python hold.py $i >& hold_${i}.log &

done

wait


