#!/bin/bash

for i in {1..11}
do
    folder="r$i"
    echo "Processing folder: $folder"
    
    cd "$folder" || { echo "Failed to enter $folder"; exit 1; }

    python Tether-Fac.py
    python Crosslink-Fac.py
    python Crosslink-Cons.py
    python IndexChanger.py
    mpirun -np 4 lmp_mpi -in in.Nuc11-1
    python DumpRead_FS.py
    python Force.py
    cd ..
done

python kspr_auc.py


