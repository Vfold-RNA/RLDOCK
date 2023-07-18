#!/bin/bash

./bin/RLDOCK -i ./Examples/3F2T/3F2T_RNA.mol2 -l ./Examples/3F2T/3F2T_ligand.mol2 -o ./Examples/3F2T/3F2T -c 10 -n 4 -s src/sphere.dat -r ./Examples/3F2T/3F2T_ref_lig.mol2

# ./bin/RLDOCK -i ./Examples/3F2T/3F2T_RNA.mol2 -l ./Examples/3F2T/3F2T_ligand.mol2 -o ./Examples/3F2T/3F2T -c 10 -n 4 -s src/sphere.dat -r ./Examples/3F2T/3F2T_ref_lig.mol2 -b 27.49,37.87,19.88,25,25,25
