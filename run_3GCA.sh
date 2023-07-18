#!/bin/bash

# ./bin/RLDOCK -i Examples/3GCA/3GCA_RNA.mol2 -l Examples/3GCA/3GCA_ligand.mol2 -o Examples/3GCA/3GCA -c 10 -n 4 -s src/sphere.dat -r Examples/3GCA/3GCA_ref_lig.mol2

./bin/RLDOCK -i Examples/3GCA/3GCA_RNA.mol2 -l Examples/3GCA/3GCA_ligand.mol2 -o Examples/3GCA/3GCA -c 10 -n 4 -s src/sphere.dat -r Examples/3GCA/3GCA_ref_lig.mol2 -b 35.05,33.23,24.74,20,20,20
