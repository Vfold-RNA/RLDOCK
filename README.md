**RLDOCK**
-----------------
RLDOCK is a source-availabe program for predicting the binding site and binding pose for ligand-RNA binding.  
Author: Li-Zhen Sun, Yangwei Jiang, & Shi-Jie Chen   
Email: chenshi@missouri.edu   
Date: Sept 27, 2020

System Requirement
------------------

 Linux/Unix  
 gcc compiler >4.8 version 
 
Compilation
----------------- 
To make `RLDOCK`, type:
```Bash
cd RLDOCK  
bash install.sh
```
If make installed, just type:
```Bash
cd RLDOCK  
make
```

Command line options
----------------- 
```Bash
-i <receptor.mol2>  # The Mol2 file of RNA.  
-l <ligand.mol2>    # The Mol2 file of ligand conformers.  
-o <output prefix>  # Path and filename for output files.  
-n <thread number>  # Number of threads used for simulation.  
-r <reference ligand file>  #(optional)The Mol2 file of the ligand for RMSD calculation.  
 ```   
  example:
 ```Bash
  ./RLDOCK -i job1_RNA.mol2 -l job1_ligand.mol2 -o job1 -n 20 -r job1_ref_lig.mol2    
```
Example
-----------------
The necessary files for the example cases are in the file `Example`.
To run the example cases, type:
```Bash
 cd RLDOCK
  bash run_4FEL.sh # An example with 0 rotatable bond.  
  bash run_3F2T.sh # An example with 7 rotatable bonds.
 ```
