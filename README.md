**RLDOCK**
-----------------
RLDOCK is a source-availabe program for predicting the binding site and binding pose for ligand-RNA binding. 
Author: Li-Zhen Sun, Yangwei Jiang, & Shi-Jie Chen Email: chenshi@missouri.edu Date: Sept 27, 2020

System Requirement
------------------

 Linux/Unix  
 gcc compiler >4.8 version 
 
Compilation
----------------- 
To make RLDOCK, type:
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
-i <receptor.mol2>  #The Mol2 file of RNA.  
-l <ligand.mol2>    #The Mol2 file of ligand conformers.  
-o <output prefix>  #Path and filename for output files.  
-n <thread number>  #Number of threads used for simulation.  
-r <reference ligand file>(optional)  #The Mol2 file of the ligand for RMSD calculation.  
 ```   
  example:
 ```Bash
  ./RLDOCK -i Example/1AJU_RNA.mol2 -l Example/1AJU_ligand.mol2 -o Example/1AJU -n 20 -r Example/1AJU_ref.mol2    
```
Examples
-----------------
The necessary files for the example case are in the file Example.
To run the example case, type:
  bash run_example.sh
