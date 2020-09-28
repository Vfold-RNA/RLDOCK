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
If `make` installed, just type:
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
-r <reference ligand file>  #(optional)The Mol2 file of the ligand for RMSD calculation. The sequence of atoms should follow that of `<ligand.mol2>`.
 ```   
  example:
 ```Bash
  ./RLDOCK -i job1_RNA.mol2 -l job1_ligand.mol2 -o job1 -n 20 -r job1_ref_lig.mol2    
```

Output files
----------------- 
There will be 4 output files for each simulation:
```#Bash
XXX_pocket.dat   # Record the information of potential binding sites.  
XXX_usepose.dat  # Record the information of selected poses for scoring step.  
XXX_SF_low.dat   # Record the scoring and ranking information by using the low resolution scoring function(SF-l).  
XXX_SF_high.dat  # Record the scoring and ranking information by using the high resolution scoring function(SF-h).
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
