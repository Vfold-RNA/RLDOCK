**RLDOCK**
-----------------
RLDOCK is a source-availabe program for predicting the binding site and binding pose for ligand-RNA binding.  
Author: Li-Zhen Sun, Yangwei Jiang, Yuanzhe Zhou, & Shi-Jie Chen   
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
Executable file
----------------- 
The executable file `RLDOCK` will be updated only if necessary.


Command line options
----------------- 
```Bash
-i <receptor.mol2>  # an RNA file in the format of MOL2  
-l <ligand.mol2>    # a ligand conformer file in the format of MOL2  
-o <output prefix>  # path and filename for output files  
-n <thread number>  # number of threads used for simulation  
-r <reference ligand file>  # (optional)the Mol2 file of the ligand for RMSD calculation.  
 ```
We recommend using [Chimera](http://www.cgl.ucsf.edu/chimera/) and [Open Babel](https://github.com/openbabel/openbabel/releases) to generate related files.  
Important Note: The order of atoms in `<reference ligand file>` should be the same as the order in `<ligand.mol2>`.
  
  example:
 ```Bash
  ./RLDOCK -i job1_RNA.mol2 -l job1_ligand.mol2 -o job1 -n 20 -r job1_ref_lig.mol2    
```

Output files
----------------- 
There will be 6 output files for each simulation:
```Bash
XXX_pocket.dat  # Record the geometric score for potential binding sites.  
XXX_pocket.xyz  # Potential binding sites are recorded in the format of .xyz for visualization. 
XXX_usepose.dat # Record the information of selected poses for scoring step.  
XXX_SF_low.dat  # Record the scoring and ranking information by using the low resolution scoring function(SF-l).  
XXX_SF_high.dat # Record the scoring and ranking information by using the high resolution scoring function(SF-h).
XXX_cluster.mol2 # Top 10 poses after clustering.
```

Example
-----------------
The necessary files for the example cases are in the file `Example`.
To run the example cases, type:
```Bash
 cd RLDOCK
 bash run_3GCA.sh # An example with 0 rotatable bond.  
 bash run_3F2T.sh # An example with 7 rotatable bonds.
 ```
