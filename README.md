**RLDOCK**
-----------------
RLDOCK is a open source program for predicting the RNA-ligand binding site and binding pose.
Author: Li-Zhen Sun, Yangwei Jiang, Yuanzhe Zhou, & Shi-Jie Chen
Email: chenshi@missouri.edu
Date: Sept 27, 2020

System Requirements
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
or, if `make` installed, just type:
```Bash
cd RLDOCK
make
```
Executable file
-----------------
The executable file `RLDOCK` will be saved in the `bin` folder.


Command line options
-----------------
```Bash
-i <receptor.mol2>  # an RNA file in MOL2 format
-l <ligand.mol2>    # a ligand conformer file in MOL2 format
-o <output prefix>  # path and filename for output files
-c <number of output poses>  # number of output poses after clustering (if this number is larger than number of clusters, the number of output poses will be same as the number of clusters)
-n <thread number>  # number of threads used for simulation
-s <path of sphere.dat>  # sphere.dat records shafts used for rotating ligands and is stored in `src`
-r <reference ligand file>  # (optional) Mol2 file of the ligand for RMSD calculation
                            # if this file is not available, the first input conformer of the ligand will be set as the conformer for RMSD calculation
-b <center_x,center_y,center_z,size_x,size_y,size_z>    #  (optional) specify the center and the size of the docking box in Å, values separated by comma and no space between each value.
                                                        #  each size_x/y/z measures the whole side length of the box in the corresponding x/y/z direction
```
We recommend using [Chimera](http://www.cgl.ucsf.edu/chimera/) and [Open Babel](https://github.com/openbabel/openbabel/releases) to generate related files.
Important Note: The order of atoms in `<reference ligand file>` should be the same as the order in `<ligand.mol2>`.

example:
```Bash
./RLDOCK -i job1_RNA.mol2 -l job1_ligand.mol2 -o job1 -c 10 -n 20 -s src/sphere.dat -r job1_ref_lig.mol2
```

Output files
-----------------
There will be 6 output files for each simulation. `.xyz` and `.mol2` are for visualization. :
```Bash
XXX_pocket.dat   # Record the geometric score for potential binding sites.
XXX_pocket.xyz   # Potential binding sites are recorded in the format of .xyz for visualization.
XXX_usepose.dat  # Record the information of selected poses for scoring step.
XXX_SF_low.dat   # Record the scoring and ranking information by using the low resolution scoring function(SF-l).
XXX_SF_high.dat  # Record the scoring and ranking information by using the high resolution scoring function(SF-h).
XXX_SF_low.mol2  # Record the poses of XXX_SF_low.dat in mol2 format
XXX_SF_high.mol2 # Record the poses of XXX_SF_high.dat in mol2 format
XXX_cluster.mol2 # Top poses (default 10 poses) after clustering.
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
