# Convert RLDOCK generated XXX_SF_low.dat and XXX_SF_high.dat to mol2 format

Run the following command in terminal to convert RLDOCK generated XXX_SF_low.dat and XXX_SF_high.dat to mol2 format:
```
python3 dump_SF_poses_to_mol2.py in_path_for_template_ligand_in_mol2_format in_path_for_SF_poses_dat_file out_path_for_SF_poses_in_mol2_format
```
**Note: the atom order in template ligand (mol2 format) should be the same as ones in SF_poses_dat_file. The template ligand file is neccessary since both the XXX_SF_low.dat and XXX_SF_high.dat do not contain enough information for mol2 format.**

Example case:
```
python3 dump_SF_poses_to_mol2.py ../3F2T/3F2T_ligand.mol2 ../3F2T/3F2T_SF_low.dat ../3F2T/3F2T_SF_low.mol2
```
