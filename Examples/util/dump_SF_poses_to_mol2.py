import sys
import os

def split_SF_poses(file):
    with open(file) as f:
        docked_lines = f.readlines()
    pose_content_list = []
    for line in docked_lines:
        if line.find('<pocket position>') != -1:
            pose_content_list.append([])
        if line[0] != '<':
            pose_content_list[-1].append(line)
    return pose_content_list

def read_SF_poses(file):
    mol_list = []
    pose_content_list = split_SF_poses(file)
    for pose_content in pose_content_list:
        mol_list.append({"atoms": [], "coords": []})
        for atom_str in pose_content:
            _, atom, _, x, y, z = atom_str.split()
            mol_list[-1]["atoms"].append(atom)
            mol_list[-1]["coords"].append((float(x), float(y), float(z)))
        assert(len(mol_list[-1]["atoms"]) == len(mol_list[-1]["coords"]))
    return mol_list

def get_mol2_template_dict(in_mol2):
    with open(in_mol2) as f:
        mol2_lines = f.readlines()
    block_dict = {'@<TRIPOS>MOLECULE': [], '@<TRIPOS>ATOM': [], '@<TRIPOS>BOND': []}
    block_id = None
    for line in mol2_lines:
        if line[0:17] == '@<TRIPOS>MOLECULE':
            if len(block_dict['@<TRIPOS>MOLECULE']) != 0:
                break
            block_id = '@<TRIPOS>MOLECULE'
            continue
        elif line[0:13] == '@<TRIPOS>ATOM':
            block_id = '@<TRIPOS>ATOM'
            continue
        elif line[0:13] == '@<TRIPOS>BOND':
            block_id = '@<TRIPOS>BOND'
            continue
        elif line[0:9] == '@<TRIPOS>':
            block_id = 'Other'
            continue

        if block_id != None:
            if line.strip() == '':
                if block_id == '@<TRIPOS>MOLECULE':
                    block_dict[block_id].append(line)
            else:
                block_dict[block_id].append(line)

    # assert(len(block_dict['@<TRIPOS>MOLECULE']) == 6)
    atom_num, bond_num = int(block_dict['@<TRIPOS>MOLECULE'][1].split()[0]), int(block_dict['@<TRIPOS>MOLECULE'][1].split()[1])
    assert(len(block_dict['@<TRIPOS>ATOM']) == (atom_num))
    assert(len(block_dict['@<TRIPOS>BOND']) == (bond_num))

    template_dict = {'@<TRIPOS>MOLECULE': [], '@<TRIPOS>ATOM': []}
    for line in block_dict['@<TRIPOS>ATOM']:
        items = line.strip().split()
        atom_id = int(items[0])
        atom_name = items[1]
        x, y, z = float(items[2]), float(items[3]), float(items[4])
        mol2_type = items[5]
        resi = items[6]
        resn = items[7]
        charge = items[8]
        if mol2_type[0] == 'H':
            continue
        template_dict['@<TRIPOS>ATOM'].append((atom_id, atom_name, x, y, z, mol2_type, resi, resn, charge))

    assert( len(template_dict['@<TRIPOS>ATOM']) == atom_num )
    template_dict['@<TRIPOS>MOLECULE'] = block_dict['@<TRIPOS>MOLECULE']
    template_dict['@<TRIPOS>BOND'] = block_dict['@<TRIPOS>BOND']
    return template_dict

def write_coords_to_mol2_template(coords, mol2_template_dict, save_path, fmode='w'):
    with open(save_path, fmode) as f:
        # first write @<TRIPOS>MOLECULE
        f.write('@<TRIPOS>MOLECULE\n')
        for line in mol2_template_dict['@<TRIPOS>MOLECULE']:
            f.write(line)

        # second write @<TRIPOS>ATOM
        f.write('@<TRIPOS>ATOM\n')
        for coord, items in zip(coords, mol2_template_dict['@<TRIPOS>ATOM']):
            f.write(f'{items[0]:>7d} {items[1]:<6s}    {coord[0]:>8.3f}  {coord[1]:>8.3f}  {coord[2]:>8.3f} {items[5]:<6s} {items[6]:<3s} {items[7]:<6s} {items[8]:>8s}\n')

        # third write @<TRIPOS>BOND
        f.write('@<TRIPOS>BOND\n')
        for line in mol2_template_dict['@<TRIPOS>BOND']:
            f.write(line)

def dump_SF_poses_to_mol2(mol_dict, mol2_template_dict, save_path, fmode='w'):
    atoms, coords = mol_dict['atoms'], mol_dict['coords']
    assert(len(atoms) == len(mol2_template_dict['@<TRIPOS>ATOM']))
    assert( '-'.join(atoms) == '-'.join([items[1] for items in mol2_template_dict['@<TRIPOS>ATOM']]) )
    write_coords_to_mol2_template(coords, mol2_template_dict, save_path, fmode=fmode)



##########################################################################
template_ligand = sys.argv[1]
in_pose_path = sys.argv[2]
out_pose_path = sys.argv[3]

out_pose_f = open(out_pose_path, 'w')

if __name__ == '__main__':
    # get mol2 template
    mol2_template_dict = get_mol2_template_dict(template_ligand)

    SF_mol_list = read_SF_poses(in_pose_path)
    print(f'load {len(SF_mol_list)} poses from {in_pose_path}')
    for mol_dict in SF_mol_list:
        dump_SF_poses_to_mol2(mol_dict, mol2_template_dict, out_pose_path, fmode='a')

    print('completed.')
out_pose_f.close()