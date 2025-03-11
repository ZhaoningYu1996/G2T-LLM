from rdkit import Chem
from rdkit.Chem import rdmolops
from utils.utils import check_valency, correct_mol, sanitize_smiles
from utils.kekulize import check_kekulize
from rdkit import rdBase

rdBase.DisableLog('rdApp.*')
ATOM_VALENCY = {6: 4, 7: 3, 8: 2, 9: 1, 15: 3, 16: 2, 17: 1, 35: 1, 53: 1}

# Define a simple function to add atoms and bonds to an RDKit molecule
# def add_atom_and_bonds(mol, atom_data):
#     atom_index = mol.AddAtom(Chem.Atom(atom_data['atom_name']))
#     if 'adjcency_atoms' in atom_data:
#         for i, adj_atom in enumerate(atom_data['adjcency_atoms']):
#             adj_index = add_atom_and_bonds(mol, adj_atom)
#             bond_type = Chem.rdchem.BondType.SINGLE  # Default to SINGLE
#             if 'bond_types' in atom_data and len(atom_data['bond_types']) > i:
#                 bond_type_str = atom_data['bond_types'][i].upper()
#                 if bond_type_str == "DOUBLE":
#                     bond_type = Chem.rdchem.BondType.DOUBLE
#                 elif bond_type_str == "TRIPLE":
#                     bond_type = Chem.rdchem.BondType.TRIPLE
#             mol.AddBond(atom_index, adj_index, bond_type)
#     return atom_index

# Define a function to add atoms and bonds to an RDKit molecule, with an atom index shows the index of the atom in the molecule
# atom_data format is like below:
# class Molecule(BaseModel):
#     atom_name: AtomEnum
#     atom_index: int
#     adjcency_atoms: List['Molecule']
#     bond_types: List[BondTypeEnum]

# def add_atom_and_bonds(mol, atom_data, atom_id_map):
#     # Check if the atom index has already in the molecule
#     if atom_data['atom_index'] in atom_id_map:
#         atom_index = atom_id_map[atom_data['atom_index']]
#     else:
#         atom_index = mol.AddAtom(Chem.Atom(atom_data['atom_name']))
#         atom_id_map[atom_data['atom_index']] = atom_index
    
#     if 'adjcency_atoms' in atom_data:
#         for i, adj_atom in enumerate(atom_data['adjcency_atoms']):
#             add_atom_and_bonds(mol, adj_atom, atom_id_map)
#             adj_index = atom_id_map[adj_atom['atom_index']]
#             bond_type = Chem.rdchem.BondType.SINGLE
#             if 'bond_types' in atom_data and len(atom_data['bond_types']) > i:
#                 bond_type_str = atom_data['bond_types'][i].upper()
#                 if bond_type_str == "DOUBLE":
#                     bond_type = Chem.rdchem.BondType.DOUBLE
#                 elif bond_type_str == "TRIPLE":
#                     bond_type = Chem.rdchem.BondType.TRIPLE
#             if atom_index != adj_index and mol.GetBondBetweenAtoms(atom_index, adj_index) is None:
#                 mol.AddBond(atom_index, adj_index, bond_type)

# # Convert a JSON molecule to SMILES
# def json2smiles(molecule_data, do_correct=False):
#     # Create a new empty RDKit molecule
#     mol = Chem.RWMol()

#     # Add atoms and bonds to the molecule
#     add_atom_and_bonds(mol, molecule_data)

#     # Check if the molecule is valid
#     flag, _ = check_valency(mol)
#     if not flag:
#         if do_correct:
#             mol, _ = correct_mol(mol)
#         else:
#             return None

#     # Convert to SMILES
#     smiles = Chem.MolToSmiles(mol)
#     mol = Chem.MolFromSmiles(smiles)
#     smiles = Chem.MolToSmiles(mol)
#     if not smiles:
#         print("Invalid molecule")
#     return smiles


# Define a function to add atoms and bonds to an RDKit molecule
# def add_atom_and_bonds(mol, atom_data, prev_atom_index=None):
#     atom_index = mol.AddAtom(Chem.Atom(atom_data['atom']))
#     if 'bond_type' in atom_data and prev_atom_index is not None:
#         bond_type_str = atom_data['bond_type_list'][i].upper()
#         if bond_type_str == "DOUBLE":
#             bond_type = Chem.rdchem.BondType.DOUBLE
#         elif bond_type_str == "TRIPLE":
#             bond_type = Chem.rdchem.BondType.TRIPLE
#         else:
#             bond_type = Chem.rdchem.BondType.SINGLE
#         mol.AddBond(prev_atom_index, atom_index, bond_type)
#     if 'bonds' in atom_data:
#         for i, adj_atom in enumerate(atom_data['bonds']):
#             add_atom_and_bonds(mol, adj_atom, atom_index)
            
#     return atom_index

# # Convert a JSON molecule to SMILES
# def json2smiles(molecule_data, do_correct=False):
#     # Create a new empty RDKit molecule
#     mol = Chem.RWMol()

#     # Add atoms and bonds to the molecule
#     add_atom_and_bonds(mol, molecule_data)

#     # Check if the molecule is valid
#     flag, _ = check_valency(mol)
#     if not flag:
#         if do_correct:
#             mol, _ = correct_mol(mol)
#         else:
#             return None

#     # Convert to SMILES
#     smiles = Chem.MolToSmiles(mol)
#     mol = Chem.MolFromSmiles(smiles)
#     smiles = Chem.MolToSmiles(mol)
#     if not smiles:
#         print("Invalid molecule")
#     return smiles

# Define a function to add atoms and bonds to an RDKit molecule
# atom_data format is like below:
# class Atom(BaseModel):
#     atom_name: AtomEnum
#     bonds: List['Bond']

# class Bond(BaseModel):
#     atom: Atom
#     bond_type: BondTypeEnum

# def add_atom_and_bonds(mol, atom_data):
#     atom_index = mol.AddAtom(Chem.Atom(atom_data['atom_name']))
#     if 'bonds' in atom_data:
#         for i, bond in enumerate(atom_data['bonds']):
#             adj_index = add_atom_and_bonds(mol, bond['atom'])
#             bond_type_str = bond['bond_type'].upper()
#             if bond_type_str == "DOUBLE":
#                 bond_type = Chem.rdchem.BondType.DOUBLE
#             elif bond_type_str == "TRIPLE":
#                 bond_type = Chem.rdchem.BondType.TRIPLE
#             else:
#                 bond_type = Chem.rdchem.BondType.SINGLE
#             mol.AddBond(atom_index, adj_index, bond_type)
#     return atom_index

# class Atom(BaseModel):
#     atom_id: int
#     atom_name: AtomEnum
#     bonds: List['Bond'] = []

# class Bond(BaseModel):
#     atom: Atom
#     bond_type: BondTypeEnum

def add_atom_and_bonds(mol, atom_data):
    atom_index = mol.AddAtom(Chem.Atom(atom_data['atom_name']))
    if 'bonds' in atom_data:
        for i, bond in enumerate(atom_data['bonds']):
            adj_index = add_atom_and_bonds(mol, bond['atom'])
            bond_type_str = bond['bond_type'].upper()
            if bond_type_str == "DOUBLE":
                bond_type = Chem.rdchem.BondType.DOUBLE
            elif bond_type_str == "TRIPLE":
                bond_type = Chem.rdchem.BondType.TRIPLE
            elif bond_type_str == "SINGLE":
                bond_type = Chem.rdchem.BondType.SINGLE
            elif bond_type_str == "AROMATIC":
                bond_type = Chem.rdchem.BondType.AROMATIC
            mol.AddBond(atom_index, adj_index, bond_type)
    return atom_index

# class Atom(BaseModel):
#     atom_id: int
#     atom_name: AtomEnum
#     bonds: List['Bond'] = []

# class Bond(BaseModel):
#     atom: Atom
#     bond_type: BondTypeEnum

def add_atom(atom_data, atom_set):
    atom_id = atom_data['atom_id']
    atom_name = atom_data['atom_name']
    if atom_name != "H":
        atom_set.add((atom_id, atom_name))
        for i, bond in enumerate(atom_data['bonds']):
            atom_set = add_atom(bond['atom'], atom_set)
    return atom_set

def add_bond(mol, atom_data, atom_id_map, w_h=False):
    if atom_data['atom_name'] == "H":
        return
    curr_index = atom_data['atom_id']
    for i, bond in enumerate(atom_data['bonds']):
        if bond['atom']['atom_name'] == "H":
            continue
        adj_index = bond['atom']['atom_id']
        bond_type_str = bond['bond_type']
        if bond_type_str == 'DOUBLE':
            bond_type = Chem.rdchem.BondType.DOUBLE
        elif bond_type_str == 'TRIPLE':
            bond_type = Chem.rdchem.BondType.TRIPLE
        elif bond_type_str == 'SINGLE':
            bond_type = Chem.rdchem.BondType.SINGLE
        elif bond_type_str == "AROMATIC":
            bond_type = Chem.rdchem.BondType.AROMATIC
        # print(f"bond_type: {bond_type}, bond_type_str: {bond_type_str}")
        # mol.AddBond(curr_index, adj_index, bond_type)
        curr_idx = atom_id_map[curr_index]
        adj_idx = atom_id_map[adj_index]
        if curr_idx != adj_idx and mol.GetBondBetweenAtoms(curr_idx, adj_idx) is None:
            mol.AddBond(curr_idx, adj_idx, bond_type)
        add_bond(mol, bond['atom'], atom_id_map)

# Convert a JSON molecule to SMILES
def json2smiles(molecule_data, do_correct=False, prompt=None):
    # Create a new empty RDKit molecule
    mol = Chem.RWMol()

    # Add atoms and bonds to the molecule
    add_atom_and_bonds(mol, molecule_data)

    # for atom in mol.GetAtoms():
    #     print(atom.GetSymbol())
    # for bond in mol.GetBonds():
    #     print(bond.GetBondType())

    # Check if the molecule is valid
    flag, _ = check_valency(mol)
    if not flag:
        print(f"The original json molecule: {molecule_data}")
        print(f"The prompt: {prompt}")
        if do_correct:
            mol, _ = correct_mol(mol)
        else:
            return None

    # Convert to SMILES
    from rdkit.Chem.Draw import MolToImage
    image = MolToImage(mol)
    image.save(f'molecule_2.png')
    smiles = Chem.MolToSmiles(mol)
    mol = Chem.MolFromSmiles(smiles)
    smiles = Chem.MolToSmiles(mol)
    if not smiles:
        print("Invalid molecule")
    return smiles

def json2smiles_id(molecule_data, do_correct=False, prompt=None):
    # Create a new empty RDKit molecule
    mol = Chem.RWMol()

    # Add atoms and bonds to the molecule
    atom_set = set()
    atom_set = add_atom(molecule_data, atom_set)

    # Sort the atom_set by atom_id
    atom_set = sorted(atom_set, key=lambda x: x[0])
    atom_id_map = {}
    for i, (atom_id, atom_name) in enumerate(atom_set):
        atom_id_map[atom_id] = i
        mol.AddAtom(Chem.Atom(atom_name))

    add_bond(mol, molecule_data, atom_id_map)

    # for atom in mol.GetAtoms():
    #     print(atom.GetSymbol())
    # for bond in mol.GetBonds():
    #     print(bond.GetBondType())
    
    # from rdkit.Chem.Draw import MolToImage
    # image = MolToImage(mol)
    # image.save(f'molecule_3.png')

    # Check if the molecule is valid
    flag, atomid_valence = check_valency(mol)
    if not flag:
        # print(f"The original json molecule: {molecule_data}")
        # print(f"The prompt: {prompt}")
        if do_correct:
            mol, _ = correct_mol(mol)
        else:
            try:
                print(atomid_valence)
                assert len(atomid_valence) == 2
                idx = atomid_valence[0]
                v = atomid_valence[1]
                an = mol.GetAtomWithIdx(idx).GetAtomicNum()
                print(f"idx: {idx}, v: {v}, an: {an}")
                if an in (7, 8, 16) and (v - ATOM_VALENCY[an]) == 1:
                    mol.GetAtomWithIdx(idx).SetFormalCharge(1)
                    # Chem.SanitizeMol(mol)
                    # smiles = Chem.MolToSmiles(mol)
                    Chem.MolToSmiles(mol, isomericSmiles=True)
                    print("Molecule fixed after adjusting formal charge.")
                    return smiles
                
            except Exception as e:
                print(f"incorrect valence check")
                pass

    # Convert to SMILES
    
    try:
        mol = mol.GetMol()
        from rdkit.Chem import AllChem
        Chem.SanitizeMol(mol)
        Chem.Kekulize(mol, clearAromaticFlags=True)
        # rdmolops.Kekulize(mol)
        from rdkit.Chem.Draw import MolToImage
        image = MolToImage(mol)
        image.save(f'molecule_5.png')
        smiles = Chem.MolToSmiles(mol, isomericSmiles=True)
        smiles = sanitize_smiles(smiles)
    except Exception as e:
        print(f"Wrong mol to smiles")
        from rdkit.Chem.Draw import MolToImage
        image = MolToImage(mol)
        image.save(f'molecule_4.png')
        # smiles = None
        # Clear existing aromaticity flags and redefine aromaticity
        smiles = check_kekulize(mol)


    if not smiles:
        print("Invalid molecule")
    return smiles

def json2mol(molecule_data, do_correct=False, prompt=None):
    # Create a new empty RDKit molecule
    mol = Chem.RWMol()

    # Add atoms and bonds to the molecule
    add_atom_and_bonds(mol, molecule_data)

    # for atom in mol.GetAtoms():
    #     print(atom.GetSymbol())
    # for bond in mol.GetBonds():
    #     print(bond.GetBondType())

    # Check if the molecule is valid
    flag, _ = check_valency(mol)
    if not flag:
        # print(f"The original json molecule: {molecule_data}")
        # print(f"The prompt: {prompt}")
        if do_correct:
            mol, _ = correct_mol(mol)
        else:
            return None
        

    # mol = mol.GetMol()

    return mol