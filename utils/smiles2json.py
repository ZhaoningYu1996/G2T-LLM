from rdkit import Chem
from rdkit.Chem import BondType
from pydantic import BaseModel
from utils.config import AtomTypeEnum, BondTypeEnum
from typing import List
import json
from enum import Enum
from collections import deque
# from utils.format import Molecule, Atom, Bond, Atom_id, Bond_id

# Bond mapping from RDKit to BondTypeEnum
rdkit_bond_mapping = {
    BondType.SINGLE: BondTypeEnum.SINGLE,
    BondType.DOUBLE: BondTypeEnum.DOUBLE,
    BondType.TRIPLE: BondTypeEnum.TRIPLE,
    # BondType.AROMATIC: BondTypeEnum.AROMATIC
}

# Define the Molecule class
# class Molecule(BaseModel):
#     atom_name: AtomEnum
#     adjacency_atoms: List['Molecule'] = []
#     bond_types: List[BondTypeEnum] = []

# Use a queue to traverse and build the molecule structure from SMILES
def smiles_to_molecule(smiles: str) -> Molecule:
    mol = Chem.MolFromSmiles(smiles)
    atom_idx_to_molecule = {}
    visited = set()
    queue = []

    # Initialize with the first atom
    for atom in mol.GetAtoms():
        if atom.GetIdx() == 0:  # Start with the first atom
            atom_molecule = Molecule(atom_name=AtomEnum(atom.GetSymbol()), adjacency_atoms=[], bond_types=[])
            atom_idx_to_molecule[atom.GetIdx()] = atom_molecule
            queue.append((None, atom.GetIdx(), None))  # (parent_idx, current_idx, bond_type)
            break

    # Process the queue using BFS
    while queue:
        parent_idx, current_idx, bond_type = queue.pop(0)
        if current_idx in visited:
            continue
        visited.add(current_idx)

        current_atom = atom_idx_to_molecule[current_idx]
        if parent_idx is not None:
            parent_atom = atom_idx_to_molecule[parent_idx]
            parent_atom.adjacency_atoms.append(current_atom)
            parent_atom.bond_types.append(bond_type)

        # Explore the neighbors
        for bond in mol.GetAtomWithIdx(current_idx).GetBonds():
            neighbor_idx = bond.GetOtherAtomIdx(current_idx)
            if neighbor_idx not in visited:
                neighbor_bond_type = BondTypeEnum(str(bond.GetBondType()).replace("BondType.", ""))
                if neighbor_idx not in atom_idx_to_molecule:
                    neighbor_atom = Molecule(atom_name=AtomEnum(mol.GetAtomWithIdx(neighbor_idx).GetSymbol()), adjacency_atoms=[], bond_types=[])
                    atom_idx_to_molecule[neighbor_idx] = neighbor_atom
                queue.append((current_idx, neighbor_idx, neighbor_bond_type))

    return atom_idx_to_molecule[0]

# smiles to molecule with different format, also store the bond that been removed
def smiles_to_atom(smiles: str, addHs=None) -> Atom:
    mol = Chem.MolFromSmiles(smiles)
    Chem.Kekulize(mol, clearAromaticFlags=True)
    if addHs:
        mol = Chem.AddHs(mol)
    atom_idx_to_atom = {}
    visited = set()
    queue = []
    added_bonds = set([])

    # Initialize with the first atom
    for atom in mol.GetAtoms():
        if atom.GetIdx() == 0:  # Start with the first atom
            atom_molecule = Atom(atom_name=AtomEnum(atom.GetSymbol()), bonds=[])
            atom_idx_to_atom[atom.GetIdx()] = atom_molecule
            queue.append((None, atom.GetIdx(), None))  # (parent_idx, current_idx, bond_type)
            break

    # Process the queue using BFS
    while queue:
        parent_idx, current_idx, bond_type = queue.pop(0)
        if current_idx in visited:
            continue
        visited.add(current_idx)

        current_atom = atom_idx_to_atom[current_idx]
        if parent_idx is not None:
            parent_atom = atom_idx_to_atom[parent_idx]
            bond = Bond(atom=current_atom, bond_type=bond_type)
            parent_atom.bonds.append(bond)
            added_bonds.add((parent_idx, current_idx))
            added_bonds.add((current_idx, parent_idx))

        # Explore the neighbors
        for bond in mol.GetAtomWithIdx(current_idx).GetBonds():
            neighbor_idx = bond.GetOtherAtomIdx(current_idx)
            if neighbor_idx not in visited:
                neighbor_bond_type = BondTypeEnum(str(bond.GetBondType()).replace("BondType.", ""))
                if neighbor_idx not in atom_idx_to_atom:
                    neighbor_atom = Atom(atom_name=AtomEnum(mol.GetAtomWithIdx(neighbor_idx).GetSymbol()), bonds=[])
                    atom_idx_to_atom[neighbor_idx] = neighbor_atom
                queue.append((current_idx, neighbor_idx, neighbor_bond_type))
        
    # Get bonds that not in added_bonds
    removed_bond = set([])
    for bond in mol.GetBonds():
        i = bond.GetBeginAtomIdx()
        j = bond.GetEndAtomIdx()
        if (i, j) not in added_bonds:
            removed_bond.add((i, j))
            removed_bond.add((j, i))
        
    return atom_idx_to_atom[0], removed_bond

# Convert a SMILES string to a JSON molecule object using only Atom class
def smiles_to_atom_json(smiles: str) -> Atom:
    mol = Chem.MolFromSmiles(smiles)
    atom_idx_to_atom = {}
    visited = set()
    queue = []
    added_bonds = set([])

    # Initialize with the first atom
    for atom in mol.GetAtoms():
        if atom.GetIdx() == 0:  # Start with the first atom
            atom_molecule = Atom(atom_name=AtomEnum(atom.GetSymbol()), neighbor_atoms=[])
            atom_idx_to_atom[atom.GetIdx()] = atom_molecule
            queue.append((None, atom.GetIdx()))  # (parent_idx, current_idx)
            break

    # Process the queue using BFS
    while queue:
        parent_idx, current_idx = queue.pop(0)
        if current_idx in visited:
            continue
        visited.add(current_idx)

        current_atom = atom_idx_to_atom[current_idx]
        if parent_idx is not None:
            parent_atom = atom_idx_to_atom[parent_idx]
            parent_atom.neighbor_atoms.append(current_atom)
            added_bonds.add((parent_idx, current_idx))
            added_bonds.add((current_idx, parent_idx))

        # Explore the neighbors
        for bond in mol.GetAtomWithIdx(current_idx).GetBonds():
            neighbor_idx = bond.GetOtherAtomIdx(current_idx)
            if neighbor_idx not in visited:
                if neighbor_idx not in atom_idx_to_atom:
                    neighbor_atom = Atom(atom_name=AtomEnum(mol.GetAtomWithIdx(neighbor_idx).GetSymbol()), neighbor_atoms=[])
                    atom_idx_to_atom[neighbor_idx] = neighbor_atom
                queue.append((current_idx, neighbor_idx))
        
    # Get bonds that not in added_bonds
    removed_bond = set([])
    for bond in mol.GetBonds():
        i = bond.GetBeginAtomIdx()
        j = bond.GetEndAtomIdx()
        if (i, j) not in added_bonds:
            removed_bond.add((i, j))
            removed_bond.add((j, i))
        
    return atom_idx_to_atom[0], removed_bond

def smiles_to_atom_w_id(smiles: str, kekulize=None, addHs=None) -> Atom:
    mol = Chem.MolFromSmiles(smiles)
    if kekulize:
        Chem.Kekulize(mol, clearAromaticFlags=True)

    if addHs:
        mol = Chem.AddHs(mol)
    atom_idx_to_atom = {}
    visited = set()
    queue = []
    added_bonds = set([])

    # print("====")
    # Initialize with the first atom
    for atom in mol.GetAtoms():
        if atom.GetIdx() == 0:  # Start with the first atom
            atom_molecule = Atom_id(atom_id=atom.GetIdx(), atom_name=AtomEnum(atom.GetSymbol()), bonds=[])
            atom_idx_to_atom[atom.GetIdx()] = atom_molecule
            queue.append((None, atom.GetIdx(), None))  # (parent_idx, current_idx, bond_type)
            break

    # Process the queue using BFS
    while queue:
        parent_idx, current_idx, bond_type = queue.pop(0)
        # if current_idx in visited:
        #     continue
        # visited.add(current_idx)

        current_atom = atom_idx_to_atom[current_idx]
        if parent_idx is not None:
            if current_idx in visited:
                current_atom = Atom_id(atom_id=current_idx, atom_name=AtomEnum(mol.GetAtomWithIdx(current_idx).GetSymbol()), bonds=[])
            parent_atom = atom_idx_to_atom[parent_idx]
            bond = Bond_id(atom=current_atom, bond_type=bond_type)
            parent_atom.bonds.append(bond)
            visited.add(current_idx)
            

        # Explore the neighbors
        for bond in mol.GetAtomWithIdx(current_idx).GetBonds():
            neighbor_idx = bond.GetOtherAtomIdx(current_idx)
            if (current_idx, neighbor_idx) not in added_bonds:
                # print(current_idx, neighbor_idx)

                neighbor_bond_type = BondTypeEnum(str(bond.GetBondType()).replace("BondType.", ""))
                # print(f"neighbor_bond_type: {neighbor_bond_type}, bond.GetBondType(): {str(bond.GetBondType())}")
                if neighbor_idx not in atom_idx_to_atom:
                    neighbor_atom = Atom_id(atom_id=neighbor_idx, atom_name=AtomEnum(mol.GetAtomWithIdx(neighbor_idx).GetSymbol()), bonds=[])
                    atom_idx_to_atom[neighbor_idx] = neighbor_atom
                
                queue.append((current_idx, neighbor_idx, neighbor_bond_type))
                # added_bonds.add((neighbor_idx, current_idx))
                added_bonds.add((current_idx, neighbor_idx))
        
    # Get bonds that not in added_bonds
    removed_bond = set([])
    for bond in mol.GetBonds():
        i = bond.GetBeginAtomIdx()
        j = bond.GetEndAtomIdx()
        if (i, j) not in added_bonds:
            removed_bond.add((i, j))
            removed_bond.add((j, i))
        
    return atom_idx_to_atom[0], removed_bond

# When serializing, convert enum to its value
def custom_json_serializer(obj):
    if isinstance(obj, Enum):
        return obj.value
    raise TypeError("Type not serializable")

# Example usage
if __name__ == "__main__":
    smiles = "C1CC1CO"  # Cyclopropane with a hydroxymethyl group
    molecule = smiles_to_molecule(smiles)
    molecule_json = json.dumps(molecule.dict(), indent=0, default=custom_json_serializer)
    print("========================================")
    print(molecule_json)