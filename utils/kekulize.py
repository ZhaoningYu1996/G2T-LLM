from rdkit import Chem
from rdkit import rdBase
rdBase.DisableLog('rdApp.*')


def check_kekulize(mol):
    try:
        Chem.SanitizeMol(mol)
        print("Molecule is initially valid with no kekulization issues.")
        return Chem.MolToSmiles(mol)
    except Exception as e:
        print("Initial kekulization failed:", e)

    try:
        Chem.SanitizeMol(mol, sanitizeOps=Chem.SanitizeFlags.SANITIZE_ALL ^ Chem.SanitizeFlags.SANITIZE_SETAROMATICITY)
        Chem.SetAromaticity(mol)
        Chem.SanitizeMol(mol)
        print("Molecule fixed after adjusting aromaticity.")
        smiles = Chem.MolToSmiles(mol)
        return smiles
    except Exception as e:
        print("Kekulization failed after aromaticity adjustment:", e)
        
    try:
        for atom in mol.GetAtoms():
            atom.SetIsAromatic(False)
        for bond in mol.GetBonds():
            bond.SetIsAromatic(False)
        Chem.SanitizeMol(mol)
        smiles = Chem.MolToSmiles(mol)
        print("Molecule fixed after change aromatic.")
        return smiles
    except Exception as e:
        print("Final attempt to fix by modifying aromatic failed:", e)

    return None