from pydantic import BaseModel
from typing import List

def create_mol_format(atom_type_enum, bond_type_enum):
    class Atom(BaseModel):
        atom_id: int
        atom_type: atom_type_enum
        bonds: List['Bond'] = []
    
    class Bond(BaseModel):
        atom: Atom
        bond_type: bond_type_enum
    
    return Atom, Bond