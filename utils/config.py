from enum import Enum, IntEnum

class AtomTypeEnumZinc(str, Enum):
    N = "N"
    P = "P"
    F = "F"
    C = "C"
    Br = "Br"
    S = "S"
    I = "I"
    Cl = "Cl"
    O = "O"

class BondTypeEnumZinc(Enum):
    SINGLE = 'SINGLE'
    DOUBLE = 'DOUBLE'
    TRIPLE = 'TRIPLE'
    AROMATIC = 'AROMATIC'

class AtomTypeEnumQM9(Enum):
    O = 'O'
    C = 'C'
    F = 'F'
    N = 'N'

class BondTypeEnumQM9(Enum):
    SINGLE = 'SINGLE'
    DOUBLE = 'DOUBLE'
    TRIPLE = 'TRIPLE'
    AROMATIC = 'AROMATIC'