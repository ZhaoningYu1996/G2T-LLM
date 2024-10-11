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

class BondTypeEnumZinc(IntEnum):
    SINGLE = 1
    DOUBLE = 2
    TRIPLE = 3
    # AROMATIC = 'AROMATIC'

class AtomTypeEnumQM9(Enum):
    O = 'O'
    C = 'C'
    F = 'F'
    N = 'N'

class BondTypeEnumQM9(IntEnum):
    SINGLE = 1
    DOUBLE = 2
    TRIPLE = 3
    # AROMATIC = 'AROMATIC'