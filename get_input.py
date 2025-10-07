import os
import json
from utils.format import create_mol_format
from utils.config import AtomTypeEnumZinc, BondTypeEnumZinc, AtomTypeEnumQM9, BondTypeEnumQM9
from utils.smiles2tree import custom_json_serializer
import pandas as pd
from tqdm import tqdm

data_path = 'data/processed/zinc/10000/json/zinc.json'

with open(data_path, 'r') as file:
    data = json.load(file)

raw_file_path = 'data/raw/zinc250k_property.csv'
valid_idx_path = 'data/raw/valid_idx_zinc250k.json'
raw_data = pd.read_csv(raw_file_path)
with open(valid_idx_path, 'r') as file:
    valid_idx = set(json.load(file))
atom_type_enum = AtomTypeEnumZinc
bond_type_enum = BondTypeEnumZinc
atom_format, bond_format = create_mol_format(atom_type_enum, bond_type_enum)
key = 'smile'

input_list = []

for json_data in tqdm(data):
    json_data = json_data['output']
    json_data = json.loads(json_data)
    input = atom_format(atom_id=0, atom_type=atom_type_enum(json_data['atom_type']), bonds=[])
    for bond in json_data['bonds']:
        input.bonds.append(bond_format(atom=atom_format(atom_id=bond['atom']['atom_id'], atom_type=atom_type_enum(bond['atom']['atom_type']), bonds=[]), bond_type=bond['bond_type']))
        for inner_bond in bond['atom']['bonds']:
            input.bonds[-1].atom.bonds.append(bond_format(atom=atom_format(atom_id=inner_bond['atom']['atom_id'], atom_type=atom_type_enum(inner_bond['atom']['atom_type']), bonds=[]), bond_type=inner_bond['bond_type']))
    input_schema = json.dumps(input.dict(), indent=0, default=custom_json_serializer)
    input_list.append(input_schema)

output_path = 'data/processed/zinc/10000/'
# Write the input set into a json file
if not os.path.exists(os.path.dirname(output_path+"/input_list/")):
    os.makedirs(os.path.dirname(output_path+"/input_list/"))
with open(output_path+"/input_list/zinc.json", 'w') as file:
    file.write(json.dumps(input_list, indent=0))