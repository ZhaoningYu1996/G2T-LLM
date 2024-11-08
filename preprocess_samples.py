import os
import argparse
import pandas as pd
import json
from utils.dataset import SampleDataset
from utils.format import create_mol_format
from utils.config import AtomTypeEnumZinc, BondTypeEnumZinc, AtomTypeEnumQM9, BondTypeEnumQM9

argparser = argparse.ArgumentParser(description='Preprocess the data')
argparser.add_argument('--data_name', type=str, default='qm9', help='Name of the dataset')
argparser.add_argument('--start_dict', type=dict, default={"C": 7080, "O": 1602, "N": 1307, "F": 11})
argparser.add_argument('--num_samples', type=int, default=14000, help='Number of samples to process')
argparser.add_argument('--index', type=int, default=0)
argparser.add_argument('--output_path', type=str, default='data/sampled/', help='Path to save the processed data')
args = argparser.parse_args()

if args.data_name == 'zinc':
    raw_file_path = 'data/raw/zinc250k_property.csv'
    valid_idx_path = 'data/raw/valid_idx_zinc250k.json'
    raw_data = pd.read_csv(raw_file_path)
    with open(valid_idx_path, 'r') as file:
        valid_idx = set(json.load(file))
    atom_format, bond_format = create_mol_format(AtomTypeEnumZinc, BondTypeEnumZinc)
    atom_type_enum = AtomTypeEnumZinc
    key = 'smile'
elif args.data_name == 'qm9':
    raw_file_path = 'data/raw/qm9.csv'
    valid_idx_path = 'data/raw/valid_idx_qm9.json'
    raw_data = pd.read_csv(raw_file_path)
    with open(valid_idx_path, 'r') as file:
        valid_idx = json.load(file)["valid_idxs"]
    valid_idx = set([int(idx) for idx in valid_idx])
    atom_format, bond_format = create_mol_format(AtomTypeEnumQM9, BondTypeEnumQM9)
    atom_type_enum = AtomTypeEnumQM9
    key = 'SMILES1'
else:
    raise ValueError(f"Invalid dataset name: {args.data_name}")

output_path = args.output_path + args.data_name + '/' + str(args.num_samples) + '/' + str(args.index) + '/'

if not os.path.exists(output_path):
    os.makedirs(output_path)

SampleDataset(atom_format=atom_format, bond_format=bond_format, atom_type_enum=atom_type_enum, start_dict=args.start_dict, num_samples=args.num_samples, output_path=output_path)
