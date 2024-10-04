import argparse
import pandas as pd
import json
from utils.dataset import G2TDataset

argparser = argparse.ArgumentParser(description='Preprocess the data')
argparser.add_argument('--data_name', type=str, default='zinc250k', help='Name of the dataset')
argparser.add_argument('--num_samples', type=int, default=100, help='Number of samples to process')
argparser.add_argument('--output_path', type=str, default='data/', help='Path to save the processed data')
argparser.add_argument('--addHs', type=bool, default=False, help='Add hydrogens to the molecule')
argparser.add_argument('--kekulize', type=bool, default=False, help='Kekulize the molecule')
args = argparser.parse_args()

if args.data_name == 'zinc250k':
    raw_file_path = 'data/raw/zinc250k_property.csv'
    valid_idx_path = 'data/raw/valid_idx_zinc250k.json'
    raw_data = pd.read_csv(raw_file_path)
    with open(valid_idx_path, 'r') as file:
        valid_idx = set(json.load(file))
elif args.data_name == 'qm9':
    raw_file_path = 'data/raw/qm9.csv'
    valid_idx_path = 'data/raw/valid_idx_qm9.json'
    raw_data = pd.read_csv(raw_file_path)
    with open(valid_idx_path, 'r') as file:
        valid_idx = json.load(file)["valid_idxs"]
    valid_idx = set([int(idx) for idx in valid_idx])
else:
    raise ValueError(f"Invalid dataset name: {args.data_name}")

G2TDataset(data_name="zinc250k", raw_data=raw_data, valid_idx=valid_idx, output_path=args.output_path, num_samples=args.num_samples, addHs=False, kekulize=False)