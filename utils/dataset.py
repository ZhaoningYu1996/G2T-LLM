import os
from torch_geometric.data import download_url
from torch.utils.data import Dataset
import requests
import pandas as pd
import json
from rdkit import Chem
from utils.smiles2tree import smiles_to_tree, custom_json_serializer
from utils.utils import sanitize_smiles
from tqdm import tqdm
from utils.config import AtomTypeEnumZinc, BondTypeEnumZinc, AtomTypeEnumQM9, BondTypeEnumQM9
from utils.format import create_mol_format
from collections import defaultdict

class G2TDataset(Dataset):
    def __init__(self, data_name, raw_data, valid_idx, output_path, num_samples:int, addHs: bool = False, kekulize: bool = False):
            self.raw_data = raw_data
            self.valid_idx = valid_idx
            if data_name == 'zinc250k':
                atom_format, _ = create_mol_format(AtomTypeEnumZinc, BondTypeEnumZinc)
                atom_type_enum = AtomTypeEnumZinc
            elif data_name == 'qm9':
                atom_format, _ = create_mol_format(AtomTypeEnumQM9, BondTypeEnumQM9)
                atom_type_enum = AtomTypeEnumQM9
            self.addHs = addHs
            self.kekulize = kekulize
            self.data_name = data_name
            self.prep_dataset()
            self.format_data(num_samples, atom_format, atom_type_enum, output_path)
            
    # Function to download the dataset if it doesn't exist
    def download_dataset(self, url: str, dest_path: str):
        if not os.path.exists(dest_path):
            print(f"Downloading dataset from {url} to {dest_path}")
            response = requests.get(url)
            response.raise_for_status()  # Will stop the script if the download fails
            with open(dest_path, 'wb') as f:
                f.write(response.content)
            print("Download complete.")
        else:
            print("Dataset already downloaded.")

    def get_atom_set(self):
        atom_set = set([])
        for data in tqdm(self.dataset):
            mol = Chem.MolFromSmiles(data)
            for atom in mol.GetAtoms():
                atom_set.add(atom.GetSymbol())
        return atom_set
    
    def format_data(self, num_samples, atom_format, atom_type_enum, output_path: str):
        print("Sampling and formatting data...")
        # Define the schema for the Atom class
        class_schema = json.dumps(atom_format.model_json_schema(), indent=0)
        data_list = []
        smiles_list = []
        start_set = defaultdict(int)

        for smiles in tqdm(self.dataset[:num_samples]):
            smiles = sanitize_smiles(smiles)
            if smiles is not None:
                json_data = smiles_to_tree(self.data_name, smiles, kekulize=self.kekulize, addHs=self.addHs)
                input_schema = json.dumps(atom_format(atom_id=0, atom_type=atom_type_enum(json_data.atom_type.value), bonds=[]).dict(), indent=0, default=custom_json_serializer)
                start_set[json_data.atom_type.value] += 1
                response_schema = json.dumps(json_data.dict(), indent=0, default=custom_json_serializer)
                data = {
                    "instruction": f"Please generate a valid tree structure molecule from the following starting structure. You must respond using JSON format, according to the following schema: {class_schema}.",
                    "input": f"{input_schema}",
                    "output": f"{response_schema}"
                }
                data_list.append(data)
                smiles_list.append(smiles)

        print(f"Number of processed data: {len(data_list)}")
        print(f"Sart set: {start_set}")

        # Write the start set into a json file
        if not os.path.exists(os.path.dirname(output_path+"/start_set/")):
            os.makedirs(os.path.dirname(output_path+"/start_set/"))
        with open(output_path+"/start_set/"+self.data_name+".json", 'w') as file:
            file.write(json.dumps(start_set, indent=0))

        # Write the data_list into a json file
        if not os.path.exists(os.path.dirname(output_path+"/json/")):
            os.makedirs(os.path.dirname(output_path+"/json/"))
        with open(output_path+"/json/"+self.data_name+".json", 'w') as file:
            file.write(json.dumps(data_list, indent=0))

        # Write the smiles_list into a txt file
        if not os.path.exists(os.path.dirname(output_path+"smiles/")):
            os.makedirs(os.path.dirname(output_path+"smiles/"))
        with open(output_path+"smiles/"+self.data_name+".txt", 'w') as file:
            for smiles in smiles_list:
                file.write(smiles + '\n')

    def prep_dataset(self):
        print("Preparing dataset...")
        self.dataset = []
        self.valid_set = []
        train_set = []
        for idx, row in tqdm(self.raw_data.iterrows(), total=self.raw_data.shape[0]):
            if self.data_name == 'zinc250k':
                smiles = row['smile']
            elif self.data_name == 'qm9':
                smiles = row['SMILES1']
            smiles = sanitize_smiles(smiles)
            if smiles is not None:
                if idx in self.valid_idx:
                    self.valid_set.append(smiles)
                else:
                    train_set.append(smiles)
        print(f"Number of valid data: {len(self.valid_set)}")
        print(f"Number of train data: {len(train_set)}")
        # Shuffle the train set
        import random
        random.shuffle(train_set)
        self.dataset = train_set


    def __len__(self):
         return len(self.data)


# if __name__ == "__main__":
#     file_path = 'datasets/zinc250k_train.json'
#     raw_file_path = 'datasets/raw/zinc250k_property.csv'
#     download_raw_url = 'https://www.dropbox.com/s/feo9qle74kg48gy/molecules.zip?dl=1'
#     valid_idx_path = 'datasets/raw/valid_idx_zinc250k.json'
#     download_index_url = None
#     dataset = ZINCDataset(file_path, raw_file_path, download_raw_url, valid_idx_path, download_index_url)