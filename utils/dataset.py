import os
from torch.utils.data import Dataset
import json
import csv
from rdkit import Chem
from tqdm import tqdm
from collections import defaultdict
from utils.smiles2tree import smiles_to_tree, custom_json_serializer
from utils.utils import sanitize_smiles
from utils.sample_fragment import sample_fragment
from torchtune.data import AlpacaInstructTemplate

class G2TDataset(Dataset):
    def __init__(self, data_name, raw_data, key, valid_idx, atom_format, bond_format, atom_type_enum, output_path, num_samples:int):
        self.data_name = data_name
        self.prep_dataset(raw_data, key, valid_idx)
        self.format_data(num_samples, atom_format, bond_format, atom_type_enum, output_path)

    def get_atom_set(self):
        atom_set = set([])
        for data in tqdm(self.dataset):
            mol = Chem.MolFromSmiles(data)
            for atom in mol.GetAtoms():
                atom_set.add(atom.GetSymbol())
        return atom_set
    
    def format_data(self, num_samples, atom_format, bond_format, atom_type_enum, output_path: str):
        print("Sampling and formatting data...")
        # Define the schema for the Atom class
        class_schema = json.dumps(atom_format.model_json_schema(), indent=0)
        data_list = []
        smiles_list = []
        start_set = defaultdict(int)

        for smiles in tqdm(self.dataset[:num_samples]):
            smiles = sanitize_smiles(smiles)
            if smiles is not None:
                json_data = smiles_to_tree(atom_format, bond_format, atom_type_enum, smiles)
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

    def prep_dataset(self, raw_data, key, valid_idx):
        print("Preparing dataset...")
        train_set = []
        for idx, row in tqdm(raw_data.iterrows(), total=raw_data.shape[0]):
            smiles = row[key]
            smiles = sanitize_smiles(smiles)
            if smiles is not None:
                if idx not in valid_idx:
                    train_set.append(smiles)
        print(f"Number of train data: {len(train_set)}")
        # Shuffle the train set
        import random
        random.shuffle(train_set)
        self.dataset = train_set


class SampleDataset(Dataset):
    def __init__(self, atom_format, bond_format, atom_type_enum, start_dict, num_samples, output_path):
        self.sample_dataset(atom_format, bond_format, atom_type_enum, start_dict, num_samples, output_path)

    def sample_dataset(self, atom_format, bond_format, atom_type_enum, start_dict, num_samples, output_path):
        print("Sampling data...")
        class_schema = json.dumps(atom_format.model_json_schema(), indent=0)
        samples_list = []
        for _ in tqdm(range(num_samples), desc='Generating molecules', total=num_samples):
            start_input = sample_fragment(start_dict, "distributed")
            input_json = smiles_to_tree(atom_format, bond_format, atom_type_enum, start_input)
            input_schema = json.dumps(input_json.dict(), indent=0, default=custom_json_serializer)
            sample = {
                "instruction": f'Please generate a valid molecule from the following starting structure. You must respond using JSON format, according to the following schema: {class_schema}.',
                "input": f"{input_schema}",
                "output": "",
            }
            prompt = AlpacaInstructTemplate.format(sample)
            samples_list.append(prompt)
        with open(output_path+'samples.json', 'w') as file:
            json_sample = json.dumps(samples_list, indent=0)
            file.write(json_sample)
        # with open(output_path+'samples.csv', 'w', newline='') as file:
        #     writer = csv.writer(file)
        #     for s in samples_list:
        #         writer.writerow([s]) 

# if __name__ == "__main__":
#     file_path = 'datasets/zinc250k_train.json'
#     raw_file_path = 'datasets/raw/zinc250k_property.csv'
#     download_raw_url = 'https://www.dropbox.com/s/feo9qle74kg48gy/molecules.zip?dl=1'
#     valid_idx_path = 'datasets/raw/valid_idx_zinc250k.json'
#     download_index_url = None
#     dataset = ZINCDataset(file_path, raw_file_path, download_raw_url, valid_idx_path, download_index_url)