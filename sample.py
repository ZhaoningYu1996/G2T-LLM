
import os
import argparse
from lmformatenforcer import JsonSchemaParser
from lmformatenforcer.integrations.transformers import build_transformers_prefix_allowed_tokens_fn
from transformers import pipeline
import json
import torch
from huggingface_hub import login
from utils.tree2smiles import tree2smiles
from tqdm import tqdm
from utils.format import create_mol_format
from utils.config import AtomTypeEnumZinc, BondTypeEnumZinc, AtomTypeEnumQM9, BondTypeEnumQM9, AtomTypeEnum, BondTypeEnum
from transformers.pipelines.pt_utils import KeyDataset
import datasets
import pandas as pd
from collections import defaultdict
from utils.sample_start import sample_start
from utils.smiles2tree import custom_json_serializer
from torchtune.data import AlpacaInstructTemplate

login(token = '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')



argparser = argparse.ArgumentParser()
argparser.add_argument('--data_name', type=str, default='hiv')
argparser.add_argument('--data_path', type=str, default='data/sampled/long_input/zinc/500/0/samples.json')
argparser.add_argument('--batch_size', type=int, default=1)
argparser.add_argument('--num_samples', type=int, default=2000)
argparser.add_argument('--max_new_tokens', type=int, default=4000)
argparser.add_argument('--device', type=int, default=0)
argparser.add_argument('--output_path', type=str, default='results/hiv/')
args = argparser.parse_args()

# Set the device
device = args.device if torch.cuda.is_available() else -1

# Create a transformers pipeline
hf_pipeline = pipeline(
    'text-generation', 
    model='tmp/Meta-Llama-3.1-8B-Instruct',
    torch_dtype=torch.bfloat16,
    device_map="auto",
    max_new_tokens=args.max_new_tokens,)

hf_pipeline.tokenizer.pad_token_id = hf_pipeline.model.config.eos_token_id[0]

if not os.path.exists(args.output_path+'smiles'):
    os.makedirs(args.output_path+'smiles')
if not os.path.exists(args.output_path+'json'):
    os.makedirs(args.output_path+'json')

# Define the schema for the molecule
if args.data_name == 'zinc':
    mol_schema, _ = create_mol_format(AtomTypeEnumZinc, BondTypeEnumZinc)
elif args.data_name == 'qm9':
    mol_schema, _ = create_mol_format(AtomTypeEnumQM9, BondTypeEnumQM9)
elif args.data_name == 'hiv':
    mol_schema, _ = create_mol_format(AtomTypeEnum, BondTypeEnum)
else:
    raise ValueError(f"Invalid data_name: {args.data_name}")

class_schema = json.dumps(mol_schema.model_json_schema(), indent=0)
# Create a character level parser and build a transformers prefix function from it
parser = JsonSchemaParser(mol_schema.model_json_schema())
prefix_function = build_transformers_prefix_allowed_tokens_fn(hf_pipeline.tokenizer, parser)

# Generate json format of the molecule
# dataset = datasets.load_dataset('json', data_files=args.data_path, split='train')
json_data_list = []
invalid_json = []
valid_molecules = []

# data_path = 'data/processed/zinc/10000/json/zinc.json'

# with open(data_path, 'r') as file:
#     data = json.load(file)

# raw_file_path = 'data/raw/zinc250k_property.csv'
# valid_idx_path = 'data/raw/valid_idx_zinc250k.json'
# raw_data = pd.read_csv(raw_file_path)
# with open(valid_idx_path, 'r') as file:
#     valid_idx = set(json.load(file))
# atom_type_enum = AtomTypeEnumZinc
# bond_type_enum = BondTypeEnumZinc
# atom_format, bond_format = create_mol_format(atom_type_enum, bond_type_enum)
# key = 'smile'

# input_dict = defaultdict(int)

# for json_data in tqdm(data):
#     json_data = json_data['output']
#     json_data = json.loads(json_data)
#     input = atom_format(atom_id=0, atom_type=atom_type_enum(json_data['atom_type']), bonds=[])
#     for bond in json_data['bonds']:
#         input.bonds.append(bond_format(atom=atom_format(atom_id=bond['atom']['atom_id'], atom_type=atom_type_enum(bond['atom']['atom_type']), bonds=[]), bond_type=bond['bond_type']))
#         break
#         # for inner_bond in bond['atom']['bonds']:
#         #     input.bonds[-1].atom.bonds.append(bond_format(atom=atom_format(atom_id=inner_bond['atom']['atom_id'], atom_type=atom_type_enum(inner_bond['atom']['atom_type']), bonds=[]), bond_type=inner_bond['bond_type']))
#     input_schema = json.dumps(input.dict(), indent=0, default=custom_json_serializer)
#     input_dict[input_schema] += 1

# Remove input counts less than 20
# input_dict = {k: v for k, v in input_dict.items() if v >= 50}

for i in tqdm(range(args.num_samples), desc='Generating molecules', total=args.num_samples):
    # input_schema = sample_start(input_dict, "distributed")
    sample = {
        "instruction": f'Please generate a valid hiv molecule with label 0. You must respond using JSON format, according to the following schema: {class_schema}.',
        # "input": f"{input_schema}",
        "input": "",
    }
    prompt = AlpacaInstructTemplate.format(sample)

    output_dict = hf_pipeline(
        prompt, prefix_allowed_tokens_fn=prefix_function, max_new_tokens=5000,
        do_sample=True, top_k=50, top_p=0.95, temperature=0.8, num_return_sequences=1, pad_token_id=hf_pipeline.tokenizer.eos_token_id)

    result = output_dict[0]['generated_text'][len(prompt):]
    try:
        molecule_data = json.loads(result)
        json_data_list.append(molecule_data)
        path = 'results/hiv/json/hiv_rand_no_input_10_sft_72_label_0/'
        if not os.path.exists(os.path.dirname(path)):
            os.makedirs(os.path.dirname(path))
        with open(path+str(i)+'.json', 'w') as file:
            json.dump(molecule_data, file, indent=0)
    except Exception as e:
        print(e)
        print('Invalid JSON')
        invalid_json.append(result)
        continue
    smiles = tree2smiles(molecule_data, do_correct=False)
    if smiles:
        valid_molecules.append(smiles)
        path = args.output_path + 'smiles/hiv_rand_no_input_10_sft_72_label_0.txt'
        with open(path, 'a') as f:
            f.write(smiles + '\n')

print(f"Number of valid generated json data: {len(json_data_list)}")
print(f"Number of valid generated smiles: {len(valid_molecules)}")

# path = 'results/qm9/json/5k_50_1000_0.json'
# with open(path, 'r') as file:
#     json_data_list = json.load(file)

print(f"Number of valid generated json data: {len(json_data_list)}")
# Save the json data into a file
path = args.output_path + 'json/hiv_rand_no_input_10_sft_72_label_0.json'
with open(path, 'w') as file:
    file.write(json.dumps(json_data_list, indent=0))

# Save the invalid json data into a file
path = args.output_path + 'json/hiv_rand_no_input_10_sft_72_label_0_invalid.json'
with open(path, 'w') as file:
    file.write(json.dumps(invalid_json, indent=0))

# for molecule_data in json_data_list:
#     smiles = tree2smiles(molecule_data, do_correct=False)
#     if smiles:
#         valid_molecules.append(smiles)
#         path = args.output_path + 'smiles/10k_17_1000_0.txt'
#         with open(path, 'a') as f:
#             f.write(smiles + '\n')