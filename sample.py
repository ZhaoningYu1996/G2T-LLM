
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
from utils.smiles2tree import custom_json_serializer, smiles_to_tree
from utils.sample_fragment import sample_fragment
from utils.format import create_mol_format
from utils.config import AtomTypeEnumZinc, BondTypeEnumZinc, AtomTypeEnumQM9, BondTypeEnumQM9
from torchtune.datasets import alpaca_dataset
from transformers.pipelines.pt_utils import KeyDataset
from torch.utils.data import DataLoader
import datasets
from torchtune.models.llama3 import llama3_tokenizer

login(token = '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')



argparser = argparse.ArgumentParser()
argparser.add_argument('--data_name', type=str, default='qm9')
argparser.add_argument('--data_path', type=str, default='data/sampled/qm9/1000/0/samples.json')
argparser.add_argument('--device', type=int, default=0)
argparser.add_argument('--output_path', type=str, default='results/qm9/')
args = argparser.parse_args()

# Set the device
device = args.device if torch.cuda.is_available() else -1

# Create a transformers pipeline
hf_pipeline = pipeline(
    'text-generation', 
    model='tmp/Meta-Llama-3.1-8B-Instruct',
    torch_dtype=torch.bfloat16,
    device_map="auto",
    max_new_tokens=1000,)

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
else:
    raise ValueError(f"Invalid data_name: {args.data_name}")

class_schema = json.dumps(mol_schema.model_json_schema(), indent=0)
# Create a character level parser and build a transformers prefix function from it
parser = JsonSchemaParser(mol_schema.schema())
prefix_function = build_transformers_prefix_allowed_tokens_fn(hf_pipeline.tokenizer, parser)

# Generate json format of the molecule
json_data = []
invalid_json = []

dataset = datasets.load_dataset('json', data_files=args.data_path, split='train')
json_data_list = []
valid_molecules = []
batch_size = 30

# for out in tqdm(hf_pipeline(KeyDataset(dataset, "text"), batch_size=batch_size, prefix_allowed_tokens_fn=prefix_function, do_sample=True, top_k=50, top_p=0.95, temperature=0.8, pad_token_id=hf_pipeline.tokenizer.eos_token_id, num_workers=batch_size), total=len(dataset)):
#     try:
#         molecule_data = json.loads(out[0]['generated_text'].split("Response:\n")[1])
#         json_data_list.append(molecule_data)
#     except Exception as e:
#         pass
    # smiles = tree2smiles(molecule_data, do_correct=False)
    # if smiles:
    #     valid_molecules.append(smiles)
    #     path = args.output_path + 'smiles/2k_94_2000.txt'
    #     with open(path, 'a') as f:
    #         f.write(smiles + '\n')

path = 'results/qm9/json/5k_50_1000_0.json'
with open(path, 'r') as file:
    json_data_list = json.load(file)

# print(f"Number of valid generated json data: {len(json_data_list)}")
# # Save the json data into a file
# path = args.output_path + 'json/5k_50_1000_0.json'
# with open(path, 'w') as file:
#     file.write(json.dumps(json_data_list, indent=0))

for molecule_data in json_data_list:
    smiles = tree2smiles(molecule_data, do_correct=False)
    if smiles:
        valid_molecules.append(smiles)
        path = args.output_path + 'smiles/5k_50_1000_0.txt'
        with open(path, 'a') as f:
            f.write(smiles + '\n')