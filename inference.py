import sys
import os
sys.path.append(os.path.abspath('/home/yu/projects/LLM-Mol/'))

from pydantic import BaseModel, fields, constr, StringConstraints
from lmformatenforcer import JsonSchemaParser
from lmformatenforcer.integrations.transformers import build_transformers_prefix_allowed_tokens_fn
from transformers import pipeline
from typing import Optional, Union, Any, ForwardRef, List
from typing_extensions import Annotated
import json
import torch
from huggingface_hub import login
from transformers import AutoTokenizer, AutoModelForCausalLM
from transformers import pipeline
# from ../config import AtomName, BondType
from utils.config import AtomEnum, BondTypeEnum
from utils.json2smiles import json2smiles, json2smiles_id
from tqdm import tqdm
from torchtune.data import AlpacaInstructTemplate
from utils.smiles2json import smiles_to_molecule, custom_json_serializer
from utils.sample_fragment import sample_fragment
from utils.smiles2json import smiles_to_atom_w_id
from utils.format import Atom_id

login(token = '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')

device = 0 if torch.cuda.is_available() else -1

# Create a transformers pipeline
hf_pipeline = pipeline(
    'text-generation', 
    model='tmp/Meta-Llama-3.1-8B-Instruct',
    torch_dtype=torch.bfloat16,
    device_map="auto")

atom_dict = {'C': 8249, 'N': 361, 'O': 1295}

path = 'results/'
if not os.path.exists(path):
    os.makedirs(path)

data_size = 500
valid_molecules = []

class_schema = json.dumps(Atom_id.model_json_schema(), indent=0)
# Create a character level parser and build a transformers prefix function from it
parser = JsonSchemaParser(Atom_id.model_json_schema())
prefix_function = build_transformers_prefix_allowed_tokens_fn(hf_pipeline.tokenizer, parser)

# Generate json format of the molecule
json_data = []
invalid_json = []

for i in tqdm(range(data_size), desc='Generating molecules', total=data_size):
    start_input = sample_fragment(atom_dict, "distributed")
    input_json, _ = smiles_to_atom_w_id(start_input)
    input_schema = json.dumps(input_json.dict(), indent=0, default=custom_json_serializer)
    sample = {
        "instruction": f'Please generate a valid molecule from the following starting structure. You must respond using JSON format, according to the following schema: {class_schema}.',
        "input": f"{input_schema}",
    }
    prompt = AlpacaInstructTemplate.format(sample)

    output_dict = hf_pipeline(
        prompt, prefix_allowed_tokens_fn=prefix_function, max_new_tokens=4000,
        do_sample=True, top_k=50, top_p=0.95, temperature=0.8, num_return_sequences=1, pad_token_id=hf_pipeline.tokenizer.eos_token_id)

    result = output_dict[0]['generated_text'][len(prompt):]

    try:
        molecule_data = json.loads(result)
        json_data.append(molecule_data)
        path = 'results/zinc/json/10k_dpo_2e5/'
        if not os.path.exists(os.path.dirname(path)):
            os.makedirs(os.path.dirname(path))
        with open(path+str(i)+'.json', 'w') as file:
            json.dump(molecule_data, file, indent=0)
    except Exception as e:
        print(e)
        print('Invalid JSON')
        invalid_json.append(i)
        continue
    smiles = json2smiles_id(molecule_data, do_correct=False)
    if smiles:
        valid_molecules.append(smiles)
        with open('results/zinc/smiles/10k_test.txt', 'a') as f:
            f.write(smiles + '\n')
    # else:
    #     print("None SMILES")

# Parallelly convert the json molecules to SMILES
# from concurrent.futures import ThreadPoolExecutor
# def process_chunk(chunk):
#     for molecule_data in chunk:
#         smiles = json2smiles(molecule_data, do_correct=False)
#         if smiles:
#             valid_molecules.append(smiles)
#             with open('results/valid_molecules_100_wo_correct_atom_format_fine_tune_5k_w_id_epoch_10_2.txt', 'a') as f:
#                 f.write(smiles + '\n')

# num_threads = 32
# chunk_size = data_size // (num_threads + 1)
# def process_parallel(data, chunk_size=chunk_size):
#     chunks = [data[i:i + chunk_size] for i in range(0, len(data), chunk_size)]
#     with ThreadPoolExecutor(max_workers=32) as executor:
#         list(tqdm(executor.map(process_chunk, chunks)))

# process_parallel(json_data)
print(invalid_json)
print(len(valid_molecules))
print(len(set(valid_molecules)))