import sys
import os
import json
import torch
import vllm
from vllm import LLM
import argparse
from typing import List
from pydantic import BaseModel
from tqdm import tqdm
from lmformatenforcer import JsonSchemaParser
from lmformatenforcer.integrations.vllm import build_vllm_logits_processor, build_vllm_token_enforcer_tokenizer_data
from transformers import AutoTokenizer
from utils.config import AtomTypeEnum, BondTypeEnum
from utils.tree2smiles import tree2smiles
# from utils.smiles2json import custom_json_serializer
# from utils.sample_fragment import sample_fragment
# from utils.smiles2json import smiles_to_atom_w_id
from utils.format import create_mol_format


# login(token='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')

argparser = argparse.ArgumentParser()
argparser.add_argument('--data_name', type=str, default='hiv')
argparser.add_argument('--data_path', type=str, default='data/sampled/long_input/zinc/500/0/samples.json')
argparser.add_argument('--batch_size', type=int, default=1)
argparser.add_argument('--num_samples', type=int, default=1000)
argparser.add_argument('--max_new_tokens', type=int, default=4000)
argparser.add_argument('--device', type=int, default=0)
argparser.add_argument('--output_path', type=str, default='results/hiv/')
args = argparser.parse_args()

def format(sample, column_map=None):
    template = {
        "prompt_input": (
            "Below is an instruction that describes a task, paired with an input that provides further context. "
            "Write a response that appropriately completes the request.\n\n"
            "### Instruction:\n{instruction}\n\n### Input:\n{input}\n\n### Response:\n"
        ),
    }
    column_map = column_map or {}
    key_input = column_map.get("input", "input")
    key_instruction = column_map.get("instruction", "instruction")

    prompt = template["prompt_input"].format(
        instruction=sample[key_instruction], input=sample[key_input]
    )
    return prompt

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

# Set the device
# device = args.device if torch.cuda.is_available() else -1

# print('loading model -------->', device)
model_path = 'tmp/Meta-Llama-3.1-8B-Instruct'

# tokenizer = AutoTokenizer.from_pretrained(model_path, use_fast=True)
# llm = vllm.LLM(model=model_path)
llm = vllm.LLM(model=model_path, tensor_parallel_size=1, max_model_len=5000, dtype=torch.bfloat16)
parser = JsonSchemaParser(mol_schema.model_json_schema())
tokenizer_data = build_vllm_token_enforcer_tokenizer_data(llm)
logit_processor = build_vllm_logits_processor(tokenizer_data, parser)
sampling_params = vllm.SamplingParams(max_tokens=5000, top_k=50, top_p=0.95, temperature=0.8, logits_processors=[logit_processor])

# Define the number of prompts to generate
max_iterations = 4000
valid_molecules = []
iteration = 0

path = 'results/'
if not os.path.exists(path):
    os.makedirs(path)

results = []
json_data_list = []
invalid_json = []
valid_molecules = []

for iter in tqdm(range(max_iterations), desc="Processing iterations"):
    # start_input = sample_fragment(atom_dict, "distributed")
    # input_json, _ = smiles_to_atom_w_id(start_input)
    # input_schema = json.dumps(input_json.dict(), indent=0, default=custom_json_serializer)
    sample = {
        "instruction": f'Please generate a valid hiv molecule with label 1. You must respond using JSON format, according to the following schema: {class_schema}.',
        "input": "",
    }
    prompt = format(sample)
    output_dict = llm.generate(prompt, sampling_params=sampling_params)
    result = output_dict[0].outputs[0].text

    # Convert the json molecule to SMILES
    try:
        molecule_data = json.loads(result)
        json_data_list.append(molecule_data)
        path = 'results/hiv/json/hiv_rand_no_input_10_sft_72_label_1/'
        if not os.path.exists(os.path.dirname(path)):
            os.makedirs(os.path.dirname(path))
        with open(path+str(iter)+'.json', 'w') as file:
            json.dump(molecule_data, file, indent=0)
    except Exception as e:
        print(e)
        print('Invalid JSON')
        invalid_json.append(result)
        continue
    smiles = tree2smiles(molecule_data, do_correct=False)
    if smiles:
        valid_molecules.append(smiles)
        path = args.output_path + 'smiles/hiv_rand_no_input_10_sft_72_label_1.txt'
        with open(path, 'a') as f:
            f.write(smiles + '\n')

print(f"Number of valid generated json data: {len(json_data_list)}")
print(f"Number of valid generated smiles: {len(valid_molecules)}")
