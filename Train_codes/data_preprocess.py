"""
Convert dataset into npz file which can be trained directly.

"""

import os
import sys
import json
import random
import pickle
import argparse
import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from collections import OrderedDict

from pahelix.utils.compound_tools import mol_to_graph_data
from pahelix.utils.protein_tools import ProteinTokenizer
from pahelix.utils.data_utils import save_data_list_to_npz

def main(args):
    """Entry for data preprocessing."""
    tokenizer = ProteinTokenizer()
    data_input = os.path.join(args.input)
    seq_input = os.path.join(args.input_seq)
    output_file = os.path.join(args.output)

    # combinate the data of KM_data and seq_df into data_seq.
    KM_data = json.load(open(data_input, 'r'))
    
    n = 0
    data_list = []
    for da in KM_data : 
        n += 1
        da['protein_list'] = n
        data_list.append(da)
        
    data_seq = []
    seq_df = pd.read_csv(seq_input, sep = "\t")
    X_seq = seq_df.values
    for i in range(len(X_seq)):
        lis = i + 1
        for protein_ in data_list: 
            if lis == protein_['protein_list']:
                protein_['protein_seq'] = X_seq[i,1:]
                data_seq.append(protein_)
    
    # 
    data_lst = []
    for data in data_seq : 

        ligand = data['Smiles']
        protein = data['Sequence']
        KM = data['Value']
        unit = data['Unit']
        protein_seq = data['protein_seq']
        if "." not in ligand and float(KM) > 0:
            smiles = Chem.MolToSmiles(Chem.MolFromSmiles(ligand),
                            isomericSmiles=True)
            mol = AllChem.MolFromSmiles(smiles)
            mol_graph = mol_to_graph_data(mol)
            data0 = {s: v for s, v in mol_graph.items()}

            seqs = []
            for seq in protein:
                seqs.extend(tokenizer.gen_token_ids(seq))
            data0['protein_token_ids'] = np.array(seqs)

            # KM data
            KM = float(KM)
            if unit == 'mg/ml':
                MW = Descriptors.MolWt(mol)
                KM = (KM / MW) * 1000
            KM = np.log10(KM)
            data0[args.label_type] = np.array([KM])
            data0['protein_seq'] = np.array([protein_seq])

            data_lst.append(data0)

    npz = os.path.join(output_file)
    save_data_list_to_npz(data_lst, npz)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-l', '--label_type', type=str, default='KM', help = 'Optional! Default KM. choose KM or kcat.')
    parser.add_argument('-i', '--input', type=str, help = 'Required! json file containing protein sequences, substrate SMILES codes, KM values etc.')
    parser.add_argument('--input_seq', type=str, help = 'Required! csv file of protein sequence.')
    parser.add_argument('-o', '--output', type=str, help = 'Required! output file')
    args = parser.parse_args()
    main(args)
