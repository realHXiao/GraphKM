import os
import sys
import json
import paddle
import numpy as np
import xgboost as xgb
import pandas as pd
import pickle
import argparse

from rdkit import Chem
from rdkit.Chem import AllChem
from pahelix.utils.compound_tools import mol_to_graph_data
from src.model_for_xgb import KMModel
from src.data_gen import DataCollateFunc
from sklearn.metrics import mean_squared_error,r2_score
from sklearn.model_selection import train_test_split
from pahelix.utils.data_utils import load_npz_to_data_list
import matplotlib.pyplot as plt

def from_csv_and_csv_file(args): 
    KM_data = pd.read_csv(os.path.join(args.csv_file))
    n = 0

    data_unirep = []
    Unirep_df = pd.read_csv(os.path.join(args.input_unirep), sep = "\t")
    X_Unirep = Unirep_df.values
    for i in range(len(X_Unirep)):
        protein_ = X_Unirep[i,1:]
        data_unirep.append(protein_)
    KM_data['protein_unirep'] = data_unirep

    data_lst = []
    for ind in range(len(KM_data.values)) : 
        unirep = KM_data['protein_unirep'][ind]
        ligand = KM_data['Substrate_SMILES'][ind]
        
        if "." not in ligand:
            smiles = Chem.MolToSmiles(Chem.MolFromSmiles(ligand),
                            isomericSmiles=True)
            mol = AllChem.MolFromSmiles(smiles)
            mol_graph = mol_to_graph_data(mol)
            data0 = {s: v for s, v in mol_graph.items()}
            data0['protein_unirep'] = np.array([unirep])

            data_lst.append(data0) 

    return data_lst    

def xgb_input(args, pred_data, model_config, model):
    compound_protein = []
    for data in pred_data:
        infer_collate_fn = DataCollateFunc(
            model_config['compound']['atom_names'],
            model_config['compound']['bond_names'],
            is_inference=True,
            label_name=args.label_type)

        graphs, proteins_unirep = infer_collate_fn([data])
        graphs = graphs.tensor()
        proteins_unirep = paddle.to_tensor(proteins_unirep)

        model.eval()
        compound_protein_ = model(graphs, proteins_unirep)
        compound_protein.append(compound_protein_.tolist())

    return compound_protein

def input_csv_output_csv_file(args, pred):
    KM_data = pd.read_csv(os.path.join(args.csv_file))
    prediction_df = pd.DataFrame(columns= ["Protein_sequence", "Substrate_SMILES","KM_pred"])
    KM_pred = 10 ** pred
    protein_sequence, substrate = [], []
    for ind in range(len(KM_data.values)):
        ligand = KM_data['Substrate_SMILES'][ind]
        protein = KM_data['Protein_sequence'][ind]
        if "." not in ligand:
            protein_sequence.append(protein)
            substrate.append(ligand) 
    prediction_df['Protein_sequence'] = protein_sequence
    prediction_df['Substrate_SMILES'] = substrate
    prediction_df['{}_pred'.format(args.label_type)] = KM_pred
    
    return prediction_df

def from_fasta_and_csv_file(args, model_config, model): 
    ligand = open(os.path.join(args.SMILES_file), 'r')
    ligand = ligand.read()
    Unirep_df = pd.read_csv(os.path.join(args.input_unirep), sep = "\t")

    smiles = Chem.MolToSmiles(Chem.MolFromSmiles(ligand),
                isomericSmiles=True)
    mol = AllChem.MolFromSmiles(smiles)
    mol_graph = mol_to_graph_data(mol)

    X_Unirep = Unirep_df.values
    compound_protein = []
    for i in range(len(X_Unirep)):
        protein_unirep = X_Unirep[i,1:]

        data = {k: v for k, v in mol_graph.items()}
        data['protein_unirep'] = np.array(protein_unirep)

        infer_collate_fn = DataCollateFunc(
            model_config['compound']['atom_names'],
            model_config['compound']['bond_names'],
            is_inference=True,
            label_name=args.label_type)

        graphs, proteins_unirep = infer_collate_fn([data])
        graphs = graphs.tensor()
        proteins_unirep = paddle.to_tensor(proteins_unirep)

        model.eval()
        compound_protein_ = model(graphs, proteins_unirep)
        compound_protein.append(compound_protein_.tolist())

    return compound_protein

def input_fasta_output_csv_file(args, pred):
    protein_seq = open(os.path.join(args.fasta_file), 'r')
    lines = protein_seq.readlines()
    organism = []
    sequence, sequences = [], []
    n = 0
    for i in lines:
        n += 1
        line = i.strip('\n').split()
        word = list(line[0])
        if word[0] == '>':
            organism.append(i.strip('\n'))
        if word[0] != '>' and n > 1:
            sequence.append(line[0])
        else:
            seq = ''.join(sequence)
    #         print(seq)
            sequences.append(seq)
            sequence.clear()
            
    seq = ''.join(sequence)
    sequences.append(seq)

    prediction_df = pd.DataFrame(columns= ["organism", "protein_sequence", "substrate", "KM_pred"])
    prediction_df['organism'] = organism
    prediction_df['protein_sequence'] = sequences[1:]
    prediction_df['substrate'] = open(os.path.join(args.SMILES_file), 'r').read()

    pred_value = 10 ** pred
    prediction_df['KM_pred'] = pred_value

    return prediction_df

def main(args): 
    model_config = json.load(open(args.model_config, 'r'))

    # load model parameters. 
    model = KMModel(model_config)
    model.set_state_dict(paddle.load(os.path.join(args.model_file)))
    
    # load trained XGBoost model. 
    model_xgb_bst =  pickle.load(open(os.path.join(args.xgboost_model_file), "rb"))
    
    if args.csv_file_input: 
        pred_data = from_csv_and_csv_file(args)
        # generate input for trained XGBoost model. 
        compound_protein = xgb_input(args, pred_data, model_config, model)
        # predict
        dpred = xgb.DMatrix(compound_protein)
    
    if args.fasta_file_input:    
        # protein sequences fasta file with one substrate SMILES code. 
        compound_protein = from_fasta_and_csv_file(args, model_config, model)
        dpred = xgb.DMatrix(compound_protein)
    
    pred = model_xgb_bst.predict(dpred)
    
    # save results in csv file.
    if args.csv_file_input: 
        prediction_df = input_csv_output_csv_file(args, pred)
        prediction_df.to_csv(args.output_csv_file)
    
    # save results in csv file.
    if args.fasta_file_input:
        prediction_df = input_fasta_output_csv_file(args, pred)
        prediction_df.to_csv(args.output_csv_file)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-l', '--label_type', type=str, default='KM', help = 'Optional! Default KM. choose KM or kcat.')
    parser.add_argument("--csv_file", type=str, help = 'load csv file containing protein sequences, substrate SMILES codes. Required after -c or --csv_file_input is specified.')
    parser.add_argument("-c", "--csv_file_input", action="store_true", default=False, help = 'Optional! Default False (Disable). ')
    parser.add_argument('-i_unirep', "--input_unirep", type=str, help = 'Required! Load csv file of unirep vectors from Unirep50 !')
    parser.add_argument("--model_config", type=str, help = 'Required!')
    parser.add_argument("-m", "--model_file", type=str, help = 'Required! Load model which is trained by train.py.')
    parser.add_argument("-xgb", "--xgboost_model_file", type=str, help = 'Required! Load model which is trained by train_xgb.py.')
    parser.add_argument("--fasta_file", type=str, help = 'load fasta file from BLAST. Required after -f or --fasta_file_input is specified.')
    parser.add_argument("-f", "--fasta_file_input", action="store_true", default=False, help = 'Optional! Default False (Disable). ')
    parser.add_argument("-S", "--SMILES_file", type=str, help = "Required! load txt file containing one type SMILES code.")
    parser.add_argument("-o", "--output_csv_file", type=str, default='prediction.csv', help = 'Optional! ')
    args = parser.parse_args()
    main(args)