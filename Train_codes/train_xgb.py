import os
import paddle
import json
import numpy as np
import xgboost as xgb
import pandas as pd
import pickle
import argparse

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from pahelix.utils.compound_tools import mol_to_graph_data
from src.model_for_xgb import KMModel
from src.data_gen import DataCollateFunc
from hyperopt import fmin, tpe, hp, Trials
from sklearn.metrics import mean_squared_error,r2_score
from sklearn.model_selection import train_test_split

def xgb_trials(dtest, dtrain, y_test):
    def cross_validation_mse_xgb_substrate_and_enzyme(param):

        num_round = param["num_rounds"]
        del param["num_rounds"]
        param["max_depth"] = int(np.round(param["max_depth"]))

        bst = xgb.train(param, dtrain, int(num_round), verbose_eval=False)
        y_pred = bst.predict(dtest)
        rmse = np.sqrt(mean_squared_error(y_test, y_pred))

        return(rmse)

    space = {"learning_rate": hp.uniform("learning_rate", args.learning_rate_min, args.learning_rate_max),
        "max_depth": hp.choice("max_depth", args.max_depth),
        "reg_lambda": hp.uniform("reg_lambda", args.reg_lambda_min, args.reg_lambda_max),
        "reg_alpha": hp.uniform("reg_alpha", args.reg_alpha_min, args.reg_alpha_max),
        "max_delta_step": hp.uniform("max_delta_step", args.max_delta_step_min, args.max_delta_step_max),
        "min_child_weight": hp.uniform("min_child_weight", args.min_child_weight_min, args.min_child_weight_max),
        "num_rounds":  hp.uniform("num_rounds", args.num_round_min, args.num_round_max),
        "tree_method": "gpu_hist"}

    trials = Trials()
    best = fmin(fn = cross_validation_mse_xgb_substrate_and_enzyme, space = space,
                algo=tpe.suggest, max_evals = 500, trials=trials)

    return best

def xgb_input_vector(args, model_config, model): 
    KM_data = json.load(open(args.input, 'r'))
    n = 0
    data_list = []
    for da in KM_data :
        n += 1
        da['protein_list'] = n
        data_list.append(da)

    data_unirep = []
    Unirep_df = pd.read_csv(os.path.join(args.input_unirep), sep = "\t")
    X_Unirep = Unirep_df.values
    for i in range(len(X_Unirep)):
        lis = i + 1
        for protein_ in data_list:
            if lis == protein_['protein_list']:
                protein_['protein_unirep'] = X_Unirep[i,1:]
                data_unirep.append(protein_)

    compound_protein, label_Y = [], []
    for data in data_unirep :

        ligand = data['Smiles']
        KM = data['Value']
        unit = data['Unit']
        unirep = data['protein_unirep']
        if "." not in ligand and float(KM) > 0:
            smiles = Chem.MolToSmiles(Chem.MolFromSmiles(ligand),
                            isomericSmiles=True)
            mol = AllChem.MolFromSmiles(smiles)
            mol_graph = mol_to_graph_data(mol)
            data = {s: v for s, v in mol_graph.items()}

            KM = float(KM)
            if unit == 'mg/ml':
                MW = Descriptors.MolWt(mol)
                KM = (KM / MW) * 1000
            KM = np.log10(KM)
            label_Y.append(KM)

            data['protein_unirep'] = np.array([unirep])

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

    return compound_protein, label_Y

def main(args): 
    model_config = json.load(open(args.model_config, 'r'))
    model = KMModel(model_config)
    model.set_state_dict(paddle.load(os.path.join(args.model_file)))

    # generate the input vectors for XGBoost model.
    compound_protein, label_Y = xgb_input_vector(args, model_config, model)

    # split the input vectors and labels. 
    X_train, X_test, y_train, y_test = train_test_split(compound_protein, label_Y, test_size=args.test_size, random_state=42)
    dtrain = xgb.DMatrix(X_train, label = y_train)
    dtest = xgb.DMatrix(X_test, label = y_test)
    evallist = [(dtest, 'eval'), (dtrain, 'train')]

    # test parameters
    best = xgb_trials(dtest, dtrain, y_test)

    # output the best parameters
    print("Trials: ", best)

    # add parameters before training XGBoost model. 
    param = best
    num_round = param["num_rounds"]
    del param["num_rounds"]
    param["tree_method"] = "gpu_hist"
    param["sampling_method"] = "gradient_based"

    # training XGBoost model. 
    bst = xgb.train(param, dtrain, int(num_round),evallist, verbose_eval=1)

    #Save model
    os.makedirs(args.model_dir, exist_ok=True)
    pickle.dump(bst, open(os.path.join(args.model_dir, "{}_xgboost_model.dat".format(model_config["compound"]["gnn_type"])), "wb"))
    
    # predict
    pred_test = bst.predict(dtest)
    
    # save true KM values and predicted KM values. 
    with open(os.path.join(args.model_dir, '{}_xgboost_test_dataset_label_pred.txt'.format(model_config["compound"]["gnn_type"])), 'w') as fi:
        for m in range(0, len(y_test)):
            fi.write('label: ' + str(y_test[m]) + ' pred: ' + str(pred_test[m]) + '\n')
    fi.close()

    # calculate MSE, RMSE, and R2. 
    GB_enzyme_sub_test_MSE = ((pred_test - np.reshape(y_test, (-1))) ** 2).mean(axis=0)
    GB_enzyme_sub_test_RMSE = np.sqrt(mean_squared_error(np.reshape(y_test, (-1)), pred_test))
    GB_enzyme_sub_test_R2 = r2_score(np.reshape(y_test, (-1)), pred_test)

    print("MSE: {}, RMSE: {}, R2: {}".format(GB_enzyme_sub_test_MSE, GB_enzyme_sub_test_RMSE, GB_enzyme_sub_test_R2))

    # plot correlation figure between true KM values and predicted KM values. 
    import matplotlib.pyplot as plt
    folder = './GNNs_xgboost_plt_image'
    if not os.path.exists(folder):
        os.makedirs(folder)
    plt.figure(figsize=(6, 6), dpi=300)
    plt.scatter(y_test, pred_test, alpha=0.1)  # scatter transparency is 90%
    plt.xlabel("label", fontdict={'size': 16})
    plt.ylabel("pred", fontdict={'size': 16})
    plt.savefig(os.path.join(folder, '{}_xgboost_plt_image.jpg'.format(model_config["compound"]["gnn_type"])))

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-l', '--label_type', type=str, default='KM', help = 'Optional! Default KM. choose KM or kcat.')
    parser.add_argument("-m", "--model_file", type=str, help = 'Required! load parameters file of trained model from train.py.')
    parser.add_argument("-i", "--input", type=str, help = 'Required! json file containing protein sequences, substrate SMILES codes, KM values etc.')
    parser.add_argument('-i_unirep', '--input_unirep', type=str, help = 'Required! csv file of unirep vectors from Unirep50 !')
    parser.add_argument("--test_size", type=float, default=0.2, help = 'Optional! Default 0.2. split whole dataset into train_dataset and test_dataset. if set 0.2, train_dataset = 80/100 whole dataset.')
    parser.add_argument("--random_state", type=int, default=42, help = 'Optional! Default 42. random assignment of train_dataset and test_dataset.')
    parser.add_argument("-lr_min", "--learning_rate_min", type=float, default=0.01, help = 'Optional! Default 0.01. the minimum learning rate for parameter testing of xgboost model.')
    parser.add_argument("-lr_max", "--learning_rate_max", type=float, default=1, help = 'Optional! Default 1. the maximum learning rate for parameter testing of xgboost model.')
    parser.add_argument("--max_depth", type=list, default=[3,4,5,6,7,8], help = 'Optional! Default [3,4,5,6,7,8]. choose several number for parameter testing. maximum depth of a tree. Increasing this value will make the model more complex and more likely to overfit. 0 indicates no limit on depth. Beware that XGBoost aggressively consumes memory when training a deep tree.')
    parser.add_argument("--reg_lambda_min", type=int, default=0, help = 'Optional! Default 0. the minimum reg_lambda value for parameter testing of xgboost model. L2 regularization term on weights. Increasing this value will make model more conservative.')
    parser.add_argument("--reg_lambda_max", type=int, default=5, help = 'Optional! Default 5. the maximum reg_lambda value for parameter testing of xgboost model. L2 regularization term on weights. Increasing this value will make model more conservative.')
    parser.add_argument("--reg_alpha_min", type=int, default=0, help = 'Optional! Default 0. the minimum reg_alpha value for parameter testing of xgboost model. L1 regularization term on weights. Increasing this value will make model more conservative.')
    parser.add_argument("--reg_alpha_max", type=int, default=5, help = 'Optional! Default 5. the maximum reg_alpha value for parameter testing of xgboost model. L1 regularization term on weights. Increasing this value will make model more conservative.')
    parser.add_argument("--max_delta_step_min", type=int, default=0, help = 'Optional! Default 0. the minimum max_delta_step for parameter testing of xgboost model. Maximum delta step we allow each leaf output to be. If the value is set to 0, it means there is no constraint. If it is set to a positive value, it can help making the update step more conservative. Usually this parameter is not needed, but it might help in logistic regression when class is extremely imbalanced. Set it to value of 1-10 might help control the update.')
    parser.add_argument("--max_delta_step_max", type=int, default=5, help = 'Optional! Default 5. the maximum max_delta_step for parameter testing of xgboost model. Maximum delta step we allow each leaf output to be. If the value is set to 0, it means there is no constraint. If it is set to a positive value, it can help making the update step more conservative. Usually this parameter is not needed, but it might help in logistic regression when class is extremely imbalanced. Set it to value of 1-10 might help control the update.')
    parser.add_argument("--min_child_weight_min", type=float, default=0.1, help = 'Optional! Default 0.1. the minimum min_child_weight for parameter testing of xgboost model. Minimum sum of instance weight (hessian) needed in a child. If the tree partition step results in a leaf node with the sum of instance weight less than min_child_weight, then the building process will give up further partitioning. In linear regression task, this simply corresponds to minimum number of instances needed to be in each node. The larger min_child_weight is, the more conservative the algorithm will be.')
    parser.add_argument("--min_child_weight_max", type=float, default=15, help = 'Optional! Default 15. the maximum min_child_weight for parameter testing of xgboost model. Minimum sum of instance weight (hessian) needed in a child. If the tree partition step results in a leaf node with the sum of instance weight less than min_child_weight, then the building process will give up further partitioning. In linear regression task, this simply corresponds to minimum number of instances needed to be in each node. The larger min_child_weight is, the more conservative the algorithm will be.')
    parser.add_argument("--num_round_min", type=int, default=40, help = 'Optional! Default 40. the minimum number of rounds for boosting for parameter testing of xgboost model.')
    parser.add_argument("--num_round_max", type=int, default=2000, help = 'Optional! Default 2000. the maximum number of rounds for boosting for parameter testing of xgboost model.')
    parser.add_argument("--model_config", type=str, help = 'Required!')
    parser.add_argument("--model_dir", type=str, default='Required! output the trained xgboost model parameters.')
    args = parser.parse_args()

    main(args)