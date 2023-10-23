import os
import numpy as np
import json
import random
import logging
import argparse
from pahelix.utils.data_utils import load_npz_to_data_list
from src.data_gen import KMDataset  

import pgl
import paddle
import numpy as np
from src.model import KMModel, KMModelCriterion
from src.data_gen import DataCollateFunc
from sklearn.metrics import mean_squared_error,r2_score
import matplotlib.pyplot as plt

def shuffle_dataset(dataset, seed):
    np.random.seed(seed)
    np.random.shuffle(dataset)
    return dataset

def split_dataset(dataset, ratio):
    n = int(ratio * len(dataset))
    dataset_1, dataset_2 = dataset[:n], dataset[n:]
    return dataset_1, dataset_2

def train(model, criterion, optimizer, dataloader):
    model.train()
    list_loss = []
    train_pred, train_label = [], []
    for graphs, proteins_seq, labels in dataloader:
        graphs = graphs.tensor()
        proteins_seq = paddle.to_tensor(proteins_seq)
        labels = paddle.to_tensor(labels)
        
        preds = model(graphs, proteins_seq)
        loss = criterion(preds, labels)
        
        loss.backward()
        optimizer.step()
        optimizer.clear_grad()
        list_loss.append(loss.numpy())
        
        train_pred.append(preds.numpy())
        train_label.append(labels)
        
    train_pred = np.concatenate(train_pred, 0).flatten()
    train_label = np.concatenate(train_label, 0).flatten()
    
    mse = ((train_label - train_pred) ** 2).mean(axis=0)  
    rmse = np.sqrt(mean_squared_error(train_label, train_pred))
    r2 = r2_score(train_label, train_pred)
    return np.mean(list_loss), rmse, r2, mse 

def evaluate(model, dataloader, prior_best_rmse):
    model.eval()
    total_pred, total_label = [], []
    for graphs, proteins_seq, labels in dataloader:
        graphs = graphs.tensor()
        proteins_seq = paddle.to_tensor(proteins_seq)
        
        preds = model(graphs, proteins_seq)
        total_pred.append(preds.numpy())
        total_label.append(labels)

    total_pred = np.concatenate(total_pred, 0).flatten()
    total_label = np.concatenate(total_label, 0).flatten()

    mse = ((total_label - total_pred) ** 2).mean(axis=0)  
    rmse = np.sqrt(mean_squared_error(total_label, total_pred))
    r2 = r2_score(total_label, total_pred)

    return rmse, r2, mse, total_label, total_pred

def main(args): 
    model_config = json.load(open(args.model_config, 'r'))

    logging.info('Load data ...')
    dataset = load_npz_to_data_list(os.path.join(args.dataset))
    dataset = shuffle_dataset(dataset, args.shuffle)
    train_data, test_data = split_dataset(dataset, args.split)

    max_protein_len = model_config["protein"]["max_protein_len"]
    train_dataset = KMDataset(train_data, max_protein_len=max_protein_len)
    test_dataset = KMDataset(test_data, max_protein_len=max_protein_len)
    print('training dataset number:', len(train_dataset))
    print('testing dataset number:', len(test_dataset))

    label_name = args.label_type
    collate_fn = DataCollateFunc(
    model_config['compound']['atom_names'],
    model_config['compound']['bond_names'],
    is_inference=False,
    label_name=label_name)

    train_dataloader = train_dataset.get_data_loader(
        batch_size=args.batch_size,
        num_workers=args.num_workers,
        collate_fn=collate_fn)

    test_dataloader = test_dataset.get_data_loader(
            batch_size=args.batch_size,
            num_workers=1,
            shuffle=False,
            collate_fn=collate_fn)

    logging.info("Data loaded.")

    model = KMModel(model_config)
    criterion = KMModelCriterion()

    if args.lr_CosineAnnealingDecay is False: 
        lr = args.lr
        optimizer = paddle.optimizer.Adam(
            learning_rate=lr,
            parameters=model.parameters())
    else: 
        lr_min = args.lr_min
        lr = args.lr_max
        l_r = paddle.optimizer.lr.CosineAnnealingDecay(learning_rate=lr, eta_min=lr_min, T_max=args.lr_cycle, verbose=True)
        optimizer = paddle.optimizer.Adam(
            learning_rate=l_r,
            parameters=model.parameters())

    os.makedirs(args.model_dir, exist_ok=True)

    best_rmse, best_r2, best_ep = np.inf, 0, 0  # np.inf: +∞
    best_model = os.path.join(args.model_dir, 'best_model_{}_{}_lr{}.pdparams'.format(model_config["compound"]["gnn_type"], model_config["protein"]["max_protein_len"], lr))
    metric = None
    for epoch_id in range(1, args.max_epoch + 1):
        print('========== Epoch {} =========='.format(epoch_id))
        train_loss, train_rmse, train_r2, train_mse = train(model, criterion, optimizer, train_dataloader)
        print('Epoch: {}, Train loss: {}'.format(epoch_id, train_loss))
        test_rmse, test_r2, test_mse, total_label, total_pred = evaluate(model, test_dataloader, best_rmse)
        metric = 'Epoch: {}, Train_RMSE: {}, Train_R2: {}, Train_MSE: {}, Test_RMSE: {}, Test_R2: {}, Test_MSE: {}'.format(
                epoch_id, train_rmse, train_r2, train_mse, test_rmse, test_r2, test_mse)
        
        # save metrices generated in training process. 
        with open(os.path.join(args.results_dir, '{}_training_{}_lr{}.txt'.format(model_config["compound"]["gnn_type"], model_config["protein"]["max_protein_len"], lr)), 'a') as f0:
            f0.write(str(metric) + '\n')
        
        if test_rmse < best_rmse:
            best_rmse, best_r2, best_mse, best_ep = test_rmse, test_r2, test_mse, epoch_id  
            paddle.save(model.state_dict(), best_model)
            metric = 'Epoch: {}, Best RMSE: {}, Best R2: {}, Best MSE: {}'.format(
                best_ep, best_rmse, best_r2, best_mse)
            
            # save better metrics. 
            with open(os.path.join(args.results_dir, '{}_results_{}_lr{}.txt'.format(model_config["compound"]["gnn_type"], model_config["protein"]["max_protein_len"], lr)), 'a') as f:
                f.write(str(metric) + '\n')
            
            # save KM values and predicted KM values. 
            with open(os.path.join(args.results_dir, '{}_{}_lr{}_test_dataset_label_pred.txt'.format(model_config["compound"]["gnn_type"], model_config["protein"]["max_protein_len"], lr)), 'w') as fi:
                fi.write('Epoch: {}'.format(best_ep) + '\n')
                for m in range(0, len(total_label)): 
                    fi.write('label: ' + str(total_label[m]) + ' pred: ' + str(total_pred[m]) + '\n')
            fi.close()
            
            # plot correlation figure between true KM values and predicted KM values. 
            folder = os.path.join(args.results_dir, '{}_{}_lr{}_plt_image'.format(model_config["compound"]["gnn_type"], model_config["protein"]["max_protein_len"], lr))
            if not os.path.exists(folder):
                os.makedirs(folder)
            plt.figure(figsize=(6, 6), dpi=300)
            plt.scatter(total_label, total_pred, alpha=0.1)  # scatter transparency is 90%
            plt.xlabel("label", fontdict={'size': 16})
            plt.ylabel("pred", fontdict={'size': 16})
            plt.savefig(os.path.join(args.results_dir, '{}_{}_lr{}_plt_image/Epoch: {}.jpg'.format(model_config["compound"]["gnn_type"], model_config["protein"]["max_protein_len"], lr, best_ep)))
        else:
            print('No improvement in epoch {}'.format(epoch_id))
        
        if args.lr_CosineAnnealingDecay is True:
            l_r.step()
        
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-l", '--label_type', type=str, default='KM', help = 'Optional! Default KM. choose KM or kcat.')
    parser.add_argument("--split", type=float, default=0.8, help = 'Optional! Default 0.8. split whole dataset into train_dataset and test_dataset. if set 0.8, train_dataset = 80/100 whole dataset.')
    parser.add_argument("--shuffle", type=int, default=1234, help = 'Optional! Default 1234. shuffle dataset.')
    parser.add_argument('--lr_CosineAnnealingDecay', action="store_true", default=False, help = 'Optional! Default False (Disable). ')
    parser.add_argument("--batch_size", type=int, default=128, help = 'Optional! Default 128.')
    parser.add_argument("--num_workers", type=int, default=4, help = 'Optional! Default 4.')
    parser.add_argument("--max_epoch", type=int, default=200, help = 'Optional! Default 200.')
    parser.add_argument("--lr", type=float, default=0.0005, help = 'Optional! Default 0.0005.')
    parser.add_argument("--lr_min", type=float, default=0, help = 'Optional! Default 0. Enable after --lr_CosineAnnealingDecay is specified.')
    parser.add_argument("--lr_max", type=float, default=0.0005, help = 'Optional! Default 0.0005. Enable after --lr_CosineAnnealingDecay is specified.')
    parser.add_argument("--lr_cycle", type=int, default=10, help = 'Optional! Default ')
    parser.add_argument("-d", "--dataset", type=str, help = 'Required!')
    parser.add_argument("--model_config", type=str, help = 'Required!')
    parser.add_argument("--model_dir", type=str, help = 'Required!')
    parser.add_argument("--results_dir", type=str, help = 'Required!')
    args = parser.parse_args()
    main(args)
