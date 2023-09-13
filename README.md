# GraphKM
![image](https://github.com/realHXiao/GraphKM/assets/71002556/cb30a0a9-3c48-432c-b1b9-173201875144)


## Introduction
The GraphKM toolbox is a Python package for prediction of KMs. 
## Requirements
Assuming that you use [Miniconda](https://docs.conda.io/en/latest/miniconda.html) or [Anaconda](https://www.anaconda.com/). In a terminal execute: 
```
conda env create -n GraphKM python=3.8
conda activate GraphKM
```
 Requirement packages: 
```
paddlehelix==1.0.1
pgl==2.2.4
paddlepaddle-gpu==2.3.2
matplotlib
scikit-learn
rdkit
PubChemPy
xgboost==1.7.5
hyperopt==0.2.7
Unirep50
```
Note: ``paddlepaddle-gpu==2.3.2`` is installed by command line ``conda install paddlepaddle-gpu==2.3.2 cudatoolkit=11.2 -c https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/Paddle/ -c conda-forge``. 

Please reference [this github site](https://github.com/EngqvistLab/UniRep50) for Unirep50 installation. 

## Input files
Before data preprocessing, a json file and a csv file should be ready. 

+ The json file is generated by codes in the `KM_data_clean/` dictionary, and is composed by dictionaries, like:
    ```
    [
        {
            "ECNumber": "2.1.1.22",
            "Organism": "Rattus norvegicus",
            "Smiles": "C[S+](CCC(C(=O)[O-])N)CC1C(C(C(O1)N2C=NC3=C(N=CN=C32)N)O)O",
            "Substrate": "S-Adenosyl-L-methionine",
            "Sequence": "MQRRRRAPPASQPAQDSGHSEEVEVQFSAGRLGSAAPAGPPVRGTAEDEERLEREHFWKVINAFRYYGTSMHERVNRTERQFRSLPDNQQKLLPQFPLHLDKIRKCVDHNQEILLTIVNDCIHMFENKEYGEDANGKIMPASTFDMDKLKSTLKQFVRDWSGTGKAERDACYKPIIKEIIKNFPKERWDPSKVNILVPGAGLGRLAWEIAMLGYACQGNEWSFFMLFSSNFVLNRCSEVDKYKLYPWIHQFSNNRRSADQIRPIFFPDVDPHSLPPGSNFSMTAGDFQEIYSECNTWDCIATCFFIDTAHNVIDYIDTIWRILKPGGIWINLGPLLYHFENLANELSIELSYEDIKNVVLQYGFQLEVEKESVLSTYTVNDLSMMKYYYECVLFVVRKPQ",
            "Type": "wildtype",
            "Value": "0.042",
            "Unit": "mM"
        }
    ]
    ```

+ The csv file is generated from Unirep50 (see `Generate Unirep50 vectors file` section for usage process) and contains Unirep50 vectors. 

    The sequences you wish to convert to UniRep embeddings need to be collected in a single FASTA file. And run following codes: 
    ```python
    >>> from unirep.run_inference import BatchInference
    >>> inf_obj = BatchInference(batch_size=256)
    >>> df = inf_obj.run_inference(filepath='my_sequences.fasta')
    >>> df.to_csv('my_sequences_embeddings.tsv', sep='\t')
    ```
## Train
### Preprocess
```
python data_preprocess.py -i my_data.json -l KM -i_unirep my_protein_sequences_embeddings.csv -o my_dataset.npz
```
### Training
The training needs big memory if you use GPU for acceleration. Suggestion that the memory of your GPU is 24 GB. 

```
python train.py -d path_to/my_dataset.npz --model_config path_to/gin_config.json -l KM -- model_dir path_to/ --results_dir path_to/

python train_xgb.py -i path_to/my_data.json -l KM -i_unirep path_to/my_protein_sequences_embeddings.csv -m path_to/best_model_gin_-1_lr0.0005.pdparams --model_config path_to/gin_config.json
```
## Prediction
The input for prediction.py:
+ If you want to predict KM values of different seuqences corresponding to different substrate SMILES codes, use csv file as input. The format of csv file please refer to the example.csv file. The commond line example for prediction:

    ```
    python prediction.py -c --csv_file example.csv -l KM -i_unirep example.tsv -m path_to/best_model_gin_-1_lr0.0005.pdparams --model_config gin_config.json -xgb path_to/gin_xgboost_model.dat
    ```
+ If you want to predict KM values of different seuqences corresponding to one type substrate SMILES codes, use FASTA file as input. 

    commond line example for prediction:
    ```
    python prediction.py -l KM -f --fasta_file example.fasta -i_unirep my_sequences_embeddings.tsv -S substrate.txt -m path_to/best_model_gin_-1_lr0.0005.pdparams --model_config path_to/gin_config.json -xgb path_to/gin_xgboost_model.dat
    ```
## tip
Enter `-h` tag for more helps. 
```
python data_preprocess.py -h
python train.py -h
python train_xgb.py -h
python prediction.py -h
```
