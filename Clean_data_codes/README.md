# Usage
Create an environment with python2.7. And download the latest [enzyme.dat](https://ftp.expasy.org/databases/enzyme/) file. 
```
conda create -n python2.7 python=2.7
conda activate python2.7
conda install SOAPpy
```
Run data cleaning codes in the rank: 
```
# clean data from BRENDA
python brenda_retrieve.py  # under python2.7. 
python brenda_download.py
python brenda_kcat_preprocess.py
python brenda_kcat_clean.py
python brenda_sequence_organism.py
python brenda_get_smiles.py  # needs pubchempy package.

# clean data from SABIO-RK
python sabio_download.py
python sabio_kcat_unisubstrate.py
python sabio_kcat_clean_unisubstrate.py
python sabio_kcat_clean.py
python sabio_kcat_unisubstrate_mutant.py
python sabio_get_smiles.py

# combine the data from BRENDA and SABIO-RK
python combination_brenda_sabio.py
python combination_database_data.py
```
