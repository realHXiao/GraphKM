#!/usr/bin/python
# coding: utf-8

# Author: LE YUAN
# Date: 2020-07-23  Run in python 3.7
# This script is to combine the KM data from BRENDA and Sabio-RK database

import csv
import json


with open('../KM_sabio_smiles.json', 'r') as infile:
    sabio_name_smiles = json.load(infile)

with open('../KM_brenda_smiles.json', 'r') as infile:
    brenda_name_smiles = json.load(infile)

with open("./KM_sabio_clean_unisubstrate_2.tsv", "r", encoding='utf-8') as file1 :
    sabio_lines = file1.readlines()[1:]

with open("../KM_brenda_clean.tsv", "r", encoding='utf-8') as file2 :
    brenda_lines = file2.readlines()[1:]

KM_data = list()
KM_data_include_value = list()
Substrate_name = dict()
Substrate_smiles = dict()
entry_uniprot = dict()

for line in brenda_lines :
    # print(line)
    data = line.strip().split('\t')
    ECNumber = data[1]
    Substrate = data[2]
    EnzymeType = set(data[3].split('/'))
    Organism =data[4]
    Value = data[5]
    Unit = data[6]
    Reference = data[7]
    smiles = brenda_name_smiles[Substrate]
    # print(smiles)
    if smiles is not None :
        # print(smiles)
        Substrate_name[Substrate.lower()] = Substrate
        Substrate_smiles[Substrate.lower()+'&smiles'] = smiles

        KM_data_include_value.append([ECNumber, Substrate.lower(), EnzymeType, Organism, Value, Unit, Reference])
#         print(KM_data_include_value)
        KM_data.append([ECNumber, Substrate.lower(), EnzymeType, Organism])

for line in sabio_lines :
    # print(line)
    data = line.strip().split('\t')
    ECNumber = data[1]
    Substrate = data[2]
    EnzymeType = set(data[3].split('/'))
    Organism =data[5]
    UniprotID = data[6]
    Value = data[7]
    Unit = data[8]
    Reference = data[9]
    smiles = sabio_name_smiles[Substrate]
    # print(smiles)
    if smiles is not None :
        # print(smiles)
        Substrate_name[Substrate.lower()] = Substrate
        Substrate_smiles[Substrate.lower()+'&smiles'] = smiles
        entry_uniprot[ECNumber + Substrate.lower() + Organism] = UniprotID

        KM_data_include_value.append([ECNumber, Substrate.lower(), EnzymeType, Organism, Value, Unit, Reference])
#         print(KM_data_include_value)
        KM_data.append([ECNumber, Substrate.lower(), EnzymeType, Organism])

print(len(KM_data))  # 49392


new_lines = list()
for line in KM_data :
    if line not in new_lines :
        new_lines.append(line)

print(len(new_lines))  # 48659 included all elements, 46165 included all except for KM value and unit,  44166 substrate lower and upper

i = 0
clean_KM = list()
for new_line in new_lines :
    # print(new_line)
    i += 1
    print(i)
    value_unit = dict()
    KM_values = list()
    for line in KM_data_include_value :
        if line[:-3] == new_line :
            value = line[-3]
            value_unit[str(float(value))] = line[-2]
            # print(type(value))  # <class 'str'>
            KM_values.append(float(value))
            reference = line[-1]
    # print(value_unit)
    # print(KM_values)
    max_value = max(KM_values) # choose the maximum one for duplication KM value under the same entry as the data what we use
    unit = value_unit[str(max_value)]
    # print(max_value)
    # print(unit)
    Substrate = Substrate_name[new_line[1]]
    Smiles = Substrate_smiles[new_line[1]+'&smiles']
    try :
        UniprotID = entry_uniprot[new_line[0]+new_line[1]+new_line[3]]
    except :
        UniprotID = ''
    new_line[1] = Substrate
    new_line[2] = '/'.join(new_line[2])
    new_line = new_line + [Smiles, UniprotID, str(max_value), unit, reference]
    # print(new_line)
    # new_line.append(str(max_value))
    # new_line.append(unit)
    if new_line[-2] == 'mM' or 'mg/ml' :
        clean_KM.append(new_line)

# print(clean_KM)
print(len(clean_KM))  # 

with open("./KM_combination_0814.tsv", "w") as outfile :
    records = ['ECNumber', 'Substrate', 'EnzymeType', 'Organism', 'Smiles', 'UniprotID', 'Value', 'Unit', 'Reference']
    outfile.write('\t'.join(records) + '\n')
    for line in clean_KM :
        outfile.write('\t'.join(line) + '\n')

# The above is to eliminate the dupliaction entries by Substrate name
print('-----------------------------------------------')

# This is to eliminate the duplication entries by Smiles
with open("./KM_combination_0814.tsv", "r", encoding='utf-8') as infile :
    KM_lines = infile.readlines()[1:]

KM_data = list()
KM_data_include_value = list()
Substrate_name = dict()
entry_uniprot = dict()

for line in KM_lines :
    # print(line)
    data = line.strip().split('\t')
    ECNumber = data[0]
    Substrate = data[1]
    EnzymeType = set(data[2].split('/'))
    Organism =data[3]
    Smiles = data[4]
    UniprotID = data[5]
    Value = data[6]
    Unit = data[7]
    Reference = data[8]
    Substrate_name[Smiles] = Substrate
    entry_uniprot[ECNumber + Smiles + Organism] = UniprotID

    KM_data_include_value.append([ECNumber, EnzymeType, Organism, Smiles, Value, Unit, Reference])
    KM_data.append([ECNumber, EnzymeType, Organism, Smiles])

print(len(KM_data))  # 44166


new_lines = list()
for line in KM_data :
    if line not in new_lines :
        new_lines.append(line)

print(len(new_lines))  # 43495 included all elements, 41558 included all except for KM value and unit

i = 0
clean_KM = list()
for new_line in new_lines :
    # print(new_line)
    i += 1
    print(i)
    value_unit = dict()
    KM_values = list()
    for line in KM_data_include_value :
        if line[:-3] == new_line :
            value = line[-3]
            value_unit[str(float(value))] = line[-2]
            # print(type(value))  # <class 'str'>
            KM_values.append(float(value))
            reference = line[-1]
    # print(value_unit)
    # print(KM_values)
    max_value = max(KM_values) # choose the maximum one for duplication KM value under the same entry as the data what we use
    unit = value_unit[str(max_value)]
    # print(max_value)
    # print(unit)
    Substrate = Substrate_name[new_line[3]]
    try :
        UniprotID = entry_uniprot[new_line[0]+new_line[3]+new_line[2]]
    except :
        UniprotID = ''
    new_line[1] = '/'.join(new_line[1])
    new_line = new_line + [Substrate, UniprotID, str(max_value), unit, reference]
    # print(new_line)
    # new_line.append(str(max_value))
    # new_line.append(unit)
    if new_line[-2] == 'mM' or 'mg/ml':
        clean_KM.append(new_line)

# print(clean_KM)
print(len(clean_KM))  # 41558, in which 13454 entries have UniprotID

with open("./KM_combination_0815.tsv", "w") as outfile :
    records = ['ECNumber', 'EnzymeType', 'Organism', 'Smiles', 'Substrate', 'UniprotID', 'Value', 'Unit', 'Reference']
    outfile.write('\t'.join(records) + '\n')
    for line in clean_KM :
        outfile.write('\t'.join(line) + '\n')