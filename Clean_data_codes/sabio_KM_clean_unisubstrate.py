#!/usr/bin/python
# coding: utf-8

# Author: Xiao He
# Date: 2023-08-08  Run in python 3.7


import csv

with open("./amendment/KM_sabio_4_unisubstrate.tsv", "r", encoding='utf-8') as file :
    lines = file.readlines()[1:]


KM_data = list()
KM_data_include_value = list()
for line in lines :
    # print(line)
    data = line.strip().split('\t')
    Type = data[1]
    ECNumber = data[2]
    Substrate = data[3]
    EnzymeType = data[4]
    PubMedID = data[5]
    Organism =data[6]
    UniprotID = data[7]
    Value = data[8]
    Unit = data[9]
    KM_data_include_value.append([Type, ECNumber, Substrate, EnzymeType, PubMedID, Organism, UniprotID, Value, Unit])
    KM_data.append([Type, ECNumber, Substrate, EnzymeType, PubMedID, Organism, UniprotID])

print(len(KM_data))  


new_lines = list()
for line in KM_data :
    if line not in new_lines :
        new_lines.append(line)

print(len(new_lines))  

i = 0
clean_KM = list()
for new_line in new_lines :
    # print(new_line)
    i += 1
    print(i)
    value_unit = dict()
    KM_values = list()
    for line in KM_data_include_value :
        if line[:7] == new_line and line[7] != '': 
            value = line[7]
            value_unit[str(float(value))] = line[8]
            # print(type(value))  # <class 'str'>
            KM_values.append(float(value))
    # print(value_unit)
    print(KM_values)
    if len(KM_values) > 0: 
        max_value = max(KM_values) # choose the maximum one for duplication KM value under the same entry as the data what we use
        unit = value_unit[str(max_value)]
    # print(max_value)
    # print(unit) # ['M', 'M^2', 'mg/ml', '-', 'katal*g^(-1)', 'l*g^(-1)*s^(-1)', 'mol', 'mol/mol', 'mol*s^(-1)*g^(-1)']

    if unit == 'M' :
        max_value = max_value * 1000
        unit = 'mM'
        new_line.append(str(max_value))
        new_line.append(unit)
    if unit == 'M^2' :
        max_value = max_value ** (1/2) * 1000
        unit = 'mM'
        new_line.append(str(max_value))
        new_line.append(unit)
    if unit == 'mg/ml' : 
        new_line.append(str(max_value))
        new_line.append(unit)
    if new_line[-1] == 'SABIO-RK' and (new_line[-2] == 'mM' or 'mg/ml') :
        clean_KM.append(new_line)

print(len(clean_KM))  


with open("./amendment/KM_sabio_clean_unisubstrate.tsv", "w") as outfile :
    records = ['Type', 'ECNumber', 'Substrate', 'EnzymeType', 'PubMedID', 'Organism', 'UniprotID', 'Value', 'Unit']
    outfile.write('\t'.join(records) + '\n')
    for line in clean_KM :
        outfile.write('\t'.join(line) + '\n')
