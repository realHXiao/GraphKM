#!/usr/bin/python
# coding: utf-8

# Author: LE YUAN
# Date: 2020-07-08  Run in python 3.7


import csv

with open("./KM_sabio_clean_unisubstrate.tsv", "r", encoding='utf-8') as file :
    lines = file.readlines()[1:]


KM_data = list()
KM_data_include_value = list()
n = 0
for line in lines :
    n += 1
    print(n)
    # print(line)
    data = line.strip().split('\t')
    print(data)
    Type = data[0]
    ECNumber = data[1]
    Substrate = set(data[2].split(';'))
    EnzymeType = data[3]
    PubMedID = data[4]
    Organism =data[5]
    UniprotID = data[6]
    Value = data[7]
    Unit = data[8]
    KM_data_include_value.append([Type, ECNumber, Substrate, EnzymeType, PubMedID, Organism, UniprotID, Value, Unit])
    KM_data.append([Type, ECNumber, Substrate, EnzymeType, PubMedID, Organism, UniprotID])

print(len(KM_data))  # 22416


new_lines = list()
for line in KM_data :
    if line not in new_lines :
        new_lines.append(line)

print(len(new_lines))  # 20344 included all elements, 16532 included all except for KM value and unit

i = 0
clean_KM = list()
for new_line in new_lines :
    # print(new_line)
    i += 1
    print(i)
    value_unit = dict()
    KM_values = list()
    for line in KM_data_include_value :
        if line[:7] == new_line :
            value = line[7]
            value_unit[str(float(value))] = line[8]
            # print(type(value))  # <class 'str'>
            KM_values.append(float(value))
    # print(value_unit)
    # print(KM_values)
    max_value = max(KM_values) # choose the maximum one for duplication KM value under the same entry as the data what we use
    unit = value_unit[str(max_value)]
    # print(max_value)
    # print(unit)

    if unit == 'mg/ml' or 'mM': 
        new_line.append(str(max_value))
        new_line.append(unit)
        reference = 'SABIO-RK'
        new_line.append(reference)
    new_line[2] = ';'.join(list(new_line[2]))
    new_line.append(str(max_value))
    new_line.append(unit)
    if new_line[-2] == 'mg/ml' or 'mM' :
        clean_KM.append(new_line)

# print(clean_KM)
print(len(clean_KM))  # 16484 after unifing the KM value unit to 'mM', in which 15114 has a specific Unipro ID


with open("./KM_sabio_clean.tsv", "w") as outfile :
    records = ['Type', 'ECNumber', 'Substrate', 'EnzymeType', 'PubMedID', 'Organism', 'UniprotID','Value', 'Unit', 'Reference']
    outfile.write('\t'.join(records) + '\n')
    for line in clean_KM :
        outfile.write('\t'.join(line) + '\n')

        