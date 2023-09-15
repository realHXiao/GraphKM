#!/usr/bin/python
# coding: utf-8

# Author: Xiao He
# Date: 2023-08-05  Run in python 3.7
# This script is to clean KM data extracted from BRENDA database

import csv

with open("./KM_brenda.tsv", "r", encoding='utf-8') as file :
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
    Organism =data[5]
    Value = data[6]
    Unit = data[7]
    KM_data_include_value.append([Type, ECNumber, Substrate, EnzymeType, Organism, Value, Unit])
    KM_data.append([Type, ECNumber, Substrate, EnzymeType, Organism])

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
        if line[:5] == new_line :
            value = line[5]
            value_unit[str(float(value))] = line[6]
            # print(type(value))  # <class 'str'>
            KM_values.append(float(value))
    # print(value_unit)
    # print(KM_values)
    max_value = max(KM_values) # choose the maximum one for duplication KM value under the same entry as the data what we use
    unit = value_unit[str(max_value)]
    # print(max_value)
    # print(unit)

    new_line.append(str(max_value))
    new_line.append(unit)
    new_line.append(reference)
    if new_line[6] == 'mM' :
        clean_KM.append(new_line)

# print(clean_KM)
print(len(clean_KM))  # 52390


with open("./KM_brenda_clean.tsv", "w") as outfile :
    records = ['Type', 'ECNumber', 'Substrate', 'EnzymeType', 'Organism', 'Value', 'Unit']
    outfile.write('\t'.join(records) + '\n')
    for line in clean_KM :
        outfile.write('\t'.join(line) + '\n')

