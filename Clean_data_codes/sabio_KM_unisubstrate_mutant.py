#     else :
#     # if 'mutant' in desc or 'mutated' in desc:
#         print(desc)
#         mutant = re.findall('[A-Z]\d+[A-Z]', desc)  # re is of great use
#         if len(mutant) >=1 :
#             enzymeType = '/'.join(mutant)

#     if enzymeType :
#         enzymeType_entries.append(enzymeType)

# # print(enzymeType_entries)
# print(len(enzymeType_entries))  

import re

with open("./KM_sabio_clean_unisubstrate.tsv", "r", encoding='utf-8') as file :
    lines = file.readlines()[1:]

clean_mutant = list()
for line in lines :
    # print(line)
    data = line.strip().split('\t')
    Type = data[0]
    ECNumber = data[1]
    Substrate = data[2]
    EnzymeType = data[3]
    PubMedID = data[4]
    Organism =data[5]
    UniprotID = data[6]
    Value = data[7]
    Unit = data[8]
    Reference = data[9]
    if 'wildtype' in EnzymeType :
        enzymeType = 'wildtype'
    else :
    # if 'mutant' in EnzymeType or 'mutated' in EnzymeType:
        print(EnzymeType)
        mutant = re.findall('[A-Z]\d+[A-Z]', EnzymeType)  # re is of great use
        enzymeType = '/'.join(mutant)

    print(enzymeType)
    if enzymeType :
        clean_mutant.append([Type, ECNumber, Substrate, enzymeType, PubMedID, Organism, UniprotID, Value, Unit, Reference])


# print(enzymeType_entries)
print(len(clean_mutant))  # 17384


with open("./KM_sabio_clean_unisubstrate_2.tsv", "w") as outfile :
    records = ['Type', 'ECNumber', 'Substrate', 'EnzymeType', 'PubMedID', 'Organism', 'UniprotID', 'Value', 'Unit', 'Reference']
    outfile.write('\t'.join(records) + '\n')
    for line in clean_mutant :
        outfile.write('\t'.join(line) + '\n')