#!/usr/bin/python
# coding: utf-8

# Author: Xiao He
# Date: 2023-08-05


# E-mail in BRENDA:
email = 'youremail'
# Password in BRENDA:
password = 'yourpassword'


# #Construct BRENDA client:
import string
import hashlib
import os
import json
import time
import random
from SOAPpy import SOAPProxy ## for usage without WSDL file

endpointURL = "https://www.brenda-enzymes.org/soap/brenda_server.php"
client      = SOAPProxy(endpointURL)
password    = hashlib.sha256(password).hexdigest()
credentials = email + ',' + password


filenames = os.listdir('./KM_brenda')
# print(len(filenames))
i = 0

EC_organisms = dict()
# print(filenames)
filenames.sort()
for filename in filenames :
    EC = filename[2:-4]
    print(EC)
    if filename != '.DS_Store' :
        with open("./KM_brenda/%s" %(filename), 'r') as file :
            lines = file.readlines()
    organisms = list()
    for line in lines[1:] :
        data = line.strip().split('\t')
        organism = data[1]
        organisms.append(organism)
    organisms = list(set(organisms))
    print(organisms)
    organism_seqcounts = dict()
    for organism in organisms :
        print(organism)
        succes = 0
        #time.sleep(random.randint(0, 10))
        while succes < 10:
            try:
                parameters = credentials+","+"ecNumber*%s#organism*%s" %(EC, organism)
                sequence = client.getSequence(parameters)
                print(sequence)
                split_sequences = sequence.strip().split('#!') #noOfAminoAcids #!
                organism_seqcounts[organism] = len(split_sequences)
                succes = 10

            except:
                #Let the server cool of for a bit. If after 10 times it
                #still fails, the query is discarded:
                time.sleep(1)
                succes += 1

        # for seq in split_sequences :
        #     print(seq)
        #     print('--------------------------------')

        # print(len(split_sequences))

    EC_organisms[EC] = organism_seqcounts

print(EC_organisms)
with open('./brenda_EC_organims_try.json', 'w') as outfile:
    json.dump(EC_organisms, outfile, indent=4)
    
