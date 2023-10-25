#!/usr/bin/python
# coding: utf-8

# Author: Xiao He
# Date: 2023-08-08

# This python script is to obtain protein sequence for each KM entries

import os
import re
import json
import requests
import time
from urllib import request
from zeep import Client
import hashlib

# This function is to obtain the protein sequence according to the protein id from Uniprot API
# https://www.uniprot.org/uniprot/A0A1D8PIP5.fasta
# https://www.uniprot.org/help/api_idmapping
def uniprot_sequence(id) :
    url = "https://www.uniprot.org/uniprot/%s.fasta" % id
    IdSeq = dict()
    succes = 0
    while succes < 10:
        try :
            data = request.urlopen(url)
            respdata = data.read().decode("utf-8").strip()
            IdSeq[id] =  "".join(respdata.split("\n")[1:])
            succes = 10
            
        except :
            #print(id, "can not find from uniprot!")
            #Let the server cool of for a bit. If after 10 times it
 
            #still fails, the query is discarded:
    
            time.sleep(1)
            succes += 1
            #IdSeq[id] = None
    print(IdSeq[id])
    return IdSeq[id]    
    
def uniprotID_entry() :
    with open("./KM_combination_0815.tsv", "r", encoding='utf-8') as file :
        combination_lines = file.readlines()[1:]

    uniprotID_list = list()
    uniprotID_seq = dict()
    uniprotID_noseq = list()

    i=0
    for line in combination_lines :
        data = line.strip().split('\t')
        uniprotID = data[5]

        if uniprotID :
        #     seq = uniprot_sequence('P49384')
            if ' ' in uniprotID :
                # i += 1  # 561
                # print(i)
                # print(uniprotID.split(' '))
                uniprotID_list += uniprotID.split(' ')
            else :
                # print(uniprotID)
                uniprotID_list.append(uniprotID)
    uniprotID_unique = list(set(uniprotID_list))
    print(len(uniprotID_unique)) # 1776

    for uniprotID in uniprotID_unique :
        i += 1
        print(i)
        sequence = uniprot_sequence(uniprotID)
        if sequence :
            uniprotID_seq[uniprotID] = sequence
        else :
            uniprotID_noseq.append(uniprotID)


    print(len(uniprotID_seq))  # 1755
    print(len(uniprotID_noseq))  # 21
    print(uniprotID_noseq)
    # ['O05783', 'P56967', 'P64192', 'Q10711', 'P96420', 'O06201', 'P0A4X4', 'Q07636', 'A9ALT1', 'P50224', 'A1VCV2', 'G4VAM6', 'O60344', 'Q47741', 'P0A5Y6', 'F1P2T2', 'A8J9S7', 'F1PGK8', 'C3WBR4', 'P0A5R0', 'P65163', 'P11163', 'P00892', 'A9AET5', 'P96240', 'P0A4Z2', 'P0A4X6', 'P96807', 'Q90Y38', 'D5A598', 'O53447', 'Q4QRI0', 'C3WBR3', 'O52310', 'Q02469', 'P51698', 'D4ZTT4', 'A4VVM9', 'P33064', 'P0C5C1']
    # check one by one

    with open('./uniprotID_entry.json', 'w') as outfile :
        json.dump(uniprotID_seq, outfile, indent=4)
    with open('./uniprotID_noseq.json', 'w') as outfile :
        json.dump(uniprotID_noseq, outfile, indent=4)

def uniprotID_noseq() :
    with open('./uniprotID_entry.json', 'r') as infile :
        uniprotID_seq = json.load(infile)

    print(len(uniprotID_seq))

    uniprotID_noseq = {'O05783':'P9WIQ2', 'P56967':'F2MMP0', 'P64192':'P9WN68', 'Q10711':'P0DPE1', 'P96420':'P9WQB2', 'O06201':'P9WMK8', 
    'P0A4X4':'P9WQ86', 'Q07636':'P0DOB5', 'P50224':'P0DMM9', 'O60344':'P0DPD6', 'Q47741':'F2MMN9', 'P0A5Y6':'P9WGR0', 
    'P0A5R0':'P9WIL4', 'P65163':'P9WKJ0', 'P11163':'P0CX77', 'P00892':'P0DP89', 'P96240':'P9WIC2', 'P0A4Z2':'P9WPY2', 
    'P0A4X6':'P9WQ80', 'P96807':'P9WNP2', 'O53447':'P9WN20', 'O52310':'P0CL72', 'Q02469':'P0C278', 'P51698':'A0A1L5BTC1', 
    'P33064':'P0DOQ5', 'P0C5C1':'P9WKD2'}
    # 'A9ALT1', 'A1VCV2', 'G4VAM6', 'F1P2T2', 'A8J9S7', 'F1PGK8', 'C3WBR4', 'A9AET5', 'Q90Y38', 'D5A598', 'Q4QRI0', 'C3WBR3', 'D4ZTT4', 'A4VVM9'  On August 9, 2023 this entry was made redundant.

    for uniprotID, mappedID in uniprotID_noseq.items() :
        sequence = uniprot_sequence(mappedID)
        print(uniprotID)
        print(sequence)
        if sequence :
            uniprotID_seq[uniprotID] = sequence
        else :
            print('No sequence found!---------------------------')

    print(len(uniprotID_seq))  # 4570   

    with open('./uniprotID_entry_all.json', 'w') as outfile :
        json.dump(uniprotID_seq, outfile, indent=4)

# You can try to retrieve sequences from uniprot using rest interface.
# Example: (ec: 1.1.1.1 , organisms: Homo sapiens)
# http://www.uniprot.org/uniprot/?query=ec:1.1.1.1+AND+organism:"Homo sapiens"&format=fasta
# full information abut syntax you can find here: http://www.uniprot.org/help/programmatic_access
def seq_by_ec_organism(ec, organism) :
    IdSeq = dict()
    # https://www.biostars.org/p/356687/
    params = {"query": "ec:%s AND organism:%s AND reviewed:yes" % (ec, organism), "format": "fasta"}
    response = requests.get("http://www.uniprot.org/uniprot/", params=params)
    # print(type(response.text)) # <class 'str'>

    try :
        respdata = response.text
        # print(respdata)
        sequence = list()
        seq = dict()
        i = 0
        for line in respdata.split('\n') :
            if line.startswith('>') :
                name=line
                seq[name] = ''
            else :
                seq[name] += line.replace('\n', '').strip()
        IdSeq[ec+'&'+organism] =  list(seq.values())

    except :
        print(ec+'&'+organism, "can not find from uniprot!")
        IdSeq[ec+'&'+organism] = None

    print(IdSeq[ec+'&'+organism])
    return IdSeq[ec+'&'+organism]

# Run in python 2.7
def seq_by_brenda(ec, organism) :
    #New method using Python 3 because using Python 2 method provided by BRENDA could just run less than 10 hits as above
    # E-mail in BRENDA:
    email = 'youremail'
    # Password in BRENDA:
    password = 'yourpassword'

    wsdl = "https://www.brenda-enzymes.org/soap/brenda_zeep.wsdl"
    password    = hashlib.sha256(password.encode("utf-8")).hexdigest()
    client = Client(wsdl)

    parameters = ( email,password,"ecNumber*%s" % ec,"organism*%s" % organism, "sequence*", "noOfAminoAcids*", "firstAccessionCode*", "source*Swiss-Prot", "id*" ) # *Swiss-Prot
    entries = client.service.getSequence(*parameters)

    sequences = list()

    if entries :
        for entry in entries :
            sequences.append(entry['sequence'])

    print(sequences)
    print(len(sequences))
    return sequences

def nouniprotID_entry_uniprot() :
    with open("./KM_combination_0815.tsv", "r", encoding='utf-8') as file :
        combination_lines = file.readlines()[1:]

    IdSeq = dict()
    entries = list()
    i=0
    for line in combination_lines :
        data = line.strip().split('\t')
        ec = data[0]
        organism = data[2]
        uniprotID = data[5]

        if not uniprotID :
            entries.append((ec,organism))

    entries_unique = set(entries)

    for entry in list(entries_unique) :
        # print(entry)
        ec, organism = entry[0], entry[1]
        i += 1
        print('This is', str(i)+'------------')
        
        succes = 0
        while succes < 10:
            try:
                IdSeq[ec+'&'+organism] = seq_by_ec_organism(ec, organism)
                succes = 10
                # print(len(IdSeq)
            except:
                time.sleep(1)
                succes += 1

    with open('./nouniprotID_entry_all.json', 'w') as outfile :
        json.dump(IdSeq, outfile, indent=4)

# Run in python 2.7
def nouniprotID_entry_brenda() :
    with open("./KM_combination_0815.tsv", "r") as file :
        combination_lines = file.readlines()[1:]

    IdSeq = dict()
    entries = list()
    i=0
    for line in combination_lines :
        data = line.strip().split('\t')
        ec = data[0]
        organism = data[2]
        uniprotID = data[5]

        if not uniprotID :
            entries.append((ec,organism))

    entries_unique = set(entries)

    for entry in list(entries_unique) :
        ec, organism = entry[0], entry[1]
        i += 1
        print('This is', str(i)+'------------')
        succes = 0
        while succes < 10:
            try:
                IdSeq[ec+'&'+organism] = seq_by_brenda(ec,organism)
                succes = 10
            except:
                time.sleep(1)
                succes += 1


    with open('./nouniprotID_entry_brenda.json', 'w') as outfile :
        json.dump(IdSeq, outfile, indent=4)

def combine_sequence() :
    with open('./uniprotID_entry_all.json', 'r') as file1:
        uniprot_file1 = json.load(file1)

    with open('./nouniprotID_entry_all.json', 'r') as file2:  # By Uniprot API
        nouniprot_file2 = json.load(file2)

    with open('./nouniprotID_entry_brenda.json', 'r') as file3:  # By BRENDA API
        nouniprot_file3 = json.load(file3)

    with open("./KM_combination_0815.tsv", "r", encoding='utf-8') as file4 :
        KM_lines = file4.readlines()[1:]

    i = 0
    j = 0
    n = 0
    entries = list()
    for line in KM_lines :
        data = line.strip().split('\t')
        ECNumber, EnzymeType, Organism, Smiles = data[0], data[1], data[2], data[3]
        Substrate, UniprotID, Value, Unit = data[4], data[5], data[6], data[7]

        RetrievedSeq = ''
        entry = dict()
        if UniprotID :
            try :  # because a few (maybe four) UniprotIDs have no ID as the key 
                if ' ' not in UniprotID :
                    RetrievedSeq = [uniprot_file1[UniprotID]]
                else :
                    RetrievedSeq1 = [uniprot_file1[UniprotID.split(' ')[0]]]
                    RetrievedSeq2 = [uniprot_file1[UniprotID.split(' ')[1]]]
                    if RetrievedSeq1 == RetrievedSeq2 :
                        RetrievedSeq = RetrievedSeq1
            except :
                continue

        else :
            if nouniprot_file2[ECNumber+'&'+Organism] :
                if len(nouniprot_file2[ECNumber+'&'+Organism]) == 1 :
                    RetrievedSeq = nouniprot_file2[ECNumber+'&'+Organism]
                else :
                    RetrievedSeq = ''

        try:  # local variable 'RetrievedSeq' referenced before assignment
            if len(RetrievedSeq) == 1 and EnzymeType == 'wildtype': 
                sequence = RetrievedSeq
                i += 1
                entry = {
                    'ECNumber': ECNumber,
                    'Organism': Organism,
                    'Smiles': Smiles,
                    'Substrate': Substrate,
                    'Sequence': sequence[0],
                    'Type': 'wildtype',
                    'Value': Value,
                    'Unit': Unit,
                }

                entries.append(entry)

            if len(RetrievedSeq) == 1 and EnzymeType != 'wildtype':
                sequence = RetrievedSeq[0]

                mutantSites = EnzymeType.split('/')
                # print(mutantSites)

                mutant1_1 = [mutantSite[1:-1] for mutantSite in mutantSites]
                mutant1_2 = [mutantSite for mutantSite in mutantSites]
                mutant1 = [mutant1_1, mutant1_2]
                mutant2 = set(mutant1[0])
                if len(mutant1[0]) != len(mutant2) :
                    print(mutant1)
                    n += 1
                    print(str(n) + '---------------------------') 

                mutatedSeq = sequence
                for mutantSite in mutantSites :
                    if mutatedSeq[int(mutantSite[1:-1])-1] == mutantSite[0] :
                        # pass
                        mutatedSeq = list(mutatedSeq)
                        mutatedSeq[int(mutantSite[1:-1])-1] = mutantSite[-1]
                        mutatedSeq = ''.join(mutatedSeq)
                        if not mutatedSeq :
                            print('-------------')
                    else :
                        mutatedSeq = ''

                if mutatedSeq :      
                    entry = {
                        'ECNumber': ECNumber,
                        'Organism': Organism,
                        'Smiles': Smiles,
                        'Substrate': Substrate,
                        'Sequence': mutatedSeq,
                        'Type': 'mutant',
                        'Value': Value,
                        'Unit': Unit,
                    }

                    entries.append(entry)

        except:
            continue
    print(i)

    print(len(entries))  

    with open('./KM_combination_0816.json', 'w') as outfile :
        json.dump(entries, outfile, indent=4)

if __name__ == "__main__" :
    uniprotID_entry()
    uniprotID_noseq()
    nouniprotID_entry_uniprot()
    nouniprotID_entry_brenda()
    combine_sequence()
