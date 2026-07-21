#!/usr/bin/env python
import re
import os
import xml.etree.ElementTree as ET
import itertools,json
from itertools import chain

def Target_mRNA_List(mirt_fn,pat_fn,p53_network,topcount=15):

    boolean_network_species = p53_network.keys()
    with open(mirt_fn,'r') as fp: mrna_table = fp.readlines()

    # CINECA: if not abs path, join relative path to PWD folder    
    # if not os.path.isabs(pat_fn):
    #    pat_fn = os.getenv('PWD') + '/' + pat_fn
    with open(pat_fn,'r') as fp:
        xmlstring = re.sub(' xmlns="[^"]+"','',fp.read(),count=1)
    root = ET.fromstring(xmlstring)

    #---select columns from CSV data
    cols = [line.split(',')[1:4] for line in mrna_table[1:]]

    #---temporary utility function
    lookdown = lambda obj : [a.tag for a in list(obj)]

    #---argsort
    def argsort(seq): return sorted(range(len(seq)),key=seq.__getitem__)

    #---get internal data
    internal_data = root.findall('Sample')[0].findall('Data-Table')[0].findall('Internal-Data')[0].text
    internal_lookup = [(j,float(k)) for j,k in [i.split('\t') for i in internal_data.split('\n')
            if not re.match('^\s*$',i)]]

    #---get the top hits
    tops = [internal_lookup[j] for j in argsort(zip(*internal_lookup)[1])[-topcount:][::-1]]

    #---collect associated genes
    gene_all = []
    gene_miRNA = {}
    for subject in zip(*tops)[0]:
        rows_match = [ii for ii,i in enumerate(cols) if re.match('^%s-?'%subject,i[0])]
        #rows_hs = [i for i in rows_match if cols[i][1]=='Homo sapiens']
        genes = list(set([cols[i][2] for i in rows_match]))
        gene_all.append(genes)
        matched_species_permirna = list(set(genes) & set(boolean_network_species))
        gene_miRNA.update({subject:matched_species_permirna})
    failures = [tops[ii] for ii,i in enumerate(gene_all) if i==[]]

    gene_all_unpacked = list(chain.from_iterable(gene_all))

    boolean_network_species = p53_network.keys()
    matched_species = list(set(gene_all_unpacked) & set(boolean_network_species))

    return matched_species, gene_miRNA, gene_all_unpacked
