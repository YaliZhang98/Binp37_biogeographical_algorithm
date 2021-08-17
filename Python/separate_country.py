# -*- coding: utf-8 -*-
"""
 Title: separate_country.py
 Date: 2021-04-20
 Author: Yali Zhang 
 
 Description: 
     This program is used to separate individuals according to their country in
     annotation file and then to create a namelist of each country individuals.
"""
from collections import defaultdict

name = open('ancient_namelist.tsv','r',encoding='UTF-8')
anno = open('10K_SNP_ind_annotation.anno','r',encoding='UTF-8')

# creat a dictionary to store ID of individuals
dic_name = {}
for line in name:
    line = line.strip()
    line = line.split()
    dic_name[line[1]] = line[0]

# create a dictionary to store country and ID of individuals 
dic_country =  defaultdict(list)
anno.readline()
for line in anno:
    line = line.strip()
    line = line.split('\t')
    if line[9] == '..': # some individuals are missing country notes
        dic_country['Missing'].append(line[1])
    else:        
        dic_country[line[9]].append(line[1])

# print IDs of individuals from the same country
for key in dic_country.keys():
    name_c = key
    name_c = name_c.strip()
    if ' ' in name_c:
        name_c = name_c.replace(' ','_') # change space in name of individual
    country = open('./Country/{}'.format(name_c),'w')
    for value in dic_country[key]:
        if value in dic_name.keys():
            print('{}\t{}'.format(dic_name[value],value),file=country)
    country.close()

name.close()
anno.close()
