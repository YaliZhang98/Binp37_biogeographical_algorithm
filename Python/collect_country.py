# -*- coding: utf-8 -*-
"""
Title: collect_country.py
 Date: 2021-04-22
 Author: Yali Zhang 
 
 Description: 
     This program is used to collect individuals in different countries on different cotinents.
     That is, create a information form with cotinent, country, year, ID and number of SNPs
     of all individuals.
 
 Usage:
     python collect_country.py cotinent 
"""
from collections import defaultdict
import sys

# script can be used directly on terminal with name of continent namelist
cotinent = sys.argv[1]

name = open('{}'.format(cotinent),'r',encoding='UTF-8')
anno = open('10K_SNP_ind_annotation.anno','r',encoding='UTF-8')
location = open('country_location.txt','r',encoding='UTF-8')
Individuals = open('{}_header.tsv'.format(cotinent),'w',encoding='UTF-8')

# create directory for serial number and ID of individuals 
name_l = {}
for line in name:
    line = line.strip()
    line = line.split()
    name_l[line[1]]=line[0]

# create a dictionary to store annotation information of each individuals
dic_anno =  defaultdict(list)
anno.readline()
for line in anno:
    line = line.strip()
    line = line.split('\t')
    if line[1] in name_l :
        line[9] = line[9].strip()
        dic_anno[line[1]].append(line[9])
        dic_anno[line[1]].append(line[6])
        dic_anno[line[1]].append(line[18])

# create a dictionary to store continent information of each country
dic_location = {}
for line in location:
    line = line.strip()
    line = line.split('\t')
    dic_location[line[0]] = line[2]

# print annotation information of each individual we need
for ind in name_l.keys():
    country_location = dic_anno[ind][0].replace(' ','_')
    print("{}\t{}\t{}\t{}\t{}\t{}".format(name_l[ind],dic_location[country_location],dic_anno[ind][0],dic_anno[ind][1],dic_anno[ind][2],ind),file=Individuals)
    
name.close()
anno.close()
location.close()
Individuals.close()