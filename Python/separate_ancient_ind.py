# -*- coding: utf-8 -*-
"""
Title: separate_ancient_ind.py
 Date: 2021-04-20
 Author: Yali Zhang 
 
 Description: 
     This program is used to separate individuals according to their year in
     annotation file and then to create a namelist of ancient individuals.
"""

name = open('10K_namelist.tsv','r',encoding='UTF-8')
anno = open('10K_SNP_ind_annotation.anno','r',encoding='UTF-8')
ancient_namelist = open('ancient_namelist.tsv','w')

# find individual IDs if they are not the present people according to the annotation file
ancient = []
anno.readline()
for line in anno:
    line = line.strip()
    line = line.split('\t')   
    if line[6] != '..' and line[6] != 'present':        
        ancient.append(line[1])

# extract name list of ancient individuals
for line in name:
    line = line.strip()
    line_split = line.split()
    if line_split[1] in ancient:
        print(line,file=ancient_namelist)
        
name.close()
anno.close()
ancient_namelist.close()
        