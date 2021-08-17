# -*- coding: utf-8 -*-
"""
 Title: annotation.py
 Date: 2021-04-15
 Author: Yali Zhang 
 
 Description: 
     This program is used to find annotation information of individuals whose 
     non-missing SNP number is higer than 10K
"""

name_SNP = open('10K_name_SNP.tsv','r',encoding='UTF-8')
annotation_origial = open('v44.3_HO_public.anno','r', encoding='UTF-8')
annotation_10K = open('10K_SNP_ind_annotation.anno','w',encoding='UTF-8')

# create dictionary of name list with number of SNPs
dic_SNP = {}
name_SNP.readline()
for line in name_SNP:
    line = line.strip()
    name = line.split()[0]
    SNP = line.split()[1]
    dic_SNP[name] = SNP

# extract header of annotation file and add a new column to store the number of SNP for each individual
header = annotation_origial.readline()
header = header.strip()
print('{}\tnon-missing SNP number'.format(header),file=annotation_10K)

# extract individuals with annotation according to the name list
for line in annotation_origial:
    line = line.strip()
    anno_id = line.split()[1]
    if anno_id in dic_SNP.keys():
        SNP_num = dic_SNP[anno_id]
        print('{}\t{}'.format(line, SNP_num),file=annotation_10K)

name_SNP.close()
annotation_origial.close()
annotation_10K.close()
    

    