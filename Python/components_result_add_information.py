# -*- coding: utf-8 -*-
"""
 Title: components_result_add_information.py

 Author: Yali Zhang 
 
 Description: 
     This program is used to add annotation information of individuals on namelist.

"""
from collections import defaultdict

name = open('36_components_values','r',encoding='UTF-8')
anno = open('10K_SNP_ind_annotation.anno','r',encoding='UTF-8')
location = open('country_location.txt','r',encoding='UTF-8')
fileout = open('36_components_results.tsv','w',encoding='UTF-8')
header = open('36_header.txt','r',encoding='UTF-8')
name_list = open('10K_SNP_namelist.txt','r',encoding='UTF-8')

filter_name = [] # get namelist for individuals with >10K SNPs
for line in name_list:
    line = line.strip()
    filter_name.append(line)

name_l = {}
for line in name:
    line = line.strip()
    line2 = line.split()
    name_l[line2[1]]=line

dic_anno =  defaultdict(list)
anno.readline()
for line in anno:
    line = line.strip()
    line = line.split('\t')
    if line[1] in name_l.keys() :
        line[9] = line[9].strip()
        dic_anno[line[1]].append(line[9]) # country
        dic_anno[line[1]].append(line[6]) # full year
        dic_anno[line[1]].append(line[5]) # Date mean in BP
        dic_anno[line[1]].append(line[10]) # Lat
        dic_anno[line[1]].append(line[11]) # Long
    

dic_location = defaultdict(list)
for line in location:
    line = line.strip()
    line = line.split('\t')
    dic_location[line[0]].append(line[1]) # continent
    line[2] = line[2].replace(' ','_')
    dic_location[line[0]].append(line[2]) # detailed continent

for line in header:
    line = line.strip()
    print(line,file=fileout)

# print header to the output file
for ind in filter_name:
    country_location = dic_anno[ind][0].replace(' ','_')
    print("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format(name_l[ind],dic_location[country_location][0],dic_location[country_location][1],dic_anno[ind][0],dic_anno[ind][1],ind,dic_anno[ind][2],dic_anno[ind][3],dic_anno[ind][4]),file=fileout)


name.close()
name_list.close()
anno.close()
location.close()
fileout.close()
header.close()