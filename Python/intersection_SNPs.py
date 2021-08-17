# -*- coding: utf-8 -*-
"""

name_file = open('SNPid_Oceanians','r')
name_list = []

for line in name_file:
    line = line.strip()
    name_list.append(line)

continent = ['Americans','Africans','Asians','Asians_subset1','Europeans_east','Europeans_north'
             ,'Europeans_west','Europeans_north','Oceanians']
dic_SNP = {}
for con in continent:
    filein = open('{}_components'.format(con),'r')
    fileout = open('{}_filter'.format(con),'w')
    header = filein.readline()
    header = header.strip()
    print(header, file=fileout)
    for line in filein:
        line = line.strip()
        SNP = line.split()[0]
        dic_SNP[SNP] = line
    for obj in name_list:
        print(dic_SNP[obj],file=fileout)
    filein.close()
    fileout.close()
   
name_file.close()