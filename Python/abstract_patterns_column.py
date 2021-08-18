# -*- coding: utf-8 -*-
"""

Author: Yali Zhang 
Date : May 20
Description: This program is used to collect pattern column from P files.

"""

from collections import defaultdict
import os

name = open('patterns_column.txt','r')

# change the components names to the location of columns
def num(alp):
    if alp == 'a':
        number = 0
    elif alp == 'b':
        number = 1
    elif alp == 'c':
        number = 2
    elif alp == 'd':
        number = 3
    elif alp == 'e':
        number = 4
    elif alp == 'f':
        number = 5
    elif alp == 'g':
        number = 6
    elif alp == 'h':
        number = 7
    elif alp == 'i':
        number = 8
    return(number)

# function that can create a new directory
def mkdir(path):
    path = path.strip()
    isExisits = os.path.exists(path)
    if not isExisits:
        os.makedirs(path)
        return True

# create a directory to store the location of column
dic_patterns = defaultdict(list)
for line in name:
    line = line.strip()
    line = line.split(':')
    continent = line[0]
    patterns = line[1].split(',')
    for pattern in patterns:
        dic_patterns[continent].append(pattern)
  
# print out selected components in Q file
for key in dic_patterns:
     print(key)
     for value in dic_patterns[key]:
        print(value)
        column = num(value[0])
        K = value[1:]
        if '_' in key: # for different constuction of file name with "_“
            filein = open('Analysis/{}/ADMIXTURE/{}_QC.{}.P'.format(key,key,K),'r')
        else:    
            filein = open('Analysis/{}/ADMIXTURE/{}_ind_QC.{}.P'.format(key,key,K),'r')
        
        path = 'Single_pattern_column/{}'.format(key)
        mkdir(path)
        fileout = open('{}/{}_{}.P'.format(path,key,value),'w')
        
        for line in filein:
            line = line.strip()
            line = line.split()
            col = line[column]
            print(col,file=fileout)

        filein.close()
        fileout.close()
            
name.close()


