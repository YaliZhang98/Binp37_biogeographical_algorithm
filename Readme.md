# Project - Map the spread of Yersinia pestis over time using ancient DNA with machine learning  
### Course: Binp37  
### Credits: 30  
### Student: Yali Zhang  
### Supervisor: Eran Elhaik   

This readme file contain work flow and all the scripts used in project. Scripts in python and r are all recorded in this file and can be found in directory Python and R.  
This project constructed an algorithm that can locate the origin of ancient individuals and map the spread of Y. Pestis.   

## 1 Construction of Gene Pool   
In this part, it needs three tools: eigensoft ( version 7.2.1), P-link (version 1.07), ADMIXTURE (version 1.3.0).  
eigensoft is used to transfer format of files which contains genotypes information of ancient individuals from EIGENSTRAT to PED format. P-link is used to do the quality control of database, filtration of ancient individuals, and testing of missingness information for individuals SNPs. ADMIXTURE is used to calculate allele frequencies patterns of ancient individuals.  
```shell
mkdir Binp37 # create a directory binp37 as a subdirectory in home directory
cd Binp37
```
### 1.1 Tools installation

#### 1.1.1 Install conda 
If you have been installed conda, please skip this step,  
```shell
mkdir ~/Binp37/Program
cd ~/Binp37/Program
wget https://repo.continuum.io.miniconda/Miniconda3-latest-Linux-x86_64.sh -o ~miniconda.sh # download conda
ls # whether there is a file named miniconda*.sh
bash Miniconda3-latest-Linux-x86_64.sh # Run the bash program
cd miniconda3/
cd etc
cd profile.d/
pwd # find the route of conda installed
export PATH=route:$PATH # add the route into your ~/.bashrc folder so that we can use it anywhere. eg.export PATH=/home/miniconda3/etc/profile.d:$PATH
source ~/.bashrc
```

#### 1.1.2  Install  eigensoft-7.2.1  
```shell
conda config --add channels bioconda
conda create --name binp37 # create a new environment to store programs in this project
conda activate binp37
conda install eigensoft # version 7.2.1
```

#### 1.1.3 Install PLINK
```shell
wget http://zzz.bwh.harvard.edu/plink/dist/plink-1.07-x86_64.zip # v1.07
unzip plink-1.07-x86_64.zip
cd /usr/local/bin # To your bin directory
sudo cp Binp37/Program/plink-1.07-x86_64/plink .
plink # can be used in any directory
```

#### 1.1.4 ADMIXTURE installation
```bash
conda install admixture
```

### 1.2 Database preparation

#### 1.2.1 Transfer file format through convertf 
Download Ancient DNA database (1240K+HO) from Allen Ancient DNA Resource ("Allen Ancient DNA Resource https://reich.hms.harvard.edu/allen-ancient-dna-resource-aadr-downloadable-genotypes-present-day-and-ancient-dna-data", version 44.3). Unzip files and store them in new directory binp37/Database/PED.  
Transfer EIGENSTRAT to PED format.  
```shell
mkdir ~/Binp37/Database/PED
cd ~/Binp37/Database/PED # database files downloaded from websit are stored here. There are three files: v44.3_HO_public.geno, v44.3_HO_public.snp, v44.3_HO_public.ind
convertf -p parfile
# parfile:
# 	genotypename:    v44.3_HO_public.geno
# 	snpname:         v44.3_HO_public.snp
# 	indivname:       v44.3_HO_public.ind
# 	outputformat:    PED
# 	genotypeoutname: v44.3_HO_public.ped
# 	snpoutname:      v44.3_HO_public.pedsnp
# 	indivoutname:    v44.3_HO_public.pedind

mv v44.3_HO_public.pedsnp v44.3_HO_public.bim # change file name for better understanding and usage
mv v44.3_HO_public.pedind v44.3_HO_public.fam
cp v44.3_HO_public.bim v44.3_HO_public.map 
```

#### 1.2.2 Filter individual with >10K SNPs
```shell
# find missingness
plink --file v44.3_HO_public --missing --out missingness

# change format of missingness file in order to use it in R
cat missingness.imiss | sed "s/^ * //g" | sed "s/ * /\t/g" > missingness.imiss.tsv # produce a tsv file
```
Name list of Individuals with >10K SNPs is extracted through R:  
```R
library(tidyverse)

missingness <- read_delim("missingness.imiss.tsv", delim = "\t")
head(missingness)
summary(missingness)

# transfer character to numeric
missingness$N_MISS <- as.numeric(missingness$N_MISS)
summary(missingness)

# calculate number of SNP of each individual
missingness$N_SNP <- missingness$N_GENO - missingness$N_MISS
head(missingness)

# order the table according to the number of SNP of each individual
sort_missingness <- missingness[order(missingness$N_SNP), ]

# extract individuals with >10K SNPs 
tenK_ind <- subset(missingness, missingness$N_SNP > 10000)
sort_tenK_ind <- tenK_ind[order(tenK_ind$N_SNP), ]

# output the namelist of individual with >10K SNP
name <- tenK_ind[,c(1,2)]
write.table(name, "10K_namelist.tsv", sep="\t", row.names = F)

# output the namelist of individual with >10K SNP and the number of non-missing SNPs
name_SNP <- tenK_ind[,c(2,7)]
write.table(name_SNP, "10K_name_SNP.tsv", sep="\t", row.names = F)
```
After getting name list of individuals with more than 10K SNPs, the database is filtered according to this list by plink:  
```shell
# extract individual
cat 10K_namelist.tsv | sed "s/\"//g" > 10K_namelist2.tsv
rm 10K_namelist.tsv
mv 10K_namelist2.tsv 10K_namelist.tsv

# create binary file of database
plink --file v44.3_HO_public --make-bed --out binary

# filtration of database
plink --bfile binary --keep 10K_namelist.tsv --recode --make-bed --out 10K_keep
```

#### 1.2.3 Preparation of annotation file
Filtration annotation file according to the 10K_namelist. Then add number of non-missing SNPs information in annotation file.  
```bash
cat 10K_name_SNP.tsv | sed "s/\"//g" > 10K_name_SNP2.tsv
rm 10K_name_SNP.tsv
mv 10K_name_SNP2.tsv 10K_name_SNP.tsv
```
Prodece Annotation file through Python (annotation.py):  
```python
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
```

#### 1.2.4 Separate ancient individuals
Through python script ( "separate_ancient_ind.py") to separate individuals according to their year and create a namelist of ancient people.  
```python
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
```
After getting name list of ancient individuals, filter database through P-link:  

```bash
plink --bfile 10K_keep --keep ancient_namelist.tsv --recode --make-bed --out ancient_ind
cat ancient_ind.fam | wc -l # 4980 ancient people
```

### 1.3 Unsupervised ADMIXTURE analysis
Admixture can be used to estimate the maximum likelihood of individual ancestors from a multi-site SNP genotype dataset, so we use this method to analyze the ancient DNA dataset to obtain individual genetic characteristics.  

#### 1.3.1 Create namelist of each country 
Through python script (separate_country.py) to separate individuals according to thier country.  
```python 
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
```
Create directories of each group that will be analyzed in ADMIXTURE.  

```bash
# create folder for each continents
mkdir ~/Binp37/ADMIXTURE/Analysis
cd ~/Binp37/ADMIXTURE/Analysis
mkdir Americans
mkdir Asians
mkdir Asians_subset1
mkdir Africans
mkdir European_east
mkdir European_north
mkdir European_south
mkdir European_west
mkdir Oceanians
```

#### 1.3.2 Admixture analysis on all Americans ( in one big dataset with all Americans)
ADMIXTURE analysis on Americans database:  

```bash
mkdir ~/Binp37/Database/Country_individual/Americans
cd ~/Binp37/Database/Country_individual/Americans
# copy all American country namelist into this folder and create a namelist having all Americans
ls | while read line; do cat ${line} >> Americans; done # 526 individuals
cp Americans ../..

# Go to folder which contain ancient_ind batabase
cd  ~/Binp37/Database/
plink --bfile ancient_ind --keep Americans --recode --make-bed --out Americans_ind 
plink --bfile Americans_ind --geno 0.999 --make-bed --out Americans_ind_QC # quality control: 0 variants removed
mv Americans* ~/Binp37/Database/Country_individual/Americans

# ADMIXTURE analysis
cd ~/Binp37/Database/Country_individual/Americans
mkdir ADMIXTURE
cd ADMIXTURE
for K in 2 3 4 5 6 7 8 9 ; do echo "nohup admixture --cv ../Americans_ind_QC.bed $K -j64 | tee log${K}.out &" | bash ; done

# Copy ADMIXTURE analysis result to continent directories in ~Binp37/ADMIXTURE
cd ..
tar -cvzf ADMIXTURE.tar.gz ADMIXTURE
cp ADMIXTURE.tar.gz ~Binp37/ADMIXTURE/Americans
```
After getting ADMIXTURE result, it need merge ID, annotation information and all ADMIXTURE value with different K value into one table for following analysis.  
Create a txt file which contains the continent and subregion information of each country named country_location.txt in following format:  
```txt
Cameroon	Africa	Central Africa
China	Asia	East Asia
Cuba	America	Middle America
Sweden	Europe	North Europe
```
Copy annotation file and namelist of Americans into this directory. Through python (collect_country.py) to extract ID, full year, number of SNPs, country of individuals for Americans from annotation file into file "Americans_header.tsv":   
```python
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
```
Usage of this python script on bash terminal:  
```bash
python collect_country.py Americans
```
Through R (admixture_table.R) to collect the values of all .Q files with K = 2-9 into one table "Americans_ADMIXTURE.csv" (in ~/Binp37/ADMIXTURE/Americans)  
```bash
tar -xvzf ADMIXTURE.tar.gz
```
```R
library(readr)

# list name of Q files in directory
myfiles <- Sys.glob("*.Q")

# check admixture pattern on K = 4
tb4 = read.table(myfiles[3]) # read file
tb4_2 <- round(tb4,digits=2) # change format of number in Q file
barplot(t(as.matrix(tb4)),col = rainbow(3),
        xlab = "Individual #",
        ylab = "Ancestry",
        border = NA)

tb2 = read.table(myfiles[1])
tb2_2 <- round(tb2,digits=2)

tb3 = read.table(myfiles[2])
tb3_2 <- round(tb3,digits=2)

tb5 = read.table(myfiles[4])
tb5_2 <- round(tb5,digits=2)

tb6 = read.table(myfiles[5])
tb6_2 <- round(tb6,digits=2)

tb7 = read.table(myfiles[6])
tb7_2 <- round(tb7,digits=2)

tb8 = read.table(myfiles[7])
tb8_2 <- round(tb8,digits=2)

tb9 = read.table(myfiles[8])
tb9_2 <- round(tb9,digits=2)

name <- read_tsv("Americans_header.tsv",col_names = FALSE) # read annotation information
order_name <- name[order(name$X1), ]
header <- order_name[,c(2,3,4,5,6)]

Country_ind <- cbind(header,tb2_2,tb3_2,tb4_2,tb5_2,tb6_2,tb7_2,tb8_2,tb9_2)
colnames(Country_ind) <- c('Continents','Country','Full year','SNPs','ID','a2','b2','a3','b3','c3','a4','b4','c4','d4','a5','b5','c5','d5','e5','a6','b6','c6','d6','e6','f6','a7','b7','c7','d7','e7','f7','g7','a8','b8','c8','d8','e8','f8','g8','h8','a9','b9','c9','d9','e9','f9','g9','h9','i9')

write.csv(Country_ind,"Americans_ADMIXTURE.csv")
```
#### 1.3.3 Admixture analysis on other continent and regions
For the ADMIXTURE analysis of other continents and regions, the same operations as the Americans will be carried out. Repeat the code in 1.3.2, just need to change all "Americans" to the continent currently analyzed.  

#### 1.3.4 Select components from ADMIXTURE analysis result
Select out columns that can reflect the allele frequencies characteristics for individuals in a particular region from all ADMIXTURE analysis result of each continent and region.  

#### 1.3.5 Extract components' columns in P-files 
The ancestry fractions in the P file corresponding to the column of the components selected should be extracted and merged for use in supervising ADMIXTURE analysis.  
```bash
cd ~/Binp37/ADMIXTURE
```
Create a txt file that contain the selected components from ADMIXTURE result table named "patterns_column.txt" as follows:   
```txt
Americans:c9,h8,h9,e9,c8,f7,f8,d8,f9
Oceanians:a2,b2
Africans:g9,b6,e7,d7,e9,d6,f8,e6
Asians:i9,e9,g9,a9,h9
Asians_subset1:i10,d10,d6,h10,e10,b6,f6
Europeans_north:a2,b2
Europeans_south:a4,e7,g7,b4,d4
Europeans_east:a4,b4,c4,d4
Europeans_west:a8,f8,d8,b8,a10,c8,c5
```
Extract columns in P-files according to identified components through Python (abstract_patterns_column.py):   
```python
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
```
Then add the ids of SNPs in P file columns and use file name as header.  
```bash
# Create a header file of all the columns according to the components 

ls | while read continent; do ls ${continent}/*.P | while read line; do name=`echo ${line} | sed "s/${continent}\///g" | sed s/.P//g ` ; echo ${name} >> ${continent}/${continent}_head ; done;done

ls | while read continent; do cat ${continent}/*_head | tr '\n' '\t' > ${continent}/${continent}_header; done

ls | while read continent; do sed -i 's/^/SNP_id\t/' ${continent}/${continent}_header; done

# add ID of SNPs of each row in components file

 ## the IDs are extracted from the second column in .bim file 
ls | while read line; do cat ${line}/*_QC.bim | cut -f 2 > SNPid_${line}; done # the site is where all the bim files in thier continents' folder
 ## move these id files to each continents component folder
 
# merge the heads and columns with ids
ls | while read line ; do dos2unix ${line}/* ; done
ls | while read line ; do paste ${line}/*.P >> ${line}/${line}_all; done
ls | while read line ; do paste ${line}/SNPid_${line} ${line}/${line}_all > ${line}/${line}_id; done
ls | while read line ; do awk '{print $0}' ${line}/${line}_header ${line}/${line}_id > ${line}/${line}_components ; done
```

#### 1.3.6 Find the intersection SNPs in all the continents

1. The IDs of SNPs are extracted from the second column in .bim file. This step will get 9 files for different continents and regions. (SNPid_Africans, SNPid_Asians, SNPid_Europeans_east, SNPid_Europeans_south, SNPid_Oceanians, SNPid_Americans, SNPid_Asians_subset1, SNPid_Europeans_north, SNPid_Europeans_west)  
```bash
ls | while read line; do cat ${line}/*_QC.bim | cut -f 2 > SNPid_${line}; done
cat SNPid_Oceanians | wc -l # 593463 (Oceania has the least number of SNPs)
```
2. Add all ids of SNPs in all header files of all components into one file.  
```bash
ls | while read line; do cat $line >> all.txt; done
```
3. Find SNPs that exist in all components file. That means in 'all.txt' file, these intersection SNPs can be found 9 times.  
```bash
cat all | sort | uniq -c | grep ' 9 ' | wc -l # 593463
```
4. Because the Oceania has the least number of SNPs and it is also 593463, that means the SNPs in Oceania are the  SNPs that overlapped in all continents and regions. It can use id list of SNPs of Oceania to extract these SNPs in other continents.  

5. Extract the rows of intersect SNPs of all continents in the order of IDs in the Oceania SNP_id file through python script intersection_SNPs.py .  
```python
name_file = open('SNPid_Oceanians','r')
name_list = []

# add all the SNP ids in Oceanians into a list
for line in name_file:
    line = line.strip()
    name_list.append(line)

# list all the continents and regions that contains the components in P-file
continent = ['Americans','Africans','Asians','Asians_subset1','Europeans_east','Europeans_north','Europeans_west','Europeans_south','Oceanians']

# fiter the valus in P-file columns according to the intersesction SNP ids in Oceanians SNPs_id file 
dic_SNP = {} 
for con in continent:
    filein = open('Before_filter/{}_components'.format(con),'r')
    fileout = open('After_filter/{}_filter_2'.format(con),'w')
    header = filein.readline()
    header = header.strip()
    print(header, file=fileout)
    # create a dictionary: SNP ids as key and this line as value
    for line in filein:
        line = line.strip()
        SNP = line.split()[0]
        dic_SNP[SNP] = line
    
    # print SNPs information according to the order in Oceanians SNP id list.
    # Here just when Oceanian SNPs exist in full in all the other datasets, the value of components 
    # can be printed with SNP ids in Oceanian as key. 
    for obj in name_list:
        if obj in dic_SNP.keys(): # double check the Oceanians SNPs also in other dataset
            print(dic_SNP[obj],file=fileout)
    filein.close()
    fileout.close()
   
name_file.close()
```

6. Merge the component files of all continents into one.  

```bash
ls | while read line; do sed -r 's/(\S+\s+)//' $line > ../After_filter_remove_header/${line}_noheader; done

sed -r 's/(\S+\s+)//' Africans_filter > Africans_filter_noheader

# check the number of line in each file
ls | while read line; do cat $line | wc -l ; done # 593646

# merge all the files
dos2unix *
paste * > components_final

# check the number of column and line in the final components file
cat components_final | wc -l # 593464 (593643 SNPs + 1 header)
cat components_final | awk 'END{print NF}' # 50 (49 compponents + 1 SNP_id)

# delete header of p files after filter
ls | while read line; do sed -i '1d' $line ; done
```

#### 1.3.7 filter the dataset according to the SNP
```bash
plink --bfile Africans_ind_QC --extract SNPid_Oceanians --make-bed --out Africans_filter
plink --bfile Americans_ind_QC --extract SNPid_Oceanians --make-bed --out Americans_filter
plink --bfile Asians_ind_QC --extract SNPid_Oceanians --make-bed --out Asians_filter
plink --bfile Asians_subset1_QC --extract SNPid_Oceanians --make-bed --out Asians_subset1_filter
plink --bfile Europeans_east_QC --extract SNPid_Oceanians --make-bed --out Europeans_east_filter
plink --bfile Europeans_west_QC --extract SNPid_Oceanians --make-bed --out Europeans_west_filter
plink --bfile Europeans_north_QC --extract SNPid_Oceanians --make-bed --out Europeans_north_filter
plink --bfile Europeans_south_QC --extract SNPid_Oceanians --make-bed --out Europeans_south_filter
plink --bfile Oceanians_ind_QC --extract SNPid_Oceanians --make-bed --out Oceanians_filter

# filter ancient dataset
plink --bfile ancient_ind --extract SNPid_Oceanians --make-bed --out ancient_filter
```

### 1.4 Supervised ADMIXTURE analysis

#### 1.4.1 ADMIXTURE analysis on filtration dataset
```bash
admixture32 ancient_filter.bed -f #gene pools
```
#### 1.4.2 Result analysis of gene pool file
When get output result files of supervised ADMIXTURE analysis, through following code can create a gene pool table which is convenient for analysis. After selecting components in gene pool make sense, remove remaining components in database and then run the supervised ADMIXTURE analysis again. Filter the components in gene pool until all the components make sense.  
```bash
mkdir ~Binp37/ADMIXTURE/Result
cd ~Binp37/ADMIXTURE/Result
mv ../value_merge_32.fam .
mv ../merge_35_small.32.Q .

# create component header from fam (merge_35_small.fam) file 
 ## copy ID column in fam file to a new txt file named 32_header_1.txt
cat 32_header_1.txt | cut -d ' ' -f 1 | uniq > 32_header_2.txt
cat 32_header_2.txt | tr "\n" "\t" > 32_header.txt

# create id file of individuals
awk '{print $1,$2}' merge_32_small.fam > 32_ind_id
sed -i 's/ /\t/g' 32_ind_id

# replace space in value file into tab
sed -i 's/ /\t/g' merge_32_small.32.Q

# merge id file and value file 
paste 35_ind_id merge_32_small.32.Q > 32_components_values
```

#### 1.4.3 Rename gene pools name
```bash
# change same components header in different files
sed -i "s/East_Africa_subset1_d7/E_Af_2.5KBCE_1.6KCE/g" `grep East_Africa_subset1_d7 -rl ./`

# rename of 36 gene pools
cat 36_gene_pools_newname.txt | while read line; do oldname=`echo ${line} | cut -d ' ' -f 1 | sed 's/\n//g' ` ; newname=`echo ${line} | cut -d ' ' -f 2 | sed 's/\n//g'`; sed -i "s/${oldname}/${newname}/g" allrefs_36.fam ; done 
sed -i "s/Oceania_subset1_a2/Oceania_1KBCE_2KCE/g" allrefs_36.fam

# rename of 32 gene pools
cat 32_gene_pools_newname.txt | while read line; do oldname=`echo ${line} | cut -d ' ' -f 1 | sed 's/\n//g' ` ; newname=`echo ${line} | cut -d ' ' -f 2 | sed 's/\n//g'`; sed -i "s/${oldname}/${newname}/g" allrefs_32.fam ; done 
sed -i "s/Oceania_subset1_a2/Oceania_1KBCE_2KCE/g" allrefs_32.fam
```

#### 1.4.4 Result analysis of final gene pool

1. Extact id of individuals that have SNP more than 10K  
```bash
cat plink_36.imiss | tr -s ' ' '\t' > plink_36_tab.imiss
```
Through Excel to calculate the number of non-missing SNPs of individuals. After sort the number of non-missing SNPs,  then extract the name list of these individuals with SNPs > 10K into a new txt file named 10K_SNP_namelist.txt  

2. Create files that can be used in Python for annotation excel result.  
```bash
# create component header from fam file (36_header_1.txt)
cat 36_header_1.txt | cut -d ' ' -f 1 | uniq > 36_header_2.txt
cat 36_header_2.txt | tr "\n" "\t" > 36_header.txt
# add '#	ID	' into head file 

# create pure individual id from fam file in a new file (merge_36_small.fam)
# create id file of individuals
awk '{print $1,$2}' merge_36_small.fam > 36_ind_id
sed -i 's/ /\t/g' 36_ind_id

# replace space in value file into tab
sed -i 's/ /\t/g' merge_36_small.36.Q

# merge id file and value file 
paste 36_ind_id merge_36_small.36.Q > 36_components_values
```
3. Add annotations of each individual into the gene pool file  
```python
from collections import defaultdict

name = open('36_components_values','r',encoding='UTF-8')
anno = open('10K_SNP_ind_annotation.anno','r',encoding='UTF-8')
location = open('country_location.txt','r',encoding='UTF-8')
fileout = open('36_gene_pool.tsv','w',encoding='UTF-8')
header = open('36_header.txt','r',encoding='UTF-8')
name_list = open('10K_SNP_namelist.txt','r',encoding='UTF-8')

filter_name = [] # get namelist for individuals with >10K SNPs
for line in name_list:
    line = line.strip()
    filter_name.append(line)

name_l = {} # create a dictionary to store ID and ADMIXTURE values
for line in name:
    line = line.strip()
    line2 = line.split()
    name_l[line2[1]]=line

# create a dictionary to store individuals ID and annotation information 
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
    
# create a dictionary to store country and corresponding continent or region 
dic_location = defaultdict(list)
for line in location:
    line = line.strip()
    line = line.split('\t')
    dic_location[line[0]].append(line[1]) # continent
    line[2] = line[2].replace(' ','_')
    dic_location[line[0]].append(line[2]) # detailed continent

# print header to the output file
for line in header:
    line = line.strip()
    print(line,file=fileout)

for ind in filter_name:
    country_location = dic_anno[ind][0].replace(' ','_')
    print("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format(name_l[ind],dic_location[country_location][0],dic_location[country_location][1],dic_anno[ind][0],dic_anno[ind][1],ind,dic_anno[ind][2],dic_anno[ind][3],dic_anno[ind][4]),file=fileout)

name.close()
name_list.close()
anno.close()
location.close()
fileout.close()
header.close()
```

#### 1.4.5 ADMIXTURE proportions plot of 36 gene pool
Scripts are recorded in file admixture_proportion.r:   
```R
Dataset <- read.table(file="36_gene_pool.tsv",header=TRUE)
Dataset <- Dataset[,-5:-9]
Dataset <- Dataset[,-1:-3]

png("admixture_proportion.png", width = 13,height = 8, units = 'in', res = 600)

par(xpd = T, mar = par()$mar + c(1,0,0,7), mgp = c(0,0.7,0), las=2)
bp <- barplot(t(as.matrix(Dataset)),col = c( "brown","red3","maroon1","slateblue1","olivedrab1","darkorchid4","orange","cyan4","cyan2","mediumorchid1","seashell2","dodgerblue3","violetred3","palevioletred1","rosybrown1","mediumspringgreen","gold2","orangered","black"),
              ylab = "",
              border = NA)

mtext(text = c(Dataset$Continent_detail), side = 1, at = bp, line = 0, padj = 1, cex = 0.5,las = 2)

par(mar=c(5, 4, 4, 2) + 0.1)
dev.off()
```


## 2 Origin prediction algorithm
Create a algorithm that can predict original coordinates of ancient individuals through supervised machine learning in R  

### 2.1 Installation of special packages in R
```r
devtools :: install_github("hrbrmstr/nominatim")
devtools :: install_github("adeckmyn/mapdataNE") # package for map plot
```
### 2.2 Original coordinates prediction algorithm in R
The result gene pool database should be preprocessed.   
Transfer 34_gene_pool.tsv into csv format in Excel  
Remove following columns from 34_gene_pool.csv and save file into a new file named 34_gene_pool_new.csv: #, ID, Full year. Then change the order of columns in the file as: Continent, Continent_detail, Country, Date mean in BP, Lat., Long., As_Ru_3K_0.9K_BCE, ...... (all components in gene pool) 34_gene_pool_new.csv  
Biogeographical algorithm is recorded in R file: Biogeographical_algorithm.R

```r
library(nominatim)
library(mapdataNE)
library(plyr)
library(dplyr)
library(sp)
library(rworldmap) 
library(caret)  
library(rpart)
library(maps)
library(MASS)
library(geosphere)
library(doParallel)
library(forcats)
library(DMwR)
library(png)

# Dataset preparation ----------------------------------------------------------

Dataset <- read.csv(file="34_gene_pool_new.csv",header=TRUE)

dim(Dataset)
sapply(Dataset,class)

Dataset$Country <- make.names(Dataset$Country)
Dataset$Continent_detail <- make.names(Dataset$Continent_detail)
Dataset$Continent <- make.names(Dataset$Continent)

Dataset$Continent_detail <- factor(Dataset$Continent_detail)
Dataset$Continent <- factor(Dataset$Continent)
Dataset$Country <- factor(Dataset$Country)

# longitude adjustment : adjust longitude near the boundary especially Oceania and east Asia
for (i in (1:length(Dataset$Long.))){
  if (Dataset[i,"Continent"] == "Asia" & Dataset[i,"Long."] < 0 ){
    Dataset[i,"adj_Long."] <- 360 + Dataset[i,"Long."]
   }else{
     Dataset[i,"adj_Long."] <- Dataset[i,"Long."]
   }
  if (Dataset[i,"Continent"] == "Oceania" & Dataset[i,"Long."] < 0){
      Dataset[i,"adj_Long."] <- 360 + Dataset[i,"Long."]
  }
}


# Extract headers of components
optimumVars <- names(Dataset)[7:40] # should contain all the components columns


# model construction -----------------------------------------------------------
ML_Location <-  function(training,testing,classTarget,variables){
  
  set.seed(1234) # create a same random numeric list
  
  # create 5 folds for 5 folds cross-validation
  folds <- createFolds(training[,classTarget], k = 5, returnTrain = T) # used in index item in parameters below (trainControl)
  
  # used in train model
  trControlClass <-  trainControl(
    method = "cv",
    number = 5,  
    verboseIter = FALSE,
    returnData = FALSE,
    search = "grid",
    savePredictions = "final",
    classProbs = T, # different with trControl
    allowParallel = T,
    index = folds )
  
  # used in train model
  trControl <-  trainControl(
    method = "cv",
    number = 5,  
    verboseIter = FALSE,
    returnData = FALSE,
    search = "grid",
    savePredictions = "final",
    allowParallel = T,
    index = folds)
  
  # used in train model
  tune_grid <- expand.grid(
    nrounds = c(400,600),
    eta = c( 0.05, 0.1),
    max_depth = c(3,6,9),
    gamma = 0,
    colsample_bytree = c(0.6,0.8),
    min_child_weight = c(1),
    subsample = (0.7)
  )
  
  
  ##### model training part ____________________
  
  # prediction model for continent
  Xgb_region <- train(x = training[,variables],y = training[,"Continent"], 
                      method = "xgbTree", # train model selection
                      trControl = trControlClass,
                      tuneGrid = tune_grid,
                      nthread = 1)
 
  print('continent training model has been done')
  
  l1_train <- data.frame(training[,c(variables)],Xgb_region[["pred"]][order(Xgb_region$pred$rowIndex),levels(training[,"Continent"]) ])
  
  # prediction model for continent subregions
  Xgb_class <- train(x = l1_train,y = training[,classTarget],
                     method = "xgbTree",
                     trControl = trControlClass,
                     tuneGrid = tune_grid,
                     nthread = 1)
  
  print('continent detail training model has been done')
  
  l2_train <- data.frame(l1_train,Xgb_class[["pred"]][order(Xgb_class$pred$rowIndex),levels(training[,classTarget]) ])

  ##### country prediction: country label can be added into chain prediction
  # Xgb_country <- train(x = l2_train,y = training[,'Country'],
  #                      method = "xgbTree",
  #                      trControl = trControlClass,
  #                      tuneGrid = tune_grid,
  #                      nthread = 1)
  # 
  # print('country training model has been done')
  # 
  # l22_train <- data.frame(l2_train,Xgb_country[["pred"]][order(Xgb_country$pred$rowIndex),levels(training[,'Country']) ])
  # 
  #
  ##### date prediction: age label (date mean in BP) can be added into chain prediction
  # Xgb_date <- train(x = l2_train,y = training[,"Date.mean.in.BP"],
  #                   method = "xgbTree",
  #                   trControl = trControl, #!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  #                   tuneGrid = tune_grid,
  #                   nthread = 1)
  # 
  # print('date training model has been done')
  # 
  # l_date_train <- data.frame(l2_train, "datePred" = Xgb_date[["pred"]][order(Xgb_date$pred$rowIndex),"pred" ])
 
  # prediction model for latitude
  Xgb_latitude <- train(x = l2_train ,y = training[,"Lat."], 
                        method = "xgbTree",
                        trControl = trControl,
                        tuneGrid = tune_grid,
                        nthread = 1)
  
  print('latitude training model has been done')
  
  l3_train <- data.frame(l2_train, "latPred" = Xgb_latitude[["pred"]][order(Xgb_latitude$pred$rowIndex),"pred" ])
  
  # prediction model for longitude
  Xgb_longitude <- train(x = l3_train ,y = training[,"adj_Long."],
                         method = "xgbTree",
                         trControl = trControl,
                         tuneGrid = tune_grid,
                         nthread = 1)
  
  print('longtitude training model has been done')
  
  
  ##### prediction part ___________________
  
  # predict continent
  regProbs <- predict(Xgb_region, newdata = testing[,variables],type ="prob")
  print('continent prediction has been done')
  
  l1_test <- data.frame(testing[,variables], regProbs)
  
  # predict continent subregion
  classPred <- predict(Xgb_class, newdata = l1_test)
  classProbs <- predict(Xgb_class, newdata = l1_test,type ="prob")
  print('continent detail prediction has been done')
  
  l2_test <-  data.frame(l1_test, classProbs) 

  ##### predict country
  # countryPred <- predict(Xgb_country, newdata = l2_test)
  # countryProbs <- predict(Xgb_country, newdata = l2_test,type ="prob")
  # print('country prediction has been done')
  # 
  # l22_test <-  data.frame(l2_test, countryProbs)
  # 
  ##### predict age (date mean in BP)
  # datePred <- predict(Xgb_date, newdata = l2_test)
  # print('date prediction has been done')
  # date_test <- data.frame(l2_test, datePred)
  
  # predict latitude
  latPred <- predict(Xgb_latitude, newdata = l2_test)
  print('latitude prediction has been done')
  
  l3_test <- data.frame(date_test, latPred)
  
  # predict longitude
  longPred <- predict(Xgb_longitude, newdata = l3_test)
  print('longtitude prediction has been done')
  
  return(list(classPred, latPred, longPred))
}


# Prediction of original coordinates for ancient individuals -------------------
set.seed(18)
trainFolds <-  createFolds(Dataset$Continent_detail, k = 5, returnTrain = T)
GeoPreds <- list()
registerDoParallel(7) 

for (i in 1:5){ 
  
  print(i)
  print('calculating ..............................................')
  
  train <- Dataset[trainFolds[[i]],] # select 4/5 individuals to model training part
  test <- Dataset[-trainFolds[[i]],] # exclude 1/5 individuals form model training part and predict thier origin
  
  # SMOTE: increase the number of individuals in region with a small sample size
  for (j in (1:18)){
    
    name <- as.data.frame(table(train$Continent_detail))
    
    if (min(name$Freq) < 350){
      min_continent <- name[name$Freq == min(name$Freq),]$Var1
      N <- sum(train$Continent_detail == min_continent)
      
      if (min(name$Freq) < 100){
        a <- 600 %/% N * 100
      }
      if (min(name$Freq) > 100 & min(name$Freq) < 300){
        a <- 500 %/% N * 100
      }
      if (min(name$Freq) > 300 & min(name$Freq) < 380){
        a <- 600 %/% N * 100
      }
      
      b <- 10000 * (length(train$Continent_detail) - N) / a / N 
      train2 <- SMOTE(Continent_detail~.,train,perc.over=a, perc.under=b)
      new_ind <- train2[train2$Continent_detail == min_continent, ]
      train_other <- train[train$Continent_detail != min_continent, ]
      train <- rbind(train_other,new_ind)
    }
  }
  print('SMOTE have been done !')
  
  # predict origin of individuals through constructed model 
  testPreds <- ML_Location(training = train, testing = test, classTarget = "Continent_detail",variables = optimumVars)
  GeoPreds[[i]] <- testPreds
  
  print(i) 
  print("finished !!!!!! -----------------------------------------")
}

# merge the prediction results of all individuals
add_preds <- list()
for (i in 1:5){
  add_preds[[i]] <- cbind(Dataset[-trainFolds[[i]],] , 
                          "regionPred"= GeoPreds[[i]][[1]],
                          # "datePred" = GeoPreds[[i]][[2]],
                          "latPred" = GeoPreds[[i]][[2]], 
                          "longPred" = GeoPreds[[i]][[3]] )
}

DataPreds <- rbind.fill(add_preds)

# adjust longitude to right expression
for (i in (1:length(DataPreds$longPred))){
  if (DataPreds[i,"longPred"] > 180){
    DataPreds[i,"adj_longPred"] <- DataPreds[i,"longPred"] - 360
  }else{
    DataPreds[i,"adj_longPred"] <- DataPreds[i,"longPred"]
  }
}

write.csv(DataPreds,"36_DataPreds.csv",row.names = T)

# coast adjustment
# get world coastlines
coastlines <- cbind("x"  = SpatialLines2map(coastsCoarse)$x ,"y" =SpatialLines2map(coastsCoarse)$y)
coastlines <- coastlines[complete.cases(coastlines),]
coastlines <- coastlines[coastlines[,1] < 180 ,]

# Function that can find the nearest land point for individuals whose predicted location is in sea
find_coast <- function(long,lat){
  distances_from_coastline <-  spDistsN1(coastlines , c(long,lat), longlat = TRUE)
  closest_point <-  which.min(distances_from_coastline)
  new_coords <- coastlines[closest_point,]
  return(new_coords)
}

# locate predicted coordinates of individuals into country
locate_country <- map.where(database = "world", DataPreds$adj_longPred, DataPreds$latPred)
adjust_land <- DataPreds[which(is.na(locate_country)),]

# adjust location in sea to nearest land
adjust_coordinate <- mapply(find_coast, long = adjust_land$adj_longPred, lat = adjust_land$latPred )
DataPreds[which(is.na(locate_country)), "latPred"] <- adjust_coordinate[2,]
DataPreds[which(is.na(locate_country)), "adj_longPred"] <- adjust_coordinate[1,]

write.csv(DataPreds,"36_DataPreds_adjustcoast.csv",row.names = T)


# result visualization ---------------------------------------------------------

# plot excavation location of individuals on world map

palette <-c( "brown","red3","maroon1","slateblue1","olivedrab1","darkorchid4","orange","cyan4","cyan2","mediumorchid1","seashell2","dodgerblue3","violetred3","palevioletred1","rosybrown1","mediumspringgreen","gold2","orangered")

map <- getMap(resolution = "coarse")

size1 <- c()
for (i in 1: length(levels(DataPreds$Continent_detail))){
  this_continent <- levels(DataPreds$Continent_detail)[i]
  size1[i] <- length(which(DataPreds$Continent_detail == this_continent))
}

png("origin_global_map.png", width = 13,height = 8, units = 'in', res = 600)
plot(map,xlim = c(-160,160), col = "grey", border = "darkgrey", bg = "lightskyblue1", xlab = "", ylab = "")
title(ylab="Latitude", mgp=c(2,1,0),cex.lab=1.2)

for (i in 1:length(levels(DataPreds$Continent_detail))){
  this_continent <- levels(DataPreds$Continent_detail)[i]
  find_lats <- DataPreds[DataPreds[,"Continent_detail"] == this_continent,]$Lat.
  find_longs <- DataPreds[DataPreds[,"Continent_detail"] == this_continent,]$Long.
  points(find_longs, find_lats, col = palette[i], pch = "+", cex = 1.1)
}

legend(-165,20,legend=c(legend=c(paste0(levels(DataPreds$Continent_detail),"  (",size1,")"))),
       col=palette,pch = "+",cex=1, bg = "lightskyblue1")  

map.axes()
dev.off()

# plot prediction location of all individuals on world map

png("prediction_global_map.png", width = 13,height = 8, units = 'in', res = 600)
plot(map,xlim = c(-160,160), col = "grey", border = "darkgrey", bg = "lightskyblue1", xlab = "", ylab = "")
title(ylab="Latitude",xlab = "Longitude", mgp=c(2,1,0),cex.lab=1.2)

for (i in 1:length(levels(DataPreds$Continent_detail))){
  this_continent <- levels(DataPreds$Continent_detail)[i]
  find_lats <- DataPreds[DataPreds[,"Continent_detail"] == this_continent,]$latPred
  find_longs <- DataPreds[DataPreds[,"Continent_detail"] == this_continent,]$adj_longPred
  points(find_longs, find_lats, col = palette[i], pch = "+", cex = 1.1)
}

legend(-165,20,legend=c(paste0(levels(DataPreds$Continent_detail),"  (",size1,")")),col=palette,pch = "+",cex=1, bg = "lightskyblue1")  

map.axes()
dev.off()


# distance between excavation coordinates and prediction coordinates

# continent subregions as unit

for (i in 1:length((DataPreds$Continent_detail))){ 
  region_lats <- DataPreds[i,]$Lat.
  region_longs <- DataPreds[i,]$Long.
  pred_lat <- DataPreds[i,]$latPred
  pred_long  <- DataPreds[i,]$adj_longPred
  distance <- c()
  for (n in 1:length(region_lats)){ 
    distance[n] <- distHaversine(c(pred_long ,pred_lat ),c(region_longs[n],region_lats[n]))/1000
  }
  DataPreds[i,"distance_from_continent"] <- min(distance, na.rm = TRUE)
}

bar_df <- data.frame(row.names = c( "Overall",levels(DataPreds$Continent_detail)))

for (i in 1: length(levels(DataPreds$Continent_detail))){
  overall_prop <- mean(DataPreds[,"distance_from_continent"] < 200)
  bar_df[1,"0 - 200km"] <- overall_prop
  
  this_continent <- levels(DataPreds$Continent_detail)[i]
  prop <- mean(DataPreds[DataPreds$Continent_detail == this_continent,][,"distance_from_continent"] < 200)
  bar_df[i+1,"0 - 200km"] <- prop
}

for (i in 1: length(levels(DataPreds$Continent_detail))){
  this_continent <- levels(DataPreds$Continent_detail)[i]
  prop <- mean(DataPreds[DataPreds$Continent_detail == this_continent,][,"distance_from_continent"] > 200 & DataPreds[DataPreds$Continent_detail == this_continent,][,"distance_from_continent"] < 500)
  bar_df[i+1,"200 - 500km"] <- prop
  
  overall_prop <- mean(DataPreds[,"distance_from_continent"] > 200 & DataPreds[,"distance_from_continent"] < 500)
  bar_df[ 1,"200 - 500km"] <- overall_prop
}

for (i in 1: length(levels(DataPreds$Continent_detail))){
  this_continent <- levels(DataPreds$Continent_detail)[i]
  prop <- mean(DataPreds[DataPreds$Continent_detail == this_continent,][,"distance_from_continent"] > 500 & DataPreds[DataPreds$Continent_detail == this_continent,][,"distance_from_continent"] < 1000)
  bar_df[i+1,"500 - 1000km"] <- prop
  
  overall_prop <- mean(DataPreds[,"distance_from_continent"] > 500 & DataPreds[,"distance_from_continent"] < 1000)
  bar_df[ 1,"500 - 1000km"] <- overall_prop
}

for (i in 1: length(levels(DataPreds$Continent_detail))){
  this_continent<- levels(DataPreds$Continent_detail)[i]
  prop <- mean(DataPreds[DataPreds$Continent_detail == this_continent,][,"distance_from_continent"] > 1000 & DataPreds[DataPreds$Continent_detail == this_continent,][,"distance_from_continent"] < 2000)
  bar_df[i+1,"1000 - 2000km"] <- prop
  
  overall_prop <- mean(DataPreds[,"distance_from_continent"] > 1000 & DataPreds[,"distance_from_continent"] < 2000)
  bar_df[ 1,"1000 - 2000km"] <- overall_prop
}

for (i in 1: length(levels(DataPreds$Continent_detail))){
  this_continent <- levels(DataPreds$Continent_detail)[i]
  prop <- mean(DataPreds[DataPreds$Continent_detail == this_continent,][,"distance_from_continent"] > 2000 & DataPreds[DataPreds$Continent_detail == this_continent,][,"distance_from_continent"] < 3000)
  bar_df[i+1,"2000 - 3000km"] <- prop
  
  overall_prop <- mean(DataPreds[,"distance_from_continent"] > 2000 & DataPreds[,"distance_from_continent"] < 3000)
  bar_df[ 1,"2000 - 3000km"] <- overall_prop
}

for (i in 1: length(levels(DataPreds$Continent_detail))){
  this_continent <- levels(DataPreds$Continent_detail)[i]
  prop <- mean(DataPreds[DataPreds$Continent_detail == this_continent,][,"distance_from_continent"] > 3000 )
  bar_df[i+1,"> 3000km"] <- prop
  
  overall_prop <- mean(DataPreds[,"distance_from_continent"] > 3000)
  bar_df[ 1,"> 3000km"] <- overall_prop
}

  ## change the order of columns

bar_df2 <- bar_df[c("Overall", "Central.Africa", "North.Africa", "East.Africa", "South.Africa", "Middle.America", "North.America", "South.America", "Central.Asia", "North.Asia", "South.Asia", "East.Asia", "Southeast.Asia", "West.Asia", "North.Europe", "South.Europe", "East.Europe", "West.Europe", "Oceania"),]

order <- c("Overall", "Central.Africa", "North.Africa", "East.Africa", "South.Africa", "Middle.America", "North.America", "South.America", "Central.Asia", "North.Asia", "South.Asia", "East.Asia", "Southeast.Asia", "West.Asia", "North.Europe", "South.Europe", "East.Europe", "West.Europe", "Oceania")

size1 <- c()
size1[1] <- length(DataPreds$Continent)
for (i in 2: length(order)){
  this_continent <- order[i]
  size1[i] <- length(which(DataPreds$Continent_detail == this_continent))
}

par(xpd = T, mar = par()$mar + c(1,0,0,7))
bp <- barplot(t(bar_df2*100), col=c("lightcyan","slategray1","lightblue", "skyblue", "royalblue3", "darkblue"), names.arg=c(paste0(order,"  (",size1,")")) ,args.legend = list(x = "topright", inset=c(-0.5,0)), las =2, cex.names=.7,ylab = "Proportion of sample predictions %")
legend("topright",inset = c(-0.2,0.4), rev(c(colnames(bar_df2))), fill = rev(c("lightcyan","slategray1","lightblue", "skyblue", "royalblue3", "darkblue")) , bty = 1, cex = 0.8)
par(mar=c(5, 4, 4, 2) + 0.1)

png("bar_distance_region.png", width = 13,height = 8, units = 'in', res = 600)

par(xpd = T, mar = par()$mar + c(1,0,0,7), mgp = c(0,0.7,0), las=2)
bp <- barplot(t(bar_df2*100), col=c("lightcyan","slategray1","lightblue", "skyblue", "royalblue3", "darkblue"), names.arg=c(paste0(order,"  (",size1,")")) , args.legend = list(x = "topright", inset=c(-0.5,0)), las =2, cex.names=.6,ylab = "", axisnames = F,axes = F, space =0)
axis(side =2, pos = 0)
mtext(text = c(paste0(order,"  (",size1,")")) , side = 1, at = bp, line = 0, padj = 1, cex = 0.7)
title(ylab="Proportion of sample predictions %", mgp=c(0,0,0),cex.lab=1.2)
legend("topright",inset = c(-0.12,0.4), rev(c(colnames(bar_df2))), fill = rev(c("lightcyan","slategray1","lightblue", "skyblue", "royalblue3", "darkblue")) , bty = 1, cex = 1)

par(mar=c(5, 4, 4, 2) + 0.1)
dev.off()


# country as unit

for (i in 1:length((DataPreds$Country))){

  country_lats <- DataPreds[i,]$Lat.
  country_longs <- DataPreds[i,]$Long.
  pred_lat <- DataPreds[i,]$latPred
  pred_long  <- DataPreds[i,]$adj_longPred
  distance <- c()
  for (n in 1:length(country_lats)){
    distance[n] <- distHaversine(c(pred_long ,pred_lat ),c(country_longs[n],country_lats[n]))/1000
  }
  DataPreds[i,"distance_from_country"] <- min(distance, na.rm = TRUE)
}

bar_df <- data.frame(row.names = c( "Overall_country",levels(DataPreds$Country)))

for (i in 1: length(levels(DataPreds$Country))){
  overall_prop <- mean(DataPreds[,"distance_from_country"] < 200)
  bar_df[1,"0 - 200km"] <- overall_prop

  this_country <- levels(DataPreds$Country)[i]
  prop <- mean(DataPreds[DataPreds$Country == this_country,][,"distance_from_country"] < 200)
  bar_df[i+1,"0 - 200km"] <- prop
}

for (i in 1: length(levels(DataPreds$Country))){
  this_country <- levels(DataPreds$Country)[i]
  prop <- mean(DataPreds[DataPreds$Country == this_country,][,"distance_from_country"] > 200 & DataPreds[DataPreds$Country == this_country,][,"distance_from_country"] < 500)
  bar_df[i+1,"200 - 500km"] <- prop

  overall_prop <- mean(DataPreds[,"distance_from_country"] > 200 & DataPreds[,"distance_from_country"] < 500)
  bar_df[ 1,"200 - 500km"] <- overall_prop
}


for (i in 1: length(levels(DataPreds$Country))){
  this_country <- levels(DataPreds$Country)[i]
  prop <- mean(DataPreds[DataPreds$Country == this_country,][,"distance_from_country"] > 500 & DataPreds[DataPreds$Country == this_country,][,"distance_from_country"] < 1000)
  bar_df[i+1,"500 - 1000km"] <- prop

  overall_prop <- mean(DataPreds[,"distance_from_country"] > 500 & DataPreds[,"distance_from_country"] < 1000)
  bar_df[ 1,"500 - 1000km"] <- overall_prop
}

for (i in 1: length(levels(DataPreds$Country))){
  this_country<- levels(DataPreds$Country)[i]
  prop <- mean(DataPreds[DataPreds$Country == this_country,][,"distance_from_country"] > 1000 & DataPreds[DataPreds$Country == this_country,][,"distance_from_country"] < 2000)
  bar_df[i+1,"1000 - 2000km"] <- prop

  overall_prop <- mean(DataPreds[,"distance_from_country"] > 1000 & DataPreds[,"distance_from_country"] < 2000)
  bar_df[ 1,"1000 - 2000km"] <- overall_prop
}

for (i in 1: length(levels(DataPreds$Country))){
  this_country <- levels(DataPreds$Country)[i]
  prop <- mean(DataPreds[DataPreds$Country == this_country,][,"distance_from_country"] > 2000 & DataPreds[DataPreds$Country == this_country,][,"distance_from_country"] < 3000)
  bar_df[i+1,"2000 - 3000km"] <- prop

  overall_prop <- mean(DataPreds[,"distance_from_country"] > 2000 & DataPreds[,"distance_from_country"] < 3000)
  bar_df[ 1,"2000 - 3000km"] <- overall_prop
}

for (i in 1: length(levels(DataPreds$Country))){
  this_country <- levels(DataPreds$Country)[i]
  prop <- mean(DataPreds[DataPreds$Country == this_country,][,"distance_from_country"] > 3000 )
  bar_df[i+1,"> 3000km"] <- prop

  overall_prop <- mean(DataPreds[,"distance_from_country"] > 3000)
  bar_df[ 1,"> 3000km"] <- overall_prop
}

size1 <- c()
for (i in 1: length(levels(DataPreds$Country))){

  this_country <- levels(DataPreds$Country)[i]
  size1[i] <- length(which(DataPreds$Country == this_country))
}

par(xpd = T, mar = par()$mar + c(1,0,0,7))
bp <- barplot(t(bar_df*100), col=c("lightcyan","slategray1","lightblue", "skyblue", "royalblue3", "darkblue"), names.arg=c("Overall",paste0(levels(DataPreds$Country),"  (",size1,")")) ,args.legend = list(x = "topright", inset=c(-0.5,0)), las =2, cex.names=.6,ylab = "Proportion of sample predictions %")
legend("topright",inset = c(-0.15,0.4), rev(c(colnames(bar_df))), fill = rev(c("lightcyan","slategray1","lightblue", "skyblue", "royalblue3", "darkblue")) , bty = 1, cex = 0.8)

par(mar=c(5, 4, 4, 2) + 0.1)

png("36_date_SMOTE_bar_distance_country.png", width = 13,height = 8, units = 'in', res = 600)

par(xpd = T, mar = par()$mar + c(1,0,0,7), mgp = c(0,0.7,0), las=2)
bp <- barplot(t(bar_df*100), col=c("lightcyan","slategray1","lightblue", "skyblue", "royalblue3", "darkblue"), names.arg=c("Overall",paste0(levels(DataPreds$Country),"  (",size1,")")),  args.legend = list(x = "topright", inset=c(-0.5,0)), las =2, cex.names=.6,ylab = "", axisnames = F,axes = F, space =0)

axis(side =2, pos = 0)
mtext(text = c("Overall",paste0(levels(DataPreds$Country),"  (",size1,")")) , side = 1, at = bp, line = 0, padj = 1, cex = 0.35)
title(ylab="Proportion of sample predictions %", mgp=c(0,0,0),cex.lab=1.2)
legend("topright",inset = c(-0.12,0.4), rev(c(colnames(bar_df))), fill = rev(c("lightcyan","slategray1","lightblue", "skyblue", "royalblue3", "darkblue")) , bty = 1, cex = 1)

par(mar=c(5, 4, 4, 2) + 0.1)
dev.off()


# If age (date mean in BP) label has been added in prediction model

# region age distance

for (i in 1:length(DataPreds$Date.mean.in.BP)){
  real_date <- DataPreds[i,]$Date.mean.in.BP
  pred_date <- DataPreds[i,]$datePred
  DataPreds[i,"date_distance"] <- abs(pred_date - real_date)
}

bar_df <- data.frame(row.names = c(order))

for (i in 2: length(order)){
  overall_prop <- mean(DataPreds[,"date_distance"] < 200)
  bar_df[1,"0 - 200 BP"] <- overall_prop

  this_continent <- order[i]
  prop <- mean(DataPreds[DataPreds$Continent_detail == this_continent,][,"date_distance"] < 200)
  bar_df[i,"0 - 200 BP"] <- prop
}

for (i in 2: length(order)){
  this_continent <- order[i]
  prop <- mean(DataPreds[DataPreds$Continent_detail == this_continent,][,"date_distance"] > 200 & DataPreds[DataPreds$Continent_detail == this_continent,][,"date_distance"] < 500)
  bar_df[i,"200 - 500 BP"] <- prop

  overall_prop <- mean(DataPreds[,"date_distance"] > 200 & DataPreds[,"date_distance"] < 500)
  bar_df[ 1,"200 - 500 BP"] <- overall_prop
}

for (i in 2: length(order)){
  this_continent <- order[i]
  prop <- mean(DataPreds[DataPreds$Continent_detail == this_continent,][,"date_distance"] > 500 & DataPreds[DataPreds$Continent_detail == this_continent,][,"date_distance"] < 1000)
  bar_df[i,"500 - 1000 BP"] <- prop

  overall_prop <- mean(DataPreds[,"date_distance"] > 500 & DataPreds[,"date_distance"] < 1000)
  bar_df[ 1,"500 - 1000 BP"] <- overall_prop
}

for (i in 2: length(order)){
  this_continent<- order[i]
  prop <- mean(DataPreds[DataPreds$Continent_detail == this_continent,][,"date_distance"] > 1000 & DataPreds[DataPreds$Continent_detail == this_continent,][,"date_distance"] < 2000)
  bar_df[i,"1000 - 2000 BP"] <- prop

  overall_prop <- mean(DataPreds[,"date_distance"] > 1000 & DataPreds[,"date_distance"] < 2000)
  bar_df[ 1,"1000 - 2000 BP"] <- overall_prop
}

for (i in 2: length(order)){
  this_continent <- order[i]
  prop <- mean(DataPreds[DataPreds$Continent_detail == this_continent,][,"date_distance"] > 2000 & DataPreds[DataPreds$Continent_detail == this_continent,][,"date_distance"] < 3000)
  bar_df[i,"2000 - 3000 BP"] <- prop

  overall_prop <- mean(DataPreds[,"date_distance"] > 2000 & DataPreds[,"date_distance"] < 3000)
  bar_df[ 1,"2000 - 3000 BP"] <- overall_prop
}

for (i in 2: length(order)){
  this_continent <- order[i]
  prop <- mean(DataPreds[DataPreds$Continent_detail == this_continent,][,"date_distance"] > 3000 )
  bar_df[i,"> 3000 BP"] <- prop

  overall_prop <- mean(DataPreds[,"date_distance"] > 3000)
  bar_df[ 1,"> 3000 BP"] <- overall_prop
}

bar_df2 <- bar_df[c("Overall", "Central.Africa", "North.Africa", "East.Africa", "South.Africa", "Middle.America", "North.America", "South.America", "Central.Asia", "North.Asia", "South.Asia", "East.Asia", "Southeast.Asia", "West.Asia", "North.Europe", "South.Europe", "East.Europe", "West.Europe", "Oceania"),]

size1 <- c()
size1[1] <- length(DataPreds$Continent)
for (i in 2: length(order)){
  this_continent <- order[i]
  size1[i] <- length(which(DataPreds$Continent_detail == this_continent))
}

par(xpd = T, mar = par()$mar + c(1,0,0,7))

bp <- barplot(t(bar_df2*100), col=c("lightcyan","slategray1","lightblue", "skyblue", "royalblue3", "darkblue"), names.arg=c(paste0(order,"  (",size1,")")) ,args.legend = list(x = "topright", inset=c(-0.5,0)), las =2, cex.names=.7,ylab = "Proportion of sample predictions %")
legend("topright",inset = c(-0.2,0.4), rev(c(colnames(bar_df2))), fill = rev(c("lightcyan","slategray1","lightblue", "skyblue", "royalblue3", "darkblue")) , bty = 1, cex = 0.8)

par(mar=c(5, 4, 4, 2) + 0.1)

png("age_bar_region.png", width = 13,height = 8, units = 'in', res = 600)

par(xpd = T, mar = par()$mar + c(1,0,0,7), mgp = c(0,0.7,0), las=2)
bp <- barplot(t(bar_df2*100), col=c("lightcyan","slategray1","lightblue", "skyblue", "royalblue3", "darkblue"), names.arg=c(paste0(order,"  (",size1,")")), args.legend = list(x = "topright", inset=c(-0.5,0)), las =2, cex.names=.6,ylab = "", axisnames = F,axes = F, space =0)

axis(side =2, pos = 0)
mtext(text = c(paste0(order,"  (",size1,")")) , side = 1, at = bp, line = 0, padj = 1, cex = 0.7)
title(ylab="Proportion of sample predictions %", mgp=c(0,0,0),cex.lab=1.2)
legend("topright",inset = c(-0.15,0.4), rev(c(colnames(bar_df2))), fill = rev(c("lightcyan","slategray1","lightblue", "skyblue", "royalblue3", "darkblue")) , bty = 1, cex = 1)

par(mar=c(5, 4, 4, 2) + 0.1)
dev.off()
```

### 2.3 Detect impact of  SNP number on the accuracy of location prediction
Use the genetic information of ancient individuals with less than 10K SNPs excluded from the gene pool to predict their original location. Separate SNP number into different level and then plot the distance of prediction location and excavation location.  
Extract the name list of individuals with SNPs < 10K and add a column named SNP to record SNP number level of each individual. Repeat 1.4.4 to create admixture patterns of remaining ancient people.  
After defining model from last algorithm, predict coordinates of remaining individuals and  plot the accuracy bar plot as before in r.  
```R
Remaining <- read.csv(file='removed_ancient.csv',header=TRUE)

sapply(Remaining,class)
Remaining$Continent_detail <- factor(make.names(Remaining$Continent_detail))
Remaining$Continent <- factor(make.names(Remaining$Continent))
Remaining$SNP <- factor(make.names(Remaining$SNP))

# longitude adjustment
for (i in (1:length(Remaining$Long.))){
  if (Remaining[i,"Continent"] == "Asia" & Remaining[i,"Long."] < 0 ){
    Remaining[i,"adj_Long."] <- 360 + Remaining[i,"Long."]
  }else{
    Remaining[i,"adj_Long."] <- Remaining[i,"Long."]
  }
  if (Remaining[i,"Continent"] == "Oceania" & Remaining[i,"Long."] < 0){
    Remaining[i,"adj_Long."] <- 360 + Remaining[i,"Long."]
  }
}

train <- Dataset
test <- Remaining

for (j in (1:18)){
  
  name <- as.data.frame(table(train$Continent_detail))
  
  if (min(name$Freq) < 350){
    min_continent <- name[name$Freq == min(name$Freq),]$Var1
    N <- sum(train$Continent_detail == min_continent)
    
    if (min(name$Freq) < 100){
      a <- 600 %/% N * 100
    }
    if (min(name$Freq) > 100 & min(name$Freq) < 300){
      a <- 500 %/% N * 100
    }
    if (min(name$Freq) > 300 & min(name$Freq) < 380){
      a <- 600 %/% N * 100
    }
    
    b <- 10000 * (length(train$Continent_detail) - N) / a / N 
    train2 <- SMOTE(Continent_detail~.,train,perc.over=a, perc.under=b)
    new_ind <- train2[train2$Continent_detail == min_continent, ]
    train_other <- train[train$Continent_detail != min_continent, ]
    train <- rbind(train_other,new_ind)
  }
}

testPreds <- ML_Location(training = train, testing = test, classTarget = "Continent_detail",variables = optimumVars)
GeoPreds <- testPreds

add_preds <- list()
add_preds <- cbind(Remaining, 
                   "countryPred"= GeoPreds[[1]],
                   # "datePred" = GeoPreds[[i]][[2]],
                   "latPred" = GeoPreds[[2]], 
                   "longPred" = GeoPreds[[3]] )

DataPreds <- rbind.fill(add_preds)

for (i in (1:length(DataPreds$longPred))){
  if (DataPreds[i,"longPred"] > 180){
    DataPreds[i,"adj_longPred"] <- DataPreds[i,"longPred"] - 360
  }else{
    DataPreds[i,"adj_longPred"] <- DataPreds[i,"longPred"]
  }
}
```

### 2.4 Detection on gene distance between continent subregions
Cluster the continent subregions. After calculating the Euclidean distance of each pair of individuals between each two continent regions, the average Euclidean distance between these two continent regions is obtained. Draw the heat map of distance between subregions.  
Before do the analysis, the input file should be processed: keep only column continent_detail and all the components values.  
```R
library(stats)
library(fastcluster)
library(gplots)

Dataset <- read.csv(file="cluster_components_results.csv",header=TRUE,encoding = 'UTF-8')
Dataset$X <- as.factor(make.names(Dataset$X))

df <- data.frame(row.names = c(levels(Dataset$X))) 

for (i in 1: (length(levels(Dataset$X))-1)){
  continent1 <- levels(Dataset$X)[i]
  continent1_db <- Dataset[Dataset$X == continent1,]
  continent1_db <- continent1_db[,-1]
  
  print(paste0("continent1 is: ", continent1,"--------------------------------"))
  
  for (j in (i+1): (length(levels(Dataset$X)))){
    continent2 <- levels(Dataset$X)[j]
    continent2_db <- Dataset[Dataset$X == continent2,]
    continent2_db <- continent2_db[,-1]
    
    print(paste0("continent2 is: ", continent2))
    
    distance_list <- c()
    
    for (k in 1: dim(continent1_db)[1]){
      for (n in 1: dim(continent2_db)[1]){
        group <- rbind(continent1_db[k,], continent2_db[n,])
        distance <- dist(group,method="euclidean")
        distance_list <- c(distance_list,distance[1])
      }
    }
    
    distance_avg <- mean(distance_list)
    df[continent1,continent1] <- NA
    df[continent2,continent1] <- distance_avg
    print(paste0(continent1," and ", continent2, " have been done"))
  }
}

write.csv(df, "region_distance.csv")

df_dist <- as.dist(df)
out.hclust <- hclust(df_dist,method="complete")

png("tree_cluster.png", width = 10,height = 8, units = 'in', res = 600)
plclust(out.hclust) 
dev.off()

df_matrix <- as.matrix(df_dist)

png("heatmap_cluster.png", width = 10,height = 8, units = 'in', res = 600)
heatmap.2(df_matrix, scale = "none", col=bluered(100), 
          trace = "none", density.info = "histogram",
          hclustfun = function(c)hclust(c,method="complete"),
          keysize = 1.2, cexRow=0.7, cexCol = 0.7,
          srtCol=45,  adjCol = c(1,1))
dev.off()
```

## 3 Locate Origins of ancient people with *Yersinia pestis*

### 3.1 Predict coordanates of ancient people with *Yersinia pestis*
After getting the namelist of ancient people who has been reported to infect *Yersinia pestis*, find its allel frequencies in admixture result table. Extract components pattern of these individuals and then predict thier origin using algorithm built in 2.2.  

### 3.2 Plot original locations of disease ancient people on world map 
Codes are in file Map_route.r  
```R
library(RColorBrewer)
library(rworldmap) 

Dataset <- read.csv(file="predicted_location.csv",header=TRUE,encoding = 'UTF-8')

palette <-c( "brown","red3","maroon1","slateblue1","olivedrab1","cyan2","orangered","cyan4","darkorchid4","green4","seashell2","dodgerblue3","violetred3","aliceblue")
map <- getMap(resolution = "coarse")

png("ind_origin.png", width = 13,height = 8, units = 'in', res = 600)

plot(map,xlim = c(40,70),ylim = c(26,68), col = "grey", border = "darkgrey", bg = "lightskyblue1", xlab = "", ylab = "")

points(Dataset$adj_longPred,Dataset$latPred,col = palette[1:9], pch = "+", cex = 1.5)
text(Dataset$adj_longPred,Dataset$latPred,paste0(Dataset$ID," (",Dataset$Date.mean.in.BP," BP)"), cex=1,pos=4,col=palette[1:9])
points(Dataset$Long,Dataset$Lat,col = palette[1:9], pch = 17, cex = 1)

dev.off()
```