# Find individuals with >10K SNPs

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
# write.table(name, "10K_namelist.tsv", sep="\t", row.names = F)

name_SNP <- tenK_ind[,c(2,7)]
write.table(name_SNP, "10K_name_SNP.tsv", sep="\t", row.names = F)
