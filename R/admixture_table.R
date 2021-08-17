# This script is used to collect all ancient frequencies with different K value in 
# ADMIXTURE into one table.

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


name <- read_tsv("Americans_header.tsv",col_names = FALSE)# read annotation information
order_name <- name[order(name$X1), ]
header <- order_name[,c(2,3,4,5,6)]


Country_ind <- cbind(header,tb2_2,tb3_2,tb4_2,tb5_2,tb6_2,tb7_2,tb8_2,tb9_2)
colnames(Country_ind) <- c('Continents','Country','Full year','SNPs','ID','a2','b2','a3','b3','c3','a4','b4','c4','d4','a5','b5','c5','d5','e5','a6','b6','c6','d6','e6','f6','a7','b7'
                       ,'c7','d7','e7','f7','g7','a8','b8','c8','d8','e8','f8','g8','h8','a9','b9','c9','d9','e9','f9','g9','h9','i9')


write.csv(Country_ind,"Americans_ADMIXTURE.csv")

