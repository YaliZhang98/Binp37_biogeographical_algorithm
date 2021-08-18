# This script is used to plot coordinates of individual on the world map.

library(RColorBrewer)
library(rworldmap) 

# input file should contain individual ID, predicted longitude, predicted latitude, burial logitude and burial latitude
Dataset <- read.csv(file="predicted_location.csv",header=TRUE,encoding = 'UTF-8')

palette <-c( "brown","red3","maroon1","slateblue1","olivedrab1","cyan2","orangered","cyan4","darkorchid4","green4","seashell2","dodgerblue3","violetred3","aliceblue")
map <- getMap(resolution = "coarse")

png("ind_origin.png", width = 13,height = 8, units = 'in', res = 600)

plot(map,xlim = c(40,70),ylim = c(26,68), col = "grey", border = "darkgrey", bg = "lightskyblue1", xlab = "", ylab = "")

points(Dataset$adj_longPred,Dataset$latPred,col = palette[1:9], pch = "+", cex = 1.5)
text(Dataset$adj_longPred,Dataset$latPred,paste0(Dataset$ID," (",Dataset$Date.mean.in.BP," BP)"), cex=1,pos=4,col=palette[1:9])
points(Dataset$Long,Dataset$Lat,col = palette[1:9], pch = 17, cex = 1)

dev.off()