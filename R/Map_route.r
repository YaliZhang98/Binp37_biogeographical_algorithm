library(RColorBrewer)
library(rworldmap) 


Dataset <- read.csv(file="1sickind_predicted_location.csv",header=TRUE,encoding = 'UTF-8')
up_Dataset <- read.csv(file="up_sickind_predicted_location.csv",header=TRUE,encoding = 'UTF-8')
down_Dataset <- read.csv(file="down_sickind_predicted_location.csv",header=TRUE,encoding = 'UTF-8')


palette <-c( "brown","red3","maroon1","slateblue1","olivedrab1","cyan2","orangered","cyan4","darkorchid4","green4","seashell2","dodgerblue3","violetred3","aliceblue")
map <- getMap(resolution = "coarse")


png("0814compare_ind_point.png", width = 13,height = 8, units = 'in', res = 600)

plot(map,xlim = c(40,70),ylim = c(26,68), col = "grey", border = "darkgrey", bg = "lightskyblue1", xlab = "", ylab = "")

points(Dataset$adj_longPred,Dataset$latPred,col = palette[1:9], pch = "+", cex = 1.5)
text(Dataset$adj_longPred,Dataset$latPred,paste0(Dataset$ID," (",Dataset$Date.mean.in.BP," BP)"), cex=1,pos=4,col=palette[1:9])
points(Dataset$Long,Dataset$Lat,col = palette[1:9], pch = 17, cex = 1)

points(up_Dataset$adj_longPred,up_Dataset$latPred,col = palette[10:12], pch = "+", cex = 1.5)
text(up_Dataset$adj_longPred,up_Dataset$latPred,paste0(up_Dataset$ID," (",up_Dataset$Date.mean.in.BP," BP)"), cex=1,pos=3,col=palette[10:12])
points(up_Dataset$Long,up_Dataset$Lat,col = palette[10:12], pch = 17, cex = 1)

points(down_Dataset$adj_longPred,down_Dataset$latPred,col = palette[13:14], pch = "+", cex = 1.5)
text(down_Dataset$adj_longPred,down_Dataset$latPred,paste0(down_Dataset$ID," (",down_Dataset$Date.mean.in.BP," BP)"), cex=1,pos=1,col=palette[13:14])
points(down_Dataset$Long,down_Dataset$Lat,col = palette[13:14], pch = 17, cex = 1)

dev.off()





