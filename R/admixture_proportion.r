Dataset <- read.csv(file="3admixture_proportions.csv",header=TRUE)

png("admixture_proportion.png", width = 13,height = 8, units = 'in', res = 600)

par(xpd = T, mar = par()$mar + c(1,0,0,7), mgp = c(0,0.7,0), las=2)
bp <- barplot(t(as.matrix(Dataset[,-1:-2])),col = c( "brown","red3","maroon1","slateblue1","olivedrab1","darkorchid4","orange","cyan4",
                                                     "cyan2","mediumorchid1","seashell2","dodgerblue3","violetred3","palevioletred1",
                                                     "rosybrown1","mediumspringgreen","gold2","orangered","black"),
              ylab = "",
              
              border = NA)

mtext(text = c(Dataset$Continent_detail), side = 1, at = bp, line = 0, padj = 1, cex = 0.5,las = 2)

par(mar=c(5, 4, 4, 2) + 0.1)
dev.off()
