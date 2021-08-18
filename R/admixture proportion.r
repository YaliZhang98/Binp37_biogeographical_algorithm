# This script is used produce a stacked graph for ADMIXTURE proportions in gene pool.

library(RColorBrewer)

Dataset <- read.csv(file="admixture_proportions.csv",header=TRUE) # This file just contain subregion and allele frequencies

palette <- colorRampPalette(brewer.pal(12, "Paired"))(36)


png("admixture_proportion.png", width = 13,height = 8, units = 'in', res = 600)

par(xpd = T, mar = par()$mar + c(1,0,0,7), mgp = c(0,0.7,0), las=2)
bp <- barplot(t(as.matrix(Dataset[,-1:-2])),col = palette,
              ylab = "",
              
              border = NA)

mtext(text = c(Dataset$Continent_detail), side = 1, at = bp, line = 0, padj = 1, cex = 0.5,las = 2)
bp
legend("topright",inset = c(-0.15,0.2),c(test),fill = palette, bty = 1, cex = 0.6)


par(mar=c(5, 4, 4, 2) + 0.1)
dev.off()

