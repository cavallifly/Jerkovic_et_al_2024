# Run clustering
library("ggdendro")
#library("data.table")
#suppressPackageStartupMessages(library(dendextend))
#suppressPackageStartupMessages(library(reshape2))
library(dendextend)
library("ggplot2")
library("gplots")
library("RColorBrewer")
#library("grid")
#library("reshape2")

# Color definition

dat <- read.table("_tmp",header=T,row.names=1)
dat <- as.matrix(dat)
dat <- as.dist(dat, diag=TRUE, upper=TRUE)
#print(dat)

par(cex=0.5)

#hc <- hclust(dat, method="ward.D2")
#plot(hc)

#quit()

dendro <- as.dendrogram(hclust(d = dat, method="ward.D2"))

#labels_colors(dendro) <- colorCodes[groupCodes][order.dendrogram(dendro)]
dendro <- dendro %>%
#      color_branches(k = 100, col = colorCodes[groupCodes][order.dendrogram(dendro)]) %>%
      set("branches_lwd", c(2)) %>%
      set("branches_lty", c(1)) 

# Heatmap
#cols = c("darkred", "red", "white", "blue", "darkblue")
cols = c("darkred", "red", "blue", "darkblue")
mypalette <- colorRampPalette(cols)(length(seq(XXXminXXX,XXXmaxXXX,0.01))-1)
print(length(seq(XXXminXXX,XXXmaxXXX,0.01)))

#layout(matrix(1:1, ncol=1), height=500, width=4, respect=T)
png(file="Rplots.png", width=1300,height=1200)
sM <- 3.5
tM <- (sM)/2	    
#par(mfrow = c(1, 1), mar = c(1, 1, 1, 1))
heatmap.2(as.matrix(dat), breaks=seq(XXXminXXX,XXXmaxXXX,0.01), Rowv=dendro, col = mypalette, symm = TRUE,margins=c(7,14),trace="none",density.info="none", keysize=0.75, key.xlab="1.0-SCC", key.title=NA, cexRow=1.5, labCol = FALSE, key.par = list(cex=1.5))

#ggsave("Rplots.pdf")
dev.off()

