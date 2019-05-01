dataFile <- commandArgs(6)

## load libraries
library(scales)

## load in the genome file
cat(paste("Reading in:", dataFile,".\n",sep=""))
d <- read.table(dataFile, header=T)
par(pch=16)

half_line <- matrix(c(.6,.4,.4,.6), 2,2) 
full_line_1 <- matrix(c(.33, .17, .5,.5),2,2)
full_line_2 <- matrix(c(.25, .25, .6,.4),2,2)

pdf(file = paste(dataFile,".pdf",sep = ""), width=6, height=4.5, family="ArialMT")
par(pch=16)
curve(1-x, ylab="Z1", xlab="Z0", main="Relatedness gold bars indicate half or \n full sibs AND blue dots are known related")
lines(half_line, col="darkgoldenrod1", cex=5)
lines(full_line_1, col="darkgoldenrod1", cex=5)
lines(full_line_2, col="darkgoldenrod1", cex=5)
with(subset(d,RT=="UN") , points(Z0,Z1,col=alpha(1,.5)))
with(subset(d,RT=="OT") , points(Z0,Z1,col=alpha("blue",.5)))
cat(paste(paste(dataFile,".pdf",sep = "")," was successfully generated.\n",sep=""))
dev.off()
