dataFile <- commandArgs()[6]
name <- commandArgs()[7]

# load the libraries
library(qvalue)

# Read in the 
cat(paste("Reading in: ", dataFile,".\n",sep=""))
d <- read.table(dataFile, header=F, stringsAsFactors=F)

# make the plot
cat(paste("Generating Beta Approximation Plot",".\n",sep=""))
pdf(file = paste(name,"_betaPlot", ".pdf",sep = ""), width=6,
    height=4.5, family="ArialMT")
plot(d$V18, d$V19, xlab="Direct method", ylab="Beta approximation", main=name)
abline(0, 1, col="red")
dev.off()

d_fil <- d[!is.na(d$V19),]
d_fil$qval <- qvalue(d_fil$V19)$qvalue
sum(d_fil$qval <= 0.05)
cat(paste("Found ", sum(d_fil$qval <= 0.05), " significant QTls",".\n",sep=""))
cat(paste("Writing significant QTls out",".\n",sep=""))
write.table(d_fil[which(d_fil$qval <= 0.05), ],
            paste("Sig",name,".txt", sep=""), quote=FALSE, 
            row.names=FALSE, col.names=FALSE)
