 # load libraries
library(data.table)
 # set options
options(scipen = 999, stringsAsFactors = F)
 # read in eQTLs in eGene bins
eQTL_bins <- read.table("eQTL_final3_eQTL_eGene_100kbTSS_1kbwin.txt")
head(eQTL_bins)

eQTL_bins_count <-cbind(table(factor(eQTL_bins[,4], levels=1:200)))
eQTL_bins_seq <-cbind(seq(1,200),eQTL_bins_count)


pdf("TSS_eQTL_eGenes.pdf", width=10, height = 7, family="ArialMT")
plot(1,1,xlim=c(0,200),ylim=c(0,50),type="n",xlab=" 1kb windows flanking the TSS",
     ylab="eQTL variants",axes=FALSE)
axis(side=1,at=c(0,100.5,200),labels=c("-100 Kb","TSS","+100 Kb"))
axis(side=2)
box()
points(rownames(eQTL_bins_seq),eQTL_bins_seq[,2],cex=.5,col="darkgreen")
lines(smooth.spline(rownames(eQTL_bins_seq),eQTL_bins_seq[,2],spar=.2),col="darkgreen")
dev.off()

# read in all variants in all gene bins
var_bins <- fread("allVar_allGene_100kbTSS_1kbwin.txt")
head(var_bins)

var_bins_count <-cbind(table(factor(var_bins$V4, levels=1:200)))
var_bins_seq <-cbind(seq(1,200),var_bins_count)

pdf("TSS_allvar_allG.pdf", width=10, height = 7, family="ArialMT")
plot(1,1,xlim=c(0,200),ylim=c(0,3000),type="n",
     xlab=" 1kb windows flanking the \nTSS of all genes",
     ylab="all variants",axes=FALSE)
axis(side=1,at=c(0,100.5,200),labels=c("-100 Kb","TSS","+100 Kb"))
axis(side=2)
box()
points(rownames(var_bins_seq),var_bins_seq[,2],cex=.5,col="darkgreen")
lines(smooth.spline(rownames(var_bins_seq),var_bins_seq[,2],spar=.2),col="darkgreen")
dev.off()

plot(1,1,xlim=c(80,120),ylim=c(0,3000),type="n",
     xlab=" 1kb windows flanking the \nTSS of all genes",
     ylab="all variants",axes=FALSE)
axis(side=1,at=c(0,100.5,200),labels=c("-100 Kb","TSS","+100 Kb"))
axis(side=2)
box()
points(rownames(var_bins_seq),var_bins_seq[,2],cex=.5,col="darkgreen")
lines(smooth.spline(rownames(var_bins_seq),var_bins_seq[,2],spar=.2),col="darkgreen")

 # What about density from QTL tools - all TSSs
D <- read.table("dens_TSS_eQTL.txt", head=FALSE, stringsAsFactors=FALSE)

pdf("allTSS_dens_eQTL.pdf", width=6, height = 4, family="ArialMT")
plot((D$V1+D$V2)/2, D$V3, type="l", xlab="Distance to QTLs", ylab="#annotations/kb",
     main="all gene TSSs around eQTLs")
dev.off()

 # What about density from QTL tools - only eGene TSSs
D_Egene <- read.table("dens_TSS_eGene_eQTL.txt", head=FALSE, stringsAsFactors=FALSE)

pdf("eGeneTSS_dens_eQTL.pdf", width=6, height = 4, family="ArialMT")
plot((D_Egene$V1+D_Egene$V2)/2, D_Egene$V3, type="l", xlab="Distance to QTLs",
     ylab="#annotations/kb", main = "Egene TSSs around eQTLs")
dev.off()

