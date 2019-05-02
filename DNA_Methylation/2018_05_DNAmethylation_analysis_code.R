# Made table with samples included as columns and CGs as rows for all samples that passed QC.

# Of primary table with  2,455,516 total CpG sites, did .na code to remove NAs. Left with 1,433,614 CpGs. File saved as all_meth_nona.txt

library(data.table)
table<- fread("all_meth_nona.txt", header = T)
table <- as.data.frame(table)

# read in clinical metadata for PCA of variance of gene expression explained by each technical and biological covariate

sampleinfo<-read.table("05-2018-HELPTaggingclinicaltable.txt", header=T)

str(sampleinfo)

sample<-sampleinfo[,1]
group<-sampleinfo[,2]
age<-sampleinfo[,3]
sex<-as.factor(sampleinfo[,4])
race<-as.factor(sampleinfo[,5])
glucse<-sampleinfo[,6]
chol<-sampleinfo[,7]
hdl<-sampleinfo[,8]
ldl<-sampleinfo[,9]
trig<-sampleinfo[,10]
insulin<-sampleinfo[,11]
homa<-(sampleinfo[,12])
batchseq<-(sampleinfo[,13])
batchprep<-(sampleinfo[,14])
clust0<-(sampleinfo[,15])
clust1<-(sampleinfo[,16])
clust2<-(sampleinfo[,17])
clust3<-(sampleinfo[,18])
clust4<-(sampleinfo[,20])
clust6<-(sampleinfo[,22])

pca1<-prcomp(table[,-1])


grp.p<-c()
for (i in seq(1:10)){
  group<-sampleinfo[,2]
  grp<-lm(pca1$rotation[,i]~group)
  grp.p<-c(grp.p, anova(grp)$Pr[1])
}

age.p<-c()
for (i in seq(1:10)){
  age<-sampleinfo[,3]
  age<-lm(pca1$rotation[,i]~age)
  age.p<-c(age.p, anova(age)$Pr[1])
}

sex.p<-c()
for (i in seq(1:10)){
  sex<-sampleinfo[,4]
  sex<-lm(pca1$rotation[,i]~sex)
  sex.p<-c(sex.p, anova(sex)$Pr[1])
}

race.p<-c()
for (i in seq(1:10)){
  race<-sampleinfo[,5]
  race<-lm(pca1$rotation[,i]~race)
  race.p<-c(race.p, anova(race)$Pr[1])
}


glucse.p<-c()
for (i in seq(1:10)){
  glucse<-sampleinfo[,6]
  glucse<-lm(pca1$rotation[,i]~glucse)
  glucse.p<-c(glucse.p, anova(glucse)$Pr[1])
}


chol.p<-c()
for (i in seq(1:10)){
  chol<-sampleinfo[,7]
  chol<-lm(pca1$rotation[,i]~chol)
  chol.p<-c(chol.p, anova(chol)$Pr[1])
}

hdl.p<-c()
for (i in seq(1:10)){
  hdl<-sampleinfo[,8]
  hdl<-lm(pca1$rotation[,i]~hdl)
  hdl.p<-c(hdl.p, anova(hdl)$Pr[1])
}

ldl.p<-c()
for (i in seq(1:10)){
  ldl<-sampleinfo[,9]
  ldl<-lm(pca1$rotation[,i]~ldl)
  ldl.p<-c(ldl.p, anova(ldl)$Pr[1])
}

trig.p<-c()
for (i in seq(1:10)){
  trig<-sampleinfo[,10]
  trig<-lm(pca1$rotation[,i]~trig)
  trig.p<-c(trig.p, anova(trig)$Pr[1])
}

insulin.p<-c()
for (i in seq(1:10)){
  insulin<-sampleinfo[,11]
  insulin<-lm(pca1$rotation[,i]~insulin)
  insulin.p<-c(insulin.p, anova(insulin)$Pr[1])
}

homa.p<-c()
for (i in seq(1:10)){
  homa<-sampleinfo[,12]
  homa<-lm(pca1$rotation[,i]~homa)
  homa.p<-c(homa.p, anova(homa)$Pr[1])
}

batchseq.p<-c()
for (i in seq(1:10)){
  batchseq<-sampleinfo[,13]
  batchseq<-lm(pca1$rotation[,i]~batchseq)
  batchseq.p<-c(batchseq.p, anova(batchseq)$Pr[1])
}

batchprep.p<-c()
for (i in seq(1:10)){
  batchprep<-sampleinfo[,14]
  batchprep<-lm(pca1$rotation[,i]~batchprep)
  batchprep.p<-c(batchprep.p, anova(batchprep)$Pr[1])
}

clust0.p<-c()
  for (i in seq(1:10)){
    clust0<-sampleinfo[,15]
    clust0<-lm(pca1$rotation[,i]~clust0)
    clust0.p<-c(clust0.p, anova(clust0)$Pr[1])
  }
  
  clust1.p<-c()
  for (i in seq(1:10)){
    clust1<-sampleinfo[,16]
    clust1<-lm(pca1$rotation[,i]~clust1)
    clust1.p<-c(clust1.p, anova(clust1)$Pr[1])
  }
  
  clust2.p<-c()
  for (i in seq(1:10)){
    clust2<-sampleinfo[,17]
    clust2<-lm(pca1$rotation[,i]~clust2)
    clust2.p<-c(clust2.p, anova(clust2)$Pr[1])
  }
  
  clust4.p<-c()
  for (i in seq(1:10)){
    clust4<-sampleinfo[,20]
    clust4<-lm(pca1$rotation[,i]~clust4)
    clust4.p<-c(clust4.p, anova(clust4)$Pr[1])
  }
  
  clust6.p<-c()
  for (i in seq(1:10)){
    clust6<-sampleinfo[,22]
    clust6<-lm(pca1$rotation[,i]~clust6)
    clust6.p<-c(clust6.p, anova(clust6)$Pr[1])
  }  
  
# make heatmap of p values associated with covariates for each PC

pvals.raw<-rbind(grp.p, age.p, sex.p, race.p, bmi.p, glucse.p, chol.p, hdl.p, ldl.p, trig.p, batch.p, insulin.p, homa.p, leptin.p, adipo.p, batch.p, totalreads.p, prtncoding.p, pctdup.p, clust0.p, clust1.p, clust2.p, clust4.p, clust6.p) 
    
logpvals.raw<- -log10(pvals.raw)
  
library(gplots)
library(RColorBrewer)

hmcol<-colorRampPalette(brewer.pal(10,"RdYlBu"))(256)
hmcol2<-colorRampPalette(brewer.pal(9,"Blues"))(256)
hmcol3<-colorRampPalette(brewer.pal(9,"Blues"))(5)
hmcol4<-hmcol3
hmcol4[1]<-"white"
max_num <- ifelse(max(logpvals.raw3)>4, max(logpvals.raw3), 5)
  
heatmap.2(logpvals.raw,Rowv=F,Colv=colnames(pvals.raw3),dendrogram='none',trace='none',margins=c(8,8),col=hmcol4,colsep=c(1:9),rowsep=c(1:20),sepwidth=c(0.05,0.025), sepcolor="black", breaks=c(0,1.30103,2,3,4,max_num))
 
  
# make bar graph of variance explained by each PC

pca<-princomp(table[,-1])
loadings = pca$loadings[]

p.variance.explained = pca$sdev^2 / sum(pca$sdev^2)

barplot(100*p.variance.explained[1:10], las=2, xlab='', ylab='% Variance Explained', 
        main="Variance explained by each PCA component")
  
 # regression analysis without and with cell proportions- Upto 5 PCs explain 99% variance in DNA methylation, included those variables

i <-1
for( i in 1:nrow(table)){
  y <- table[i,2:99]
  y <-as.numeric(y)
  lmfit <- lm(log(y+1) ~group+age+sex+race+ldl+insulin+batchseq+batchprep) 
  write.table(t(c(as.character(table[i,1]),summary(lmfit)$coefficients[,1],summary(lmfit)$coefficients[,4],1-pf(summary(lmfit)$fstat[1], summary(lmfit)$fstat[2],summary(lmfit)$fstat[3]))), "05-2018-helptaggingmodel_baseline.txt", col.names=F, row.names=F, quote=F, sep="\t", append=T)
}
        
fdrtable<-read.table("05-2018-helptaggingmodel_baseline.txt")

 fdrtable_baseline<-cbind(fdrtable,p.adjust(fdrtable[,19]))

for( i in 1:nrow(table))
{
  y <- table[i,2:104]
  y <-as.numeric(y)
  lmfit <- lm(log(y+1) ~group+age+sex+race+ldl+insulin+batchseq+batchprep+clust4) 
  write.table(t(c(as.character(table[i,1]),summary(lmfit)$coefficients[,1],summary(lmfit)$coefficients[,4],1-pf(summary(lmfit)$fstat[1], summary(lmfit)$fstat[2],summary(lmfit)$fstat[3]))), "05-2018-helptaggingmodel_cellprop.txt", col.names=F, row.names=F, quote=F, sep="\t", append=T) 
}

fdrtable<-read.table("05-2018_helptaggingmodel_cellprop.txt")
head(fdrtable)
fdrtable_cellprop<-cbind(fdrtable,p.adjust(fdrtable[,21]))

# bring in annotations

anno<- read.table("CpGID_Prom_GB_Enh2_annot.txt", header = T)
head(anno)

anno_fdrtable_baseline<- merge(anno, fdrtable_baseline, by="X.tid", all.y=T)

anno_fdrtable_cellprop<- merge(anno, fdrtable_cellprop, by="X.tid", all.y=T)

# Subset based on FDR<0.05 for entire table and between group difference in methylation of p<0.05

fdrtable_sig<- fdrtable_id[ which(fdrtable_id$fdr_p<0.05),]

fdrtable_sig_grp <-fdrtable_sig[ which(fdrtable_sig$p_grp<0.05),]

# calculate difference in percent methylation between groups and subset for percent methylation =>10

a_oa_p=rowttests(percentmeth,classes)
head(a_oa_p)
colnames(a_oa_p)<-c("statistic", "dm", "p.value")

X.tid <-table[,1]
a_oa_p_xtid<-cbind(X.tid, a_oa_p)

fdrtable_sig_grp_meth<- merge(fdrtable_sig_grp, a_oa_p, by="X.tid", all.x=T)

# To find overlap between hcgenes from RNASeq analysis and differential methylation

hcgenes<- read.table("05-2018-hcgenes.txt", header = T)

meth_expressn<-merge(hcgenes, fdrtable_sig_grp_meth, by="Ensemblid", all.x= T)

# to find overlap between meQTL and differentially methylated genes

mqtl<-read.table("mqtl_xtid.txt", header = T)
mqtl_ensembl <-merge(mqtl, anno, by = "X.tid", all.x = T)

mqtl_differentialmeth<-merge(fdrtable_sig_grp_meth, mqtl_ensembl, by = "X.tid", all.x = T)

# to find overlap between meQTL and differentially expressed genes (hcgenes from RNASeq analysis)

mqtl_differentialexp<-merge(hcgenes, mqtl_ensembl, by = "Ensemblid", all.x = T)
