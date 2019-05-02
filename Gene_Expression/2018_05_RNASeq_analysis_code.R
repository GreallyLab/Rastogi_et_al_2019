# starting with table labeled "seq" which included samples in columns and ensembl ID in rows

#generated the normalized counts table

require(DESeq)
cds = newCountDataSet( seq, condition )
cds = estimateSizeFactors( cds )
sizeFactors( cds )
head( counts( cds, normalized=TRUE ))
normalizedcounts <-counts( cds, normalized=TRUE )
rownames(normalizedcounts)=rownames( counts )

write.table(normalizedcounts, file="normcounts_allexcept_AJqc_final",sep="\t",row.names=TRUE,col.names=FALSE,quote=FALSE)

read.table("normcounts_allexcept_AJqc_final") -> table

head(table)
dim(table)

# read in clinical metadata for PCA of variance of gene expression explained by each technical and biological covariate

sampleinfo<-read.table("2018-05-RNAseqclinicaltable_AJqcfinal_cellprop.txt", header=T)
str(sampleinfo)

sample<-sampleinfo[,1]
group<-sampleinfo[,2]
age<-sampleinfo[,3]
sex<-as.factor(sampleinfo[,4])
race<-as.factor(sampleinfo[,5])
bmi<-sampleinfo[,6]
glucse<-sampleinfo[,7]
chol<-sampleinfo[,8]
hdl<-sampleinfo[,9]
ldl<-sampleinfo[,10]
trig<-sampleinfo[,11]
insulin<-sampleinfo[,12]
homa<-sampleinfo[,13]
totalreads<-sampleinfo[,16]
prtncoding<-sampleinfo[,18]
pctdup<-sampleinfo[,19]
batch<-sampleinfo[,20]
clust0<-(sampleinfo[,21])
clust1<-(sampleinfo[,22])
clust2<-(sampleinfo[,23])
clust3<-(sampleinfo[,24])
clust4<-(sampleinfo[,26])
clust5<-(sampleinfo[,27])
clust6<-(sampleinfo[,28])


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

bmi.p<-c()
for (i in seq(1:10)){
  bmi<-sampleinfo[,6]
  bmi<-lm(pca1$rotation[,i]~bmi)
  bmi.p<-c(bmi.p, anova(bmi)$Pr[1])
}

glucse.p<-c()
for (i in seq(1:10)){
  glucse<-sampleinfo[,7]
  glucse<-lm(pca1$rotation[,i]~glucse)
  glucse.p<-c(glucse.p, anova(glucse)$Pr[1])
}

chol.p<-c()
for (i in seq(1:10)){
  chol<-sampleinfo[,8]
  chol<-lm(pca1$rotation[,i]~chol)
  chol.p<-c(chol.p, anova(chol)$Pr[1])
}

hdl.p<-c()
for (i in seq(1:10)){
  hdl<-sampleinfo[,9]
  hdl<-lm(pca1$rotation[,i]~hdl)
  hdl.p<-c(hdl.p, anova(hdl)$Pr[1])
}

ldl.p<-c()
for (i in seq(1:10)){
  ldl<-sampleinfo[,10]
  ldl<-lm(pca1$rotation[,i]~ldl)
  ldl.p<-c(ldl.p, anova(ldl)$Pr[1])
}

trig.p<-c()
for (i in seq(1:10)){
  trig<-sampleinfo[,11]
  trig<-lm(pca1$rotation[,i]~trig)
  trig.p<-c(trig.p, anova(trig)$Pr[1])
}


insulin.p<-c()
for (i in seq(1:10)){
  insulin<-sampleinfo[,12]
  insulin<-lm(pca1$rotation[,i]~insulin)
  insulin.p<-c(insulin.p, anova(insulin)$Pr[1])
}

homa.p<-c()
for (i in seq(1:10)){
  homa<-sampleinfo[,13]
  homa<-lm(pca1$rotation[,i]~homa)
  homa.p<-c(homa.p, anova(homa)$Pr[1])
}


totalreads.p<-c()
for (i in seq(1:10)){
  totalreads<-sampleinfo[,16]
  totalreads<-lm(pca1$rotation[,i]~totalreads)
  totalreads.p<-c(totalreads.p, anova(totalreads)$Pr[1])
}

prtncoding.p<-c()
for (i in seq(1:10)){
  prtncoding<-sampleinfo[,18]
  prtncoding<-lm(pca1$rotation[,i]~prtncoding)
  prtncoding.p<-c(prtncoding.p, anova(prtncoding)$Pr[1])
}

pctdup.p<-c()
for (i in seq(1:10)){
  pctdup<-sampleinfo[,19]
  pctdup<-lm(pca1$rotation[,i]~pctdup)
  pctdup.p<-c(pctdup.p, anova(pctdup)$Pr[1])
}

batch.p<-c()
for (i in seq(1:10)){
  batch<-sampleinfo[,20]
  batch<-lm(pca1$rotation[,i]~batch)
  batch.p<-c(batch.p, anova(batch)$Pr[1])
}

clust0.p<-c()
for (i in seq(1:10)){
  clust0<-sampleinfo[,21]
  clust0<-lm(pca1$rotation[,i]~clust0)
  clust0.p<-c(clust0.p, anova(clust0)$Pr[1])
}

clust1.p<-c()
for (i in seq(1:10)){
  clust1<-sampleinfo[,22]
  clust1<-lm(pca1$rotation[,i]~clust1)
  clust1.p<-c(clust1.p, anova(clust1)$Pr[1])
}


clust2.p<-c()
for (i in seq(1:10)){
  clust2<-sampleinfo[,25]
  clust2<-lm(pca1$rotation[,i]~clust2)
  clust2.p<-c(clust2.p, anova(clust2)$Pr[1])
}

clust4.p<-c()
for (i in seq(1:10)){
  clust4<-sampleinfo[,26]
  clust4<-lm(pca1$rotation[,i]~clust4)
  clust4.p<-c(clust4.p, anova(clust4)$Pr[1])
}

clust6.p<-c()
for (i in seq(1:10)){
  clust6<-sampleinfo[,28]
  clust6<-lm(pca1$rotation[,i]~clust6)
  clust6.p<-c(clust6.p, anova(clust6)$Pr[1])
}

# make heatmap of p values associated with covariates for each PC

pvals.raw<-rbind(grp.p, age.p, sex.p, race.p, bmi.p, glucse.p, chol.p, hdl.p, ldl.p, trig.p, batch.p, insulin.p, homa.p, leptin.p, adipo.p, batch.p, totalreads.p, prtncoding.p, pctdup.p, clust0.p, clust1.p, clust2.p, clust4.p, clust6.p)


write.table(pvals.raw,"pvals_pca_rnaseq_allexceptAJqc_final_allvariables.txt", sep="\t",row.names=TRUE,col.names=TRUE,quote=FALSE)

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


# regression analysis without and with cell proportions- PCA upto 5 explains 99% variance, included those variables

for( i in 1:nrow(table))
{
  y <- table[i,2:104]
  y <-as.numeric(y)
  lmfit <- lm(log(y+1) ~group+age+sex+race+insulin+batch+pctdup+prtncoding+totalreads) 
  write.table(t(c(as.character(table[i,1]),summary(lmfit)$coefficients[,1],summary(lmfit)$coefficients[,4],1-pf(summary(lmfit)$fstat[1], summary(lmfit)$fstat[2],summary(lmfit)$fstat[3]))), "05-2018_rnaseqmodel_AJqcfinal_baseline.txt", col.names=F, row.names=F, quote=F, sep="\t", append=T) 
}

fdrtable<-read.table("05-2018_rnaseqmodel_AJqcfinal_baseline.txt")
head(fdrtable)
fdrtable<-cbind(fdrtable,p.adjust(fdrtable[,22]))
head(fdrtable)
colnames(fdrtable)<-c("id",	"b_int",	"b_grp",	"b_age",	"b_sex",	"b_race",	"b_insulin",	"b_batch",	"b_pctdup",	"b_prtncoding", "b_totalreads", "p_int",	"p_grp",	"p_age",	"p_sex",	"p_race",	"p_insulin",	"p_batch", "p_pctdup", "p_prtncoding", "p_totalreads",	"p_model",	"p_fdr")

write.table(fdrtable,"05_2018_fdrtable_AJqcfinal_baseline.txt", col.names=T, row.names=F, quote=F, sep="\t")

for( i in 1:nrow(table))
{
  y <- table[i,2:104]
  y <-as.numeric(y)
  lmfit <- lm(log(y+1) ~group+age+sex+race+insulin+batch+pctdup+prtncoding+totalreads+clust0+clust1+clust2+clust4) 
  write.table(t(c(as.character(table[i,1]),summary(lmfit)$coefficients[,1],summary(lmfit)$coefficients[,4],1-pf(summary(lmfit)$fstat[1], summary(lmfit)$fstat[2],summary(lmfit)$fstat[3]))), "05-2018_rnaseqmodel_AJqcfinal_cellprop.txt", col.names=F, row.names=F, quote=F, sep="\t", append=T) 
}

fdrtable<-read.table("05-2018_rnaseqmodel_AJqcfinal_cellprop.txt")
head(fdrtable)
fdrtable<-cbind(fdrtable,p.adjust(fdrtable[,28]))
head(fdrtable)
colnames(fdrtable)<-c("id",	"b_int", "b_grp", "b_age", "b_sex", "b_race", "b_insulin",	"b_batch", "b_pctdup", "b_prtncoding", "b_totalreads", "b_clust0", "b_clust2", "b_clust4", "p_int", "p_grp",	"p_age", "p_sex", "p_race",	"p_insulin",	"p_batch", "p_pctdup", "p_prtncoding", "p_totalreads", "p_clust0", "p_clust2", "p_clust4",	"p_model",	"p_fdr")

write.table(fdrtable,"05_2018_fdrtable_AJqcfinal_cellprop.txt", col.names=T, row.names=F, quote=F, sep="\t")

# subsetted genes with FDR<0.05, between group p value<0.05 = hcgenes, collated genes differentially expressed at baseline not influenced by cell subtype proportions, and those identified after cell subtype proportions. 

#subset normcounts to hc genes to compare gene expression between study groups for hcgenes

head(table)

genes_normcounts<-merge(hcgenes, table, by="id", all.x = T)


hcgenes<- read.table("05-2018-hcgenes.txt", header = T)

# overlap hcgenes with eQTLs to identify genes differentially expressed that overlap with eQTL

eqtls<-read.table("08_2018_eqtls.txt", header=F)

eqtl_hc<-merge(hcgenes, eqtlgenes, by="id", all.x = T)
write.table(eqtl_hc,"08_2018_hcgenes_egenes_AJqcfinal.txt", col.names=T, row.names=F, quote=F, sep="\t")

#overlap hcgenes with meQTLs  to identify genes differentially expressed that overlap with meQTL

mqtl<-read.table("08-2018_mqtlensemblmarks.txt", header = T)

mqtl_hc<-merge(hcgenes, mqtlgenes, by= "id", all.x=T)

# correlate pulm fn with gene expression

clininfo<- read.table("08_2018_RNASeq_clininfo_transposed.txt", header = T)

normcounts_clin <-rbind(clininfo, normalizedcounts)

write.table(normcounts_clin,"08_2018_normcounts_clininfo.txt", col.names=T, row.names=T, quote=F, sep="\t")

