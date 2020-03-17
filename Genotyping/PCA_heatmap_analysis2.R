# let's look at the influences on the quant PCs via metadata 

 # load libraries and set options
library(data.table)
library(RColorBrewer)
library(gplots)
options(scipen = 999, stringsAsFactors = F)

 # load in the quant data
quant <- fread("../quants2/Filtered_matrix_gene_rpkm.bed", header=T)
quant <- as.data.frame(quant)
quant <- quant[,-c(1:6)]
head(quant)
dim(quant) # 12960    91
pca1 <- prcomp(na.omit(quant), scale=F)
head(rownames(pca1$rotation))

options(scipen = 999, stringsAsFactors = F)
gene_exp_PC <- read.table("../quants2/genes_90percent.pca", header=T, 
                          colClasses = c("character", rep("numeric",76)))
head(gene_exp_PC)
str(gene_exp_PC)
dim(gene_exp_PC)
colnames(gene_exp_PC)[-1] <- substr(colnames(gene_exp_PC)[-1], 2,8)
gene_exp_PC_t <- as.data.frame(t(gene_exp_PC))
str(gene_exp_PC_t)
head(gene_exp_PC_t)
colnames(gene_exp_PC_t) <- gene_exp_PC_t[1,]
gene_exp_PC_t <- gene_exp_PC_t[-1,]
head(gene_exp_PC_t)
for(i in 1:ncol(gene_exp_PC_t)){
  gene_exp_PC_t[,i] <- as.numeric(gene_exp_PC_t[,i])
}
class(gene_exp_PC_t$genes_90percent_1_1_svd_PC1)
head(gene_exp_PC_t)
rownames(gene_exp_PC_t)


# read in the sample metadata that Deepa provided
sampleinfo <- read.csv("sample_metadata_eQTL_4.csv", header=T, stringsAsFactors = F)
head(sampleinfo)
sampleinfo[42,]
str(sampleinfo)
dim(sampleinfo) #  91 21

 # making qualitative data factors
sampleinfo[,3] <- as.factor(sampleinfo[,3]) #group
 # sampleinfo[,4] <- as.factor(sampleinfo[,4]) #age
sampleinfo[,5] <- as.factor(sampleinfo[,5]) #sex
sampleinfo[,6] <- as.factor(sampleinfo[,6]) #ethnicity
# sampleinfo[,7] <- as.factor(sampleinfo[,7]) #insulin 
sampleinfo[,15] <- as.factor(sampleinfo[,15]) # seq batch
sampleinfo[,16] <- as.factor(sampleinfo[,16]) # geno batch

sampleinfo <- sampleinfo[,-c(17:21)]
dim(sampleinfo) #  91 16

# I had already switched the RNAseq data for the sample metadata sheet

# ordering the samples to sync with the gene expr data as well as remove uneeded rows 
idx_match_meta_pca <- na.omit(match(do.call(rbind, strsplit(rownames(pca1$rotation), "_"))[,2],
                                    sampleinfo$Sample))
sampleinfo_order <- sampleinfo[idx_match_meta_pca,]
head(sampleinfo_order)
dim(sampleinfo_order) # 91 16
str(sampleinfo_order)

# read in the genotyping PCs
geno_PC <- read.table("../../Genotype2/genotypes.pca.pca", header = T)
head(geno_PC)
colnames(geno_PC)[-1] <- substr(colnames(geno_PC)[-1], 2,8)
geno_PC_t <- as.data.frame(t(geno_PC))
head(geno_PC_t)
colnames(geno_PC_t) <- geno_PC_t[1,]
geno_PC_t <- geno_PC_t[-1,]
for(i in 1:ncol(geno_PC_t)){
  geno_PC_t[,i] <- as.numeric(geno_PC_t[,i])
}
head(geno_PC_t)
class(geno_PC_t$genotypes.pca_1_1_svd_PC1)

colnames(geno_PC_t) <- paste("Geno_",do.call(rbind, strsplit(colnames(geno_PC_t), "_"))[,5], sep="")
rownames(geno_PC_t)
sampleinfo_order$Sample
# create metadata combined with geno PCs
sampleinfo_geno <- cbind(sampleinfo_order, geno_PC_t)
rownames(sampleinfo_geno) <- sampleinfo_geno[,2]
sampleinfo_geno <- sampleinfo_geno[,-c(1:2)]
head(sampleinfo_geno)
ncol(sampleinfo_geno) #105
colnames(sampleinfo_geno)

 # how many expression PCs should be anaylzed?
exp_pc_stats <- read.table("../quants2/genes_90percent.pca_stats")
head(exp_pc_stats)

barplot(t(t(exp_pc_stats[2,-1])))

# perform regression of metadata against the 
# RNAseq PCs vs. Covariates
Exp_meta_p <- NULL
Exp_meta_r_sq <- NULL
for (i in 1:21){
  sub_meta_p <- NULL
  sub_meta_r_sq <- NULL
  sub_meta <- sampleinfo_geno[,i]
  for (j in 1:20){
    #sub_meta = pca_scRNA_res_ov[,j]
    
    sub_lm <- lm(gene_exp_PC_t[,j]~sub_meta)
    
    sub_r_sq <- summary(sub_lm)$adj.r.squared
    sub_p <- anova(sub_lm)$Pr[1]
    
    sub_meta_p <- c(sub_meta_p,sub_p)
    sub_meta_r_sq <- c(sub_meta_r_sq,sub_r_sq)
  }
  Exp_meta_p <- rbind(Exp_meta_p,sub_meta_p)
  Exp_meta_r_sq <- rbind(Exp_meta_r_sq,sub_meta_r_sq)
}
rownames(Exp_meta_p) <- colnames(sampleinfo_geno)[1:21]
rownames(Exp_meta_r_sq) <- colnames(sampleinfo_geno)[1:21]
colnames(Exp_meta_p) <- colnames(pca1$rotation)[1:20]
colnames(Exp_meta_r_sq) <- colnames(pca1$rotation)[1:20]


 # plot the heatmap
Exp_meta_log_p <- -log10(Exp_meta_p)
saveRDS(Exp_meta_log_p, "Exp_meta_log_p.rds")


range(Exp_meta_log_p)
breaks <- unique(c(seq(0,10,0.025),seq(10,21,0.2)))
length(breaks)
hmcol2<-colorRampPalette(brewer.pal(9,"BuGn"))(455)

heatmap.2(Exp_meta_log_p, Rowv=FALSE, Colv=FALSE, dendrogram='none', trace='none',
          margins=c(10,10), col=hmcol2[1:which(breaks==round(max(Exp_meta_log_p)))-1],
          colsep=c(1:10), rowsep=c(1:20), sepwidth=c(0.05, 0.025),
          sub="Association between PCs of Expression and Covariates",
          symm=F,symkey=F,symbreaks=F,breaks=breaks[1:which(breaks==round(max(Exp_meta_log_p)))],
          cexRow = 0.8,cexCol = 0.8)

hmcol3<-colorRampPalette(brewer.pal(9,"Blues"))(5)
hmcol4<-hmcol3
hmcol4[1]<-"white"
max_num <- ifelse(max(Exp_meta_log_p)>4, max(Exp_meta_log_p), 5)
pdf("PCA_heatmap_gene_expression3.pdf", width=10, height = 7, family="ArialMT")
max_num <- ifelse(max(Exp_meta_log_p)>4, max(Exp_meta_log_p), 5)
heatmap.2(Exp_meta_log_p,Rowv=F,Colv=FALSE,dendrogram='none',trace='none',
          margins=c(10,10),colsep=c(1:20),rowsep=c(1:20),
          sepwidth=c(0.0025,0.0025),sepcolor="black", col=hmcol4, breaks=c(0,1.30103,2,3,4,max_num),
          key.xlab = "-log10 pvalue", sub="PCA heatmap of Gene Expression PC associations",
          cexRow = 0.9,cexCol = 0.9,symm=F,symkey=F,symbreaks=F, key=T)
dev.off()

gene_pca_stats <- read.table("../quants2/genes_90percent.pca_stats")
head(gene_pca_stats)
rownames(gene_pca_stats) <- gene_pca_stats$V1
gene_pca_stats <- gene_pca_stats[,-1]

PCs <- paste("PC",seq(1, 20,1), sep="_")
p.variance.explained10 <- unlist(c(gene_pca_stats[2,1:20]))
names(p.variance.explained10) <- paste("PC",seq(1, 20,1))
pdf("PCA_heatmap_gene_expression_CompVar.pdf", width=6, height = 4, family="ArialMT")
barplot(100*p.variance.explained10, las=2, xlab='', ylab='% Variance Explained', 
        main="Variance Explained by component", 
        ylim = c(0,100))
dev.off()

# Associate genotype PCs with covariates
sampleinfo_4geno <- sampleinfo_order[,-1]

geno_meta_p <- NULL
geno_meta_r_sq <- NULL
for (i in 2:12){
  sub_meta_p <- NULL
  sub_meta_r_sq <- NULL
  sub_meta <- sampleinfo_4geno[,i]
  for (j in 1:10){
    #sub_meta = pca_scRNA_res_ov[,j]
    
    sub_lm <- lm(geno_PC_t[,j]~sub_meta)
    
    sub_r_sq <- summary(sub_lm)$adj.r.squared
    sub_p <- anova(sub_lm)$Pr[1]
    
    sub_meta_p <- c(sub_meta_p,sub_p)
    sub_meta_r_sq <- c(sub_meta_r_sq,sub_r_sq)
  }
  geno_meta_p <- rbind(geno_meta_p,sub_meta_p)
  geno_meta_r_sq <- rbind(geno_meta_r_sq,sub_meta_r_sq)
}
rownames(geno_meta_p) <- colnames(sampleinfo_4geno)[2:12]
rownames(geno_meta_r_sq) <- colnames(sampleinfo_4geno)[2:12]
colnames(geno_meta_p) <- colnames(pca1$rotation)[1:10]
colnames(geno_meta_r_sq) <- colnames(pca1$rotation)[1:10]

geno_meta_log_p <- -log(geno_meta_p)

max_num2 <- ifelse(max(geno_meta_log_p)>4, max(geno_meta_log_p), 5)
pdf("PCA_heatmap_genotype2.pdf", width=10, height = 7, family="ArialMT")
heatmap.2(geno_meta_log_p,Rowv=F,Colv=FALSE,dendrogram='none',trace='none',
          margins=c(10,10),colsep=c(1:20),rowsep=c(1:20),
          sepwidth=c(0.025,0.025),sepcolor="black", col=hmcol4, breaks=c(0,1.30103,2,3,4,max_num2),
          key.xlab = "-log10 pvalue", sub="PCA Heatmap of Genotype PC associations",
          cexRow = 0.9,cexCol = 0.9,symm=F,symkey=F,symbreaks=F, key=T)
dev.off()

geno_pca_stats <- read.table("../../Genotype2/genotypes.pca.pca_stats")
head(geno_pca_stats)
rownames(geno_pca_stats) <- geno_pca_stats$V1
geno_pca_stats <- geno_pca_stats[,-1]

PCs <- paste("PC",seq(1, 10,1), sep="_")
p.variance.explained10_geno <- unlist(c(geno_pca_stats[2,1:10]))
names(p.variance.explained10_geno) <- paste("PC",seq(1, 10,1))
pdf("PCA_heatmap_genotype_CompVar.pdf", width=6, height = 4, family="ArialMT")
barplot(100*p.variance.explained10_geno, las=2, xlab='', ylab='% Variance Explained', 
        main="Variance Explained by component", 
        ylim = c(0,100))
dev.off()
