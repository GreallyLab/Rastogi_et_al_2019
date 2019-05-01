### creating the covariate tables for the various runs of QTL analysis

# set options
options(scipen = 999, stringsAsFactors = F)

# load in the gene expression PCs
gene_exp_PC <- read.table("../Bams/PCA/genes_90percent.pca", header=T, 
                          colClasses = c("character", rep("numeric",76)))
colnames(gene_exp_PC)[-1] <- substr(colnames(gene_exp_PC)[-1], 2,7)
head(gene_exp_PC)

# load in the genotyping PCs
geno_PC <- read.table("../Genotype/genotypes.pca.pca", header = T)
head(geno_PC)
colnames(geno_PC)[-1] <- substr(colnames(geno_PC)[-1], 2,7)
str(geno_PC)

# load in the sample covariate data
sampleinfo <- read.table("../clinicaltable_for_eqtl.txt", header=T, stringsAsFactors = F)
head(sampleinfo)
str(sampleinfo)
dim(sampleinfo)

sampleinfo$Sample
# making qualitative data factors
sampleinfo[,2] <- as.factor(sampleinfo[,2]) #group
# sampleinfo[,3] <- as.factor(sampleinfo[,3]) #age
sampleinfo[,4] <- as.factor(sampleinfo[,4]) #sex
sampleinfo[,5] <- as.factor(sampleinfo[,5]) #ethnicity
# sampleinfo[,6] <- as.factor(sampleinfo[,6]) #insulin 
sampleinfo[,12] <- as.factor(sampleinfo[,12]) # seq batch

# I need to switch up the columns for the samples that were swapped
# this was already done for the gene_exp data before pca analysis
sampleinfo <- sampleinfo[-which(sampleinfo$Sample=="B23"),]
sampleinfo <- sampleinfo[-which(sampleinfo$Sample=="B26"),]
dim(sampleinfo)
sampleinfo$Sample
sampleinfo$Sample[sampleinfo$Sample=="A24"] <- "B23"
sampleinfo$Sample[sampleinfo$Sample=="A25"] <- "B26"
sampleinfo$Sample[sampleinfo$Sample=="A23"] <- "A25"
sampleinfo$Sample[sampleinfo$Sample=="A21"] <- "A23"
sampleinfo$Sample[sampleinfo$Sample=="A22"] <- "A24"

# ordering the samples to sync with the gene expr data as well as remove uneeded rows 
idx_match_meta_pca <- na.omit(match(do.call(rbind, strsplit(colnames(geno_PC)[-1], "_"))[,2],
                                    sampleinfo$Sample))
sampleinfo_order <- sampleinfo[idx_match_meta_pca,]
length(do.call(rbind, strsplit(colnames(geno_PC)[-1], "_"))[,2])
length(idx_match_meta_pca)

head(sampleinfo_order)
dim(sampleinfo_order)
tail(sampleinfo_order)
colnames(sampleinfo_order)
sampleinfo_order_t <- as.data.frame(t(sampleinfo_order[,-1]))
head(sampleinfo_order_t)
str(sampleinfo_order_t)
colnames(sampleinfo_order_t) <- sampleinfo_order[,1]
sampleinfo_order_t_fil <- sampleinfo_order_t[c(8,11),]
head(sampleinfo_order_t_fil)
sampleinfo_order_t_fil <- cbind(rownames(sampleinfo_order_t_fil), sampleinfo_order_t_fil)
head(sampleinfo_order_t_fil)
colnames(sampleinfo_order_t_fil) <- colnames(geno_PC)
sampleinfo_order_t_fil[2,] <- paste0("Batch_",sampleinfo_order_t_fil[2,])

cov_1235_batch_pctCode_12 <- rbind(gene_exp_PC[c(1:3, 5),], geno_PC[1:2,], 
                                   sampleinfo_order_t_fil)
head(cov_1235_batch_pctCode_12)

write.table(cov_1235_batch_pctCode_12, 
            "no_normal_123n5_batch_prtCode_12/Cov_1235_batch_pctCode_12.txt",
            row.names = F, col.names = T, sep = "\t", append = F, quote = F)


# without the clin. covariates
cov_1235_12 <- rbind(gene_exp_PC[c(1:3, 5),], geno_PC[1:2,])
head(cov_1235_12)

write.table(cov_1235_12, 
            "no_normal_123n5_12/Cov_123n5_12.txt",
            row.names = F, col.names = T, sep = "\t", append = F, quote = F)

# without the 5th PC in expression
cov_123_batch_pctCode_12 <- rbind(gene_exp_PC[c(1:3),], geno_PC[1:2,], 
                                   sampleinfo_order_t_fil)
head(cov_123_batch_pctCode_12)

write.table(cov_123_batch_pctCode_12, 
            "no_normal_123_batch_prtCode_12/Cov_123_batch_pctCode_12.txt",
            row.names = F, col.names = T, sep = "\t", append = F, quote = F)

# with the 4th PC in expression
cov_12345_12 <- rbind(gene_exp_PC[c(1:5),], geno_PC[1:2,])
head(cov_12345_12)

write.table(cov_12345_12, 
            "no_normal_12345_12/Cov_12345_12.txt",
            row.names = F, col.names = T, sep = "\t", append = F, quote = F)

# with the 9th PC in expression
cov_123459_12 <- rbind(gene_exp_PC[c(1:5,9),], geno_PC[1:2,])
head(cov_123459_12)

write.table(cov_123459_12, 
            "no_normal_123459_12/Cov_123459_12.txt",
            row.names = F, col.names = T, sep = "\t", append = F, quote = F)

# with the 9th PC in expression and 3rd geno PC 
cov_123459_123 <- rbind(gene_exp_PC[c(1:5,9),], geno_PC[1:3,])
head(cov_123459_123)

write.table(cov_123459_123, 
            "no_normal_123459_123/Cov_123459_123.txt",
            row.names = F, col.names = T, sep = "\t", append = F, quote = F)

# no expression PCs
cov_batch_pctCode_12 <- rbind(geno_PC[1:2,], sampleinfo_order_t_fil)
head(cov_batch_pctCode_12)

write.table(cov_batch_pctCode_12, 
            "no_normal_batch_prtCode_12/Cov_batch_pctCode_12.txt",
            row.names = F, col.names = T, sep = "\t", append = F, quote = F)

# with the 9th PC in expression and 3rd geno PC 
cov_12359_12 <- rbind(gene_exp_PC[c(1:3,5,9),], geno_PC[1:2,])
head(cov_12359_12)

write.table(cov_12359_12, 
            "no_normal_12359_123/Cov_12359_12.txt",
            row.names = F, col.names = T, sep = "\t", append = F, quote = F)


# with the 9th PC in expression and 3rd geno PC 
Cov_12459_12 <- rbind(gene_exp_PC[c(1:3,4:5,9),], geno_PC[1:2,])
head(Cov_12459_12)

write.table(Cov_12459_12, 
            "no_normal_12459_12/Cov_12459_12.txt",
            row.names = F, col.names = T, sep = "\t", append = F, quote = F)
