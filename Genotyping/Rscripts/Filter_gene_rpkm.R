setwd("/Volumes/home/greally-lab/Deepa_Andrew/eQTL/Bams/quants/")

# filtering the genes that are not well expressed in majority of samples 
options(stringsAsFactors = F)
options(scipen = 999)
# load libraries
library(data.table)
mat <- read.table("matrix_gene_rpkm.txt", header=T, sep="\t")
head(mat,1)
mat2 <- NULL
for (i in 1:nrow(mat)){
  if (mean(mat[i,-c(1:6)]==0)<.5){
    mat2 <- rbind(mat2, mat[i,])
  }
}
dim(mat2) 
# 16607    82
head(mat2)

# these are the genes where > 50% have nonzero rpkms.
# I chose this because Geavardis paper did this :
"In general, we used only elements quantified in >50% of individuals 
unless mentioned otherwise (Table S3)."

library(genefilter)
hist(rowMeans(log(mat2[,-c(1:6)]+2)))
means_row <- rowMeans(mat2[,-c(1:6)])
summary(means_row)

# Table S3 says the following:
# "All eQTL counts are for autosomal genes, with a filter of 
# quantification in >90% of samples for genes, exons, transcripts"
mat2 <- mat
table(mat2$chr)
mat2_auto <- mat2[-which(mat2$chr == "M" | mat2$chr == "X" | mat2$chr == "Y"),]
table(mat2_auto$chr)
dim(mat2_auto)
# 16047    82

mat3 <- NULL
for (i in 1:nrow(mat2_auto)){
  if (mean(mat2_auto[i,-c(1:6)]==0)<.1){
    mat3 <- rbind(mat3, mat2_auto[i,])
  }
}
dim(mat3)
# 13210    82
head(mat3)
# this number of genes is similar to the 13,703 genes used in the Geavardis 
# paper for eQTL analysis (see Table S3).

# Now I'll make the sample swaps with the IDs
# geno orig
# B23	A24
# B26	A25
# A23	A21
# A24	A22
# A25	A23
colnames(mat3) <- colnames(mat2_auto)

colnames(mat3)[colnames(mat3)=="A24"] <- "B23"
colnames(mat3)[colnames(mat3)=="A25"] <- "B26"
colnames(mat3)[colnames(mat3)=="A23"] <- "A25"
colnames(mat3)[colnames(mat3)=="A21"] <- "A23"
colnames(mat3)[colnames(mat3)=="A22"] <- "A24"

# Now I need to make sure that I am calling my samples the same name and get 
# them in the same order
geno_pca <- read.table("../../Genotype/genotypes.pca.pca", header=T)
head(geno_pca)
geno_IDs <- substring(colnames(geno_pca)[-1], 2, 7)

idx_match <- match(colnames(mat3[,-c(1:6)]), do.call(rbind, strsplit(geno_IDs, "_"))[,2])
same_ids <- which(colnames(mat3) %in% do.call(rbind, strsplit(geno_IDs, "_"))[,2])
length(same_ids)
colnames(mat3)[-same_ids]

idx_match2 <- match(do.call(rbind, strsplit(geno_IDs, "_"))[,2],colnames(mat3[,-c(1:6)]))


mat4 <- mat3
colnames(mat4)[-c(1:6)] <- geno_IDs[idx_match]
order(colnames(mat4)[-c(1:6)])
mat4_lite <- mat4[,-c(1:6)]
mat4_lite <- mat4_lite[,idx_match2]
head(mat4_lite)

mat5 <- cbind(mat4[,c(1:6)], mat4_lite)
dim(mat5)
head(mat5)

write.table(mat5, "Filtered_matrix_gene_rpkm.txt", col.names = T, row.names = F,
            append = F, quote = F, sep = "\t")
