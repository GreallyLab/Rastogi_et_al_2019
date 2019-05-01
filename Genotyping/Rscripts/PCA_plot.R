setwd("/Volumes/home/greally-lab/Deepa_Andrew/eQTL/Genotype/")

options(stringsAsFactors = F, scipen = 999)
library(ggplot2)

pca <- read.table("genotypes.pca.pca", header=T)
rownames(pca) <- pca$SampleID
dim(pca) # 76   77
head(pca)
pca <- pca[,-1]
dim(pca) # 76   76
pca_t <- t(pca)
head(pca_t)
colnames(pca_t) <- do.call(rbind,strsplit(colnames(pca_t), split="_"))[,5]
pca_t <- as.data.frame(pca_t)

pca_t$group <- "normal wt"
pca_t$group[grep("B",rownames(pca_t))] <- "obese/ow wt"

metadata <- read.table("../clinicaltable_for_eqtl.txt", header=T)
head(metadata)

pca_t_samples <-  data.frame(ID=do.call(rbind,strsplit(rownames(pca_t), split="_"))[,2])

idx_samples <- match(pca_t_samples$ID, metadata$Sample)
metadata_ordered <- metadata[idx_samples,]

idx_female <- which(metadata_ordered$sex==1)
?match
pca_t$sex <- "Male"
pca_t$sex[idx_female] <- "Female"

pdf(file = "PCA_1_2.pdf", height = 4.75, width=7)
p_12 <- ggplot(pca_t, aes(PC1, PC2)) +
  geom_point(aes(colour=factor(group), shape=factor(sex)), alpha = .7) +
  scale_color_discrete(name="Group") +
  scale_shape_discrete(name="Sex") +
  ylab("PC2")+
  xlab("PC1")+
  theme(axis.text.x = element_blank(),axis.text.y = element_blank(),
        axis.ticks = element_blank())+
  theme_classic()+
  ggtitle("Genotype Data PC1 and PC2")
p_12
dev.off()

pdf(file = "PCA_1_3.pdf", height = 4.75, width=7)
p_13 <- ggplot(pca_t, aes(PC1, PC3)) +
  geom_point(aes(colour=factor(group), shape=factor(sex)), alpha = .7) +
  scale_color_discrete(name="Group") +
  scale_shape_discrete(name="Sex") +
  ylab("PC3")+
  xlab("PC1")+
  theme(axis.text.x = element_blank(),axis.text.y = element_blank(),
        axis.ticks = element_blank())+
  theme_classic()+
  ggtitle("Genotype Data PC1 and PC3")
p_13
dev.off()

pdf(file = "PCA_3_2.pdf", height = 4.75, width=7)
p_32 <- ggplot(pca_t, aes(PC3, PC2)) +
  geom_point(aes(colour=factor(group), shape=factor(sex)), alpha = .7) +
  scale_color_discrete(name="Group") +
  scale_shape_discrete(name="Sex") +
  ylab("PC2")+
  xlab("PC3")+
  theme(axis.text.x = element_blank(),axis.text.y = element_blank(),
        axis.ticks = element_blank())+
  theme_classic()+
  ggtitle("Genotype Data PC3 and PC2")
p_32
dev.off()

pdf(file = "PCA_1_4.pdf", height = 4.75, width=7)
p_14 <- ggplot(pca_t, aes(PC1, PC4)) +
  geom_point(aes(colour=factor(group), shape=factor(sex)), alpha = .7) +
  scale_color_discrete(name="Group") +
  scale_shape_discrete(name="Sex") +
  ylab("PC4")+
  xlab("PC1")+
  theme(axis.text.x = element_blank(),axis.text.y = element_blank(),
        axis.ticks = element_blank())+
  theme_classic()+
  ggtitle("Genotype Data PC1 and PC4")
p_14
dev.off()

pdf(file = "PCA_1_5.pdf", height = 4.75, width=7)
p_15 <- ggplot(pca_t, aes(PC1, PC5)) +
  geom_point(aes(colour=factor(group), shape=factor(sex)), alpha = .7) +
  scale_color_discrete(name="Group") +
  scale_shape_discrete(name="Sex") +
  ylab("PC4")+
  xlab("PC1")+
  theme(axis.text.x = element_blank(),axis.text.y = element_blank(),
        axis.ticks = element_blank())+
  theme_classic()+
  ggtitle("Genotype Data PC1 and PC5")
p_15
dev.off()

pdf(file = "PCA_1_6.pdf", height = 4.75, width=7)
p_16 <- ggplot(pca_t, aes(PC1, PC6)) +
  geom_point(aes(colour=factor(group), shape=factor(sex)), alpha = .7) +
  scale_color_discrete(name="Group") +
  scale_shape_discrete(name="Sex") +
  ylab("PC4")+
  xlab("PC1")+
  theme(axis.text.x = element_blank(),axis.text.y = element_blank(),
        axis.ticks = element_blank())+
  theme_classic()+
  ggtitle("Genotype Data PC1 and PC6")
p_16
dev.off()

## outputting the first 3 PCs to add to covariate table
head(pca,1)
rownames(pca) <- paste("Genotype_", 
                       do.call(rbind,strsplit(rownames(pca), split="_"))[,5],
                       sep="")
colnames(pca) <- substring(colnames(pca), 2, 7)
cov3 <- pca[1:3,]
head(cov3)
write.table(cov3, "Geno_cov_3PCs.txt", col.names = T, row.names = T,
            append = F, quote = F, sep = "\t")
