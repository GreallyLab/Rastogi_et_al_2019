args <- commandArgs(trailingOnly = T)
print(args)
gene_name <- as.character(args[1])
geneID <- as.character(args[2])
eqtlFile <- as.character(args[3])
vcfFile <- as.character(args[4])
covFile <- as.character(args[5])
expFile <- as.character(args[6])

# vcfFile <- "SplX_converted_noHH_eQTL_snpfilt_hg38_sort.vcf"

# load in libraries and set options
library(ggplot2)
library(data.table)
library(scales)
options(stringsAsFactors = F)

# read in the expression data
cat(paste("Reading in expression table: ", expFile,".\n",sep=""))
exprs_dat <- fread(expFile, header=T)
exprs_dat <- as.data.frame(exprs_dat)
exprs_ID <- grep(geneID, exprs_dat$gene)

# read in the eQTL data
cat(paste("Reading in selected eQTL table: ", eqtlFile,".\n",sep=""))
eQTL_dat <- fread(eqtlFile, header=F)
eQTL_dat <- as.data.frame(eQTL_dat)
rsID <- eQTL_dat$V8[grep(geneID, eQTL_dat$V1)]

if(length(rsID)> 1){
  cat(paste("There is more than 1 variant associated with changes in ", gene_name,
            "'s gene expression", sep=""))
}

# read in vcf file
cat(paste("Reading in vcf: ", vcfFile,".\n",sep=""))
vcf <- fread(vcfFile, header=T, skip=27)
colnames(vcf)[1] <- "CHROM"
vcf <- as.data.frame(vcf)

# read in clinical covariates (ethnicity and Normal vs. Obese)
cat(paste("Reading in Covariate Table: ", covFile,".\n",sep=""))
Cov_table <- read.table(covFile, header=T)

# order the VCF information to match the covariate table
sort_vcf <- match(Cov_table$Sample, do.call(rbind, strsplit(colnames(vcf)[10:ncol(vcf)], "_"))[,2])+9

for (i in 1:length(rsID)){

idx_rsID <- which(vcf$ID == rsID[i])

vcf_ID <- t(vcf[idx_rsID, sort_vcf])

# combine the vcf file information with the covariate table
Cov_table_genotype <- cbind(Cov_table,vcf_ID)

# do the expression data columns match the Covariate table order?
if (any(!(Cov_table_genotype$Sample==do.call(rbind,
                                             strsplit(colnames(exprs_dat)[7:ncol(exprs_dat)],
                                                      "_"))[,2]))) {
  sort_exp <- match(Cov_table_genotype$Sample, do.call(rbind, strsplit(colnames(exprs_dat)[7:ncol(exprs_dat)], "_"))[,2])+6
  exprs_dat_gene <- t(exprs_dat[exprs_ID, sort_exp])
} else {
  exprs_dat_gene <- t(exprs_dat[exprs_ID, 7:ncol(exprs_dat)])
}

# combine expression with geno and covariate data
Cov_table_genotype_expr <- cbind(Cov_table_genotype,exprs_dat_gene)

# fix the column names
colnames(Cov_table_genotype_expr)[(ncol(Cov_table)+1):(ncol(Cov_table)+2)] <- c("geno", "expression")

cat(paste("Generating Plots",".\n",sep=""))
pdf(file = paste(geneID,"_", rsID[i],"_QTL_by_Ethnicity", ".pdf",sep = ""), width=6,
    height=4.5,family="ArialMT")
print(ggplot(Cov_table_genotype_expr, aes(x=geno, y=expression))+
  geom_boxplot(fill='#A4A4A4', color="black") +
  geom_point(position=position_jitterdodge(),aes(group=geno, color=factor(ethnicity)))+
  scale_color_discrete(labels=c("AA", "Hispanic"), name="Ethnicity")+
  xlab(paste("Genotype (", vcf$REF[idx_rsID], "->", vcf$ALT[idx_rsID], ")", sep = "")) +
  ylab("Expression (RPKM)") +
  ggtitle(paste(gene_name, "-", geneID, "-", rsID[i], " QTL plot", sep="")))
dev.off()

pdf(file = paste(geneID,"_", rsID[i],"_QTL_by_batch", ".pdf",sep = ""), width=6,
    height=4.5, family="ArialMT")
print(ggplot(Cov_table_genotype_expr, aes(x=geno, y=expression))+
  geom_boxplot(fill='#A4A4A4', color="black") +
  geom_point(position=position_jitterdodge(jitter.width=.75),aes(group=geno, color=factor(rnaseq_batch)))+
  scale_color_discrete(name="RNAseq batch") +
  xlab(paste("Genotype (", vcf$REF[idx_rsID], "->", vcf$ALT[idx_rsID], ")", sep = "")) +
  ylab("Expression (RPKM)") +
  ggtitle(paste(gene_name, "-", geneID, "-", rsID[i], " QTL plot", sep="")))
dev.off()

pdf(file = paste(geneID,"_", rsID[i],"_QTL_by_group", ".pdf",sep = ""), width=6,
    height=4.5, family="ArialMT")
print(ggplot(Cov_table_genotype_expr, aes(x=geno, y=expression))+
  geom_boxplot(fill='#A4A4A4', color="black") +
  geom_point(position=position_jitterdodge(),aes(group=geno, color=factor(group)))+
  scale_color_discrete(labels=c("NormalWT", "Obese"),name="Group") +
  xlab(paste("Genotype (", vcf$REF[idx_rsID], "->", vcf$ALT[idx_rsID], ")", sep = "")) +
  ylab("Expression (RPKM)") +
  ggtitle(paste(gene_name, "-", geneID, "-", rsID[i], " QTL plot", sep="")))
dev.off()

pdf(file = paste(geneID,"_", rsID[i],"_QTL_by_homa", ".pdf",sep = ""), width=6,
    height=4.5, family="ArialMT")
midp_homa <- range(Cov_table_genotype_expr$homa)[2]-((range(Cov_table_genotype_expr$homa)[2] - range(Cov_table_genotype_expr$homa)[1])/2)
print(ggplot(Cov_table_genotype_expr, aes(x=geno, y=expression))+
  geom_boxplot(fill='#A4A4A4', color="black") +
  geom_point(position=position_jitterdodge(),aes(group=geno, color=homa))+
  scale_colour_gradient2(low = muted("red"), mid = "white",
                         high = muted("blue"), midpoint = midp_homa, space = "Lab",
                         na.value = "grey50", guide = "colourbar", name="HOMA") +
  xlab(paste("Genotype (", vcf$REF[idx_rsID], "->", vcf$ALT[idx_rsID], ")", sep = "")) +
  ylab("Expression (RPKM)") +
  ggtitle(paste(gene_name, "-", geneID, "-", rsID[i], " QTL plot", sep="")))
dev.off()

pdf(file = paste(geneID,"_", rsID[i],"_QTL_by_homa_colorMedian", ".pdf",sep = ""), width=6,
    height=4.5, family="ArialMT")
print(ggplot(Cov_table_genotype_expr, aes(x=geno, y=expression))+
  geom_boxplot(fill='#A4A4A4', color="black") +
  geom_point(position=position_jitterdodge(),aes(group=geno, color=homa))+
  scale_colour_gradient2(low = muted("red"), mid = "white",
                         high = muted("blue"), midpoint = median(Cov_table_genotype_expr$homa), space = "Lab",
                         na.value = "grey50", guide = "colourbar", name="HOMA") +
  xlab(paste("Genotype (", vcf$REF[idx_rsID], "->", vcf$ALT[idx_rsID], ")", sep = "")) +
  ylab("Expression (RPKM)") +
  ggtitle(paste(gene_name, "-", geneID, "-", rsID[i], " QTL plot", sep="")))
dev.off()
}
