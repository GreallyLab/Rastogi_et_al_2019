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
library(ggthemes)
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
i<-1
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
# colors for graph
cols_rastogi <- c(rgb(4,51,255, maxColorValue=255), rgb(255,38,0, maxColorValue=255))

# create table to include combined violin plot
Cov_table_genotype_expr2 <- Cov_table_genotype_expr
Cov_table_genotype_expr2$geno <- "Combined"
Cov_table_genotype_expr_all <- rbind(Cov_table_genotype_expr, Cov_table_genotype_expr2)
Cov_table_genotype_expr_all$geno <- factor(Cov_table_genotype_expr_all$geno, 
                                           levels=c("Combined", "0/0", "0/1", "1/1"))
violin_plot <- ggplot(Cov_table_genotype_expr_all, aes(x=geno, y=expression, 
                                                       fill=factor(group), color=factor(group)))+
  geom_violin(color="black", alpha=0.6) +
  geom_point(position=position_jitterdodge()) +
  stat_summary(fun.y="mean", colour="black", geom="point", 
               shape =95, show.legend = FALSE, size = 15, 
               aes(group=factor(group)), alpha = .9, position=position_dodge(.9)) +
  scale_fill_manual(values =cols_rastogi, guide=FALSE) +
  scale_color_manual(values = cols_rastogi, labels=c("NwA", "OA"),name="Group") +
  xlab(paste("Genotype (", vcf$REF[idx_rsID], "->", vcf$ALT[idx_rsID], ")", sep = "")) +
  ylab("Expression (RPKM)") +
  ggtitle(paste(gene_name, "-", geneID, "-", rsID[i], " QTL plot", sep="")) +
  theme_tufte()
ggsave(file = paste(gene_name,"_", rsID[i], "_violin.pdf",sep = ""),
       plot=violin_plot, width=6, height=4, useDingbats=FALSE)
}
