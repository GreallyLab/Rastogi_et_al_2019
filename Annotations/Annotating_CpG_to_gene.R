# load the libraries
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(GenomicFeatures)
library(RMariaDB)

# set options
options(scipen=999, stringsAsFactors = FALSE)

enh_CpG <- readRDS("../enh_cpg.rds")

head(enh_CpG[c(1:10)])
dim(enh_CpG) # 57939    30

# unfortunately, these annotations apply to the enhancers, so not as high resolution
# as annotating the CpGs themselves to specific genes
dens <- density(enh_CpG$distanceToTSS)
plot(dens)

plot(dens, xlim= c(-10000,10000), main = "Density of enhancer CpGs around TSS")
plot(dens, xlim= c(-5000,5000), main = "Density of enhancer CpGs around TSS")
sum(enh_CpG$distanceToTSS==0) # 15714, overlap the TSS
sum(enh_CpG$distanceToTSS>0) #
sum(enh_CpG$distanceToTSS<0)
head(enh_CpG[enh_CpG$distanceToTSS==0,1:10])

# excluding gene body CpGs
enh_CpG_GB <- enh_CpG[!is.na(enh_CpG$GeneBody),]
nrow(enh_CpG_GB) # 36220
enh_CpG_noGB <- enh_CpG[is.na(enh_CpG$GeneBody),]
nrow(enh_CpG_noGB) # 21719


 ### CREATING CpG Annotation
# create CpG bed file:
enh_CpG_bed <- enh_CpG[,c(3,4,4,1,2)]
enh_CpG_bed$pos.x.x.1 <- enh_CpG_bed$pos.x.x.1+1
enh_CpG_bed$strand <- "."
head(enh_CpG_bed)
# grep("198279", enh_CpG_bed$X.tid)
write.table(enh_CpG_bed, "enh_CpG_bed.txt", row.names = F, col.names = TRUE,
            sep="\t", append = F, quote = F)
cpg_bed_file <- "Randoms/enh_CpG_bed.txt"

# read in the bed as a peak file for ChIPSeeker
CpG_enh <- readPeakFile(cpg_bed_file)
CpG_enh

# grab the refseq genes
hg38.refseq.db <- makeTxDbFromUCSC(genome="hg38", table="refGene")
 # UCSC data anomaly in 581 transcript(s): the cds cumulative length is not a multiple of 3 for transcripts
 # ignore

# annotate the CpGs
enhAnno <- annotatePeak(CpG_enh, tssRegion=c(-1000, 1000),
                        TxDb=hg38.refseq.db,
                        annoDb="org.Hs.eg.db",overlap='all', ignoreOverlap=TRUE,
                        sameStrand = FALSE)
saveRDS(enhAnno, "CpG_in_enhAnno_obj_2.rds")
enhAnno <- readRDS("../Deepa-helptagging/CpG_Annotations/Randoms/CpG_in_enhAnno_obj.R")
enhAnno_df <- as.data.frame(enhAnno)
head(enhAnno_df)
summary(enhAnno_df$distanceToTSS)
sum(enhAnno_df$distanceToTSS==0) #0, overlap the TSS
sum(enhAnno_df$distanceToTSS>0) #34918, downstream of the TSS
sum(enhAnno_df$distanceToTSS<0) #23021, upstream the TSS
sum(enhAnno_df$distanceToTSS<=50000 & enhAnno_df$distanceToTSS>=-50000) # 48426
sum(is.na(enhAnno_df$distanceToTSS)) # 0
range(enhAnno_df$distanceToTSS)

dens_anno <- density(enhAnno_df$distanceToTSS)

pdf(file = "Density_plots_TSS_CpG_enh.pdf", width = 7, height = 5, family="ArialMT")
plot(dens_anno, main="CpG density around TSSs")
plot(dens_anno, main="CpG density around TSSs\n+/- 10Kb zoom", xlim=c(-10000,10000))
plot(dens_anno, main="CpG density around TSSs\n+/- 50Kb zoom", xlim=c(-50000,50000))
dev.off()

pdf(file = "CpG_enh_annotation_proportions.pdf", width = 7, height = 5, family="ArialMT")
plotAnnoPie(enhAnno, main="\n\nProportion of Enhancer CpG Annotations")
plotAnnoBar(enhAnno)
upsetplot(enhAnno)
dev.off()

# add these specific enhancer annotations to CpG Annotation
library(data.table)
CpG_Anno <- fread("CpGID_Prom_GB_Enh2_annot.txt")
CpG_Anno
CpG_Anno <- as.data.frame(CpG_Anno)
head(CpG_Anno)
 # only add the gene information if enhancer is within 50kb
enhAnno_df_50kb <- enhAnno_df[enhAnno_df$distanceToTSS<=50000 &
                                enhAnno_df$distanceToTSS>=(-50000),]
dim(enhAnno_df_50kb) # 48426
head(enhAnno_df_50kb)
 # add some additional information
enhAnno_df_50kb_min <- enhAnno_df_50kb[,c(6,7,9,17,18,19)]
head(enhAnno_df_50kb_min)
CpG_Anno_wEnh <- merge(x= CpG_Anno, y= enhAnno_df_50kb_min, by="X.tid", all.x=TRUE)
head(CpG_Anno_wEnh)

sum(CpG_Anno_wEnh$Enhancer.x==CpG_Anno_wEnh$Enhancer.y, na.rm = T) # 48,426
# don't need the second enhancer name column
CpG_Anno_wEnh <- CpG_Anno_wEnh[,-7]
colnames(CpG_Anno_wEnh)[6] <- "Enhancer"
dim(CpG_Anno_wEnh)
head(CpG_Anno_wEnh)

test <- fread("Enhancer_annotations.txt", header=T)
head(test)
summary(test$distanceToTSS)

write.table(CpG_Anno_wEnh, "CpG_Anno_wEnh50kb_info.txt", append = F, quote = F,
                        col.names = T, row.names = F, sep = "\t", na = "NA")
