Generating the cell subtype proportion gene expression table.

I decided to use the fpkm counts for the star mapped files because they are normalized within sample and the cell subtype proportion is only generated within sample anyways. Of note, the GTF used for these quantifications is GenCode which is different from the ensembl GTF used for the differential gene expression. It is highly doubtful that this will alter the genes used to determine cell subtype proportion.

First I generated the count columns for the files:

```bash
for f1 in *.exon.rpkm.bed;
do
SAMPLE="$(echo ${f1} | cut -d '.' -f1)"
cut -f7  $f1 > sample_exon_rpkm_$SAMPLE
done

for f1 in *.gene.rpkm.bed;
do
SAMPLE="$(echo ${f1} | cut -d '.' -f1)"
echo $SAMPLE
cut -f7  $f1 > sample_gene_rpkm_$SAMPLE
done
```

easy concat:

```bash

paste exon_rpkm_annotation.txt sample_exon_rpkm_* > CPDeepa_exon_rpkm.txt
paste gene_rpkm_annotation.txt sample_gene_rpkm_* > CPDeepa_gene_rpkm.txt

 # Sanity check - answer should correspond to the number of samples +6
awk '{print NF}' CPDeepa_exon_rpkm.txt | sort -nu | tail -n 1 # 109 = 102+6 +1
awk '{print NF}' CPDeepa_gene_rpkm.txt | sort -nu | tail -n 1 # 109 = 102+6 +1
```

Now use R to clean up the data `Deepa_cellProp_table.R`:

```R
 # To create an expression table to be used to generate cell subtype proportions
 # set the directory
setwd("/Volumes/home/greally-lab/Deepa_Andrew/eQTL/Bams/quants2/")

 # filtering the genes that are not well expressed in majority of samples
options(stringsAsFactors = F)
options(scipen = 999)
 # load libraries
library(data.table)
library(genefilter)
 # read in quant matrix
mat <- fread("CPDeepa_gene_rpkm.txt", header=T, sep="\t")
mat <- as.data.frame(mat)
dim(mat) # 27518    109
colnames(mat)
 # remove gene information and B26 bam file
mat <- mat [,-c(1:3,5:6,83)]
 # switch A25 bam to B26
colnames(mat)[colnames(mat)=="A25"] <- "B26"
head(mat)

 # remove sub ENSG ID numbers
mat[,1] <- do.call(rbind, strsplit(mat[,1], "[.]"))[,1]

 # write out the table
write.table(mat, "CPDeepa_matrix_gene_rpkm_ready.txt", col.names = T, row.names = F,
            append = F, quote = F, sep = "\t")
```
