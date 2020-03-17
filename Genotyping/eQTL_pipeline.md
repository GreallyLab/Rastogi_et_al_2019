# eQTL Anaylsis
## Andrew D. Johnston
## 02/27/18

### eQTL preparation and performing eQTL analysis

####Table of Contents

1. [Preparing Bam files](#prepbam)
2. [QC of .bam files](#bamqc)
3. [Preparing covariate file](#Covs)
4. [Preparaing phenotype (expression) file](#)
5. [Preparing gtf](#gtf)
6. [Preparing the genotype data](#genoprep)
7. [Explore PCA covariates](#PCA)
8. [](#)
9.


<a name="prepbam">

1. Preparing the BAM files

I moved all of the bam files previously mapped using STAR onto hg38 with a ensembl 83 release GTF to `/home/greally-lab/Deepa_Andrew/eQTL/Bams/raw_bams`. The reads were originally mapped using this command:

```bash
for f1 in *_trimmed.fq.gz   
do
SAMPLE="$(echo ${f1} | cut -d '_' -f5 | cut -d '.' -f1)"
echo ${SAMPLE}
echo ${f1}
qsub -S /bin/bash -N ${SAMPLE}_alignEns -cwd -l h_vmem=5.6G -j y -pe smp 20 << EOF
module load samtools
module load STAR
STAR --runThreadN 20 --genomeDir /home/greally-lab/indexes/hg38/Star --readFilesIn ${f1} --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --chimSegmentMin 20 --quantMode TranscriptomeSAM GeneCounts --quantTranscriptomeBan Singleend --sjdbGTFfile /home/greally-lab/indexes/hg38/ensembl/Homo_sapiens.GRCh38.83.gtf --sjdbGTFchrPrefix chr --sjdbOverhang 99 --outFileNamePrefix Mapped_Ensembl/${SAMPLE}
EOF
done
```

Samples A01 and B06 were ran twice (denoted by A01a and B06a). I chose the bam file with more uniquely mapped reads - A01a and B06. The unused bams were placed in `/home/greally-lab/Deepa_Andrew/eQTL/Bams/unsused_bams`. Samples A49 and A50 were also ran twice but both runs were of similar quality (denoted by A49.1 and A50.1). I chose to merge the filtered bam files (see below).

Use the following to filter for only uniquely aligned reads in the bam files.
```bash
for f1 in *.bam;
do
SAMPLE="$(echo ${f1} | cut -c1-3)"
qsub -S /bin/bash -N filter_${SAMPLE} -cwd -q highmem.q -l h_vmem=10G -j y << EOF
module load samtools
samtools view ${f1} | awk '\$5=="255"' > ${f1}tmp.sam
samtools view -H ${f1} | cat - ${f1}tmp.sam | samtools view -Sb - > ../filt_bams/${f1%.bam}_filt.bam
rm ${f1}tmp.sam
EOF
done
```

Merging the A49 and A50 files
```bash
qsub -S /bin/bash -N merge_49 -cwd -q highmem.q -l h_vmem=10G -j y << EOF
module load samtools
samtools merge A49_merged_Aligned.sortedByCoord.out_filt.bam A49.1Aligned.sortedByCoord.out_filt.bam A49Aligned.sortedByCoord.out_filt.bam
EOF
qsub -S /bin/bash -N merge_50 -cwd -q highmem.q -l h_vmem=10G -j y << EOF
module load samtools
samtools merge A50_merged_Aligned.sortedByCoord.out_filt.bam A50.1Aligned.sortedByCoord.out_filt.bam A50Aligned.sortedByCoord.out_filt.bam
EOF
```
I placed the original files in `merged_bams` and then changed the merged file names to the original e.g. A49_merged to A49

Then I indexed all of the bam files

```bash
for f1 in *.bam;
do
SAMPLE="$(echo ${f1} | cut -c1-3)"
qsub -S /bin/bash -N index_${SAMPLE} -cwd -q highmem.q -l h_vmem=5G -j y << EOF
module load samtools
samtools index $f1
EOF
done
```

<a name="bamqc">

I also wanted to assess the transcript integrity of each sample as a QC measure. I used three measures (1) the proportion of protein coding reads (out of total reads), (2) RSeQC's TIN (trasncript integrity) score, and (3) mRIN score, an approach similar to TIN. I had already calculated the proportion of reads counted toward protein coding genes, using the ensembl GTF file `Homo_sapiens.GRCh38.83.gtf`.


To calculate both TIN and mRIN, I needed to generate a bed12 file containing only the longest transcripts for each gene though. This is the GTF used for gene expression quantification `/home/greally-lab/indexes/hg38/GenCode/gencode.v24.primary_assembly.annotation.lincNprot.gtf.gz`. I used this code to obtain the longest transcript.

```bash
 # I want protein_coding and lincRNA
egrep '(protein_coding|lincRNA)' gencode.v24.primary_assembly.annotation.gtf > gencode.v24.primary_assembly.annotation.lincNprot.gtf

 #sanity check
awk '$3=="gene" {print $12}' gencode.v24.primary_assembly.annotation.gtf | tr -d '";' | sort | uniq -c
awk '$3=="gene"' gencode.v24.primary_assembly.annotation.lincNprot.gtf | wc -l # 27,518 is correct
awk '$3=="transcript" {print $10}' gencode.v24.primary_assembly.annotation.lincNprot.gtf | tr -d '";' | sort | uniq | wc -l # 27,748 is not correct
 # turns out that some processed transcript genes have transcripts that are considered lincRNAs
 # looking at the composition of transcript types:
awk '$3=="transcript" {print $14}' gencode.v24.primary_assembly.annotation.lincNprot.gtf | tr -d '";' | sort | uniq -c
 # 12656 lincRNA
 #    48 polymorphic_pseudogene
 #   867 processed_transcript
 # 143984 protein_coding

 # I tried to use cgat to extract only protein coding and lincRNA genes
 # (and their transcripts) but this also extracted the trasncripts of processed transcripts and pseudogenes
 # observe:
source /gs/gsfs0/users/anjohnst/cgat-install/conda-install/bin/activate cgat-s
cgat gtf2gtf --method=filter --filter-method=proteincoding -I gencode.v24.primary_assembly.annotation.gtf> gencode.v24.primary_assembly.annotation.prot.gtf
cgat gtf2gtf --method=filter --filter-method=lincrna -I gencode.v24.primary_assembly.annotation.gtf> gencode.v24.primary_assembly.annotation.lincrna.gtf
 # then when one checks the lincrna gtf
awk '$3=="transcript" {print $10}' gencode.v24.primary_assembly.annotation.lincrna.gtf | sort | uniq | wc -l # 7868
awk '$3=="gene" {print $10}' gencode.v24.primary_assembly.annotation.lincrna.gtf | sort | uniq | wc -l # 7674

 # therefore I just grep'd to remove any lines with polymorphic_pseudogene or processed_transcript
grep -v "polymorphic_pseudogene" gencode.v24.primary_assembly.annotation.lincNprot.gtf > gencode.v24.primary_assembly.annotation.lincNprot_temp.gtf
grep -v "processed_transcript" gencode.v24.primary_assembly.annotation.lincNprot_temp.gtf > gencode.v24.primary_assembly.annotation.lincNprot_fix.gtf
rm gencode.v24.primary_assembly.annotation.lincNprot_temp.gtf

 # now use cgat to get the longest transcript per gene
cgat gtf2gtf --method=filter --filter-method=longest-transcript -I gencode.v24.primary_assembly.annotation.lincNprot_fix.gtf > gencode.v24.primary_assembly.annotation.lincNprot_fix.filterLongestTrans.gtf

 #sanity checks
awk '$3=="gene"' gencode.v24.primary_assembly.annotation.lincNprot_fix.gtf | wc -l # 27,518
awk '$3=="transcript"' gencode.v24.primary_assembly.annotation.lincNprot_fix.filterLongestTrans.gtf  | wc -l # 27,518

 #create bed file from gtf
module load ucsc/080613 # aka kent tools
gtfToGenePred gencode.v24.primary_assembly.annotation.lincNprot_fix.filterLongestTrans.gtf gencode.v24.primary_assembly.annotation.lincNprot_fix.filterLongestTrans.genePred
genePredToBed gencode.v24.primary_assembly.annotation.lincNprot_fix.filterLongestTrans.genePred gencode.v24.primary_assembly.annotation.lincNprot_fix.filterLongestTrans.bed
```

To calculate TIN, I used RSeQC's tin.py script to run the analysis on each of the bam files.
```bash
for f1 in *bam;
do
SAMPLE="$(echo $f1 | cut -c1-3)"
echo $SAMPLE
qsub -S /bin/bash -N tin_$SAMPLE -cwd -l h_vmem=10G -j y << EOF
module load RSeQC/2.6.4/python.2.7.8
tin.py --input=${f1} --minCov=10 --sample-size=100 --refgene=/home/greally-lab/indexes/hg38/GenCode/gencode.v24.primary_assembly.annotation.lincNprot_fix.filterLongestTrans.bed
EOF
done
```

I wanted to look at the RNA integrity (mRIN) of each of the libraries. I used mRIN software to look at the quality of the libraries.

```bash
 # convert the bam files to bed
for f1 in *.bam;
do
SAMPLE="$(echo ${f1} | cut -c 1-3)"
echo $SAMPLE
qsub -S /bin/bash -N BamToBed_${SAMPLE} -cwd -l h_vmem=7G -j y << EOF
module load bedtools2/2.26.0/gcc.4.4.7
bedtools bamtobed -i $f1 > ${f1%.bam}.bed
EOF
done

mkdir mRIN

 # run the mRIN analysis, starting with conversion of bed to bedgraph
module load mRIN/012218
for f1 in *bed;
do
SAMPLE="$(echo ${f1} | cut -c 1-3)"
echo $SAMPLE
qsub -S /bin/bash -N tag2profile_${SAMPLE} -cwd -l h_vmem=7G -j y << EOF
module load mRIN/012218
module load R/3.4.0/gcc.4.7.4
PERL5LIB=~/Programs/mrin-master/czplib-master
perl $(which tag2profile.pl) -big -exact -of bedgraph -v $f1 mRIN/${SAMPLE}.bedGraph
EOF
done

cd mRIN
mkdir cdf #we will save them in a separate directory to be organized
mkdir ks #we will save them in a separate directory to be organized

module load mRIN/012218
for f1 in *.bedGraph;
do
SAMPLE="$(echo ${f1} | cut -c 1-3)"
echo $SAMPLE
qsub -S /bin/bash -N ks_${SAMPLE} -cwd -l h_vmem=7G -j y << EOF
module load mRIN/012218
module load R/3.4.0/gcc.4.7.4
PERL5LIB=~/Programs/mrin-master/czplib-master
perl $(which gen_transcript_cdf.pl) -v /home/greally-lab/indexes/hg38/GenCode/gencode.v24.primary_assembly.annotation.lincNprot_fix.filterLongestTrans.bed ${SAMPLE}.bedGraph cdf/${SAMPLE}.cdf.bedGraph
perl $(which ks_test_uniform.pl) -v cdf/${SAMPLE}.cdf.bedGraph ks/${SAMPLE}.ks.txt
EOF
done

 # making the .conf file
ls -lah ks/ | awk 'NR>3 {print $9}' | awk -F . '{print "ks/"$0"\t"$1}'  >  all_samples.conf

 # generate the mRIN for genes and samples
qsub -S /bin/bash -N ks_all -cwd -l h_vmem=7G -j y << EOF
module load mRIN/012218
module load R/3.4.0/gcc.4.7.4
PERL5LIB=~/Programs/mrin-master/czplib-master
perl $(which gen_ks_matrix.pl) -v -base ./ --min-avg-cov 10 -v all_samples.conf all_samples.KS.mat.txt
Rscript ~/Programs/mrin-master/cal_mrin.R -k all_samples.KS.mat.txt -s 0.05 -e 0.5 -b a -m out.mRIN.txt -v 1 -G out.GIS.txt
EOF
```

I then looked for samples that were flagged as being poor quality in more than 1 of the methods. The following samples were identified: A19, A21?, A47, B03, B04, B17, B25, and B52.

<a name="Covs">

3. Preparing the covariate file

I'm going to be using the first 3 PCs from the genotype data and then a few of the PCs from the expression data. The COV file contains the covariate data in simple TXT format.

1. The file is TAB delimited
2. First row gives the sample ID and each additional one corresponds to a single covariate
3. First column gives the covariate ID and each additional one corresponds to a sample
4. The file should have S+1 rows and C+1 columns where S and C are the numbers of samples and covariates, respectively.

Both quantitative and qualitative covariates are supported. Quantitative covariates are assumed when only numeric values are provided. Qualitative covariates are assumed when only non-numeric values are provided. In practice, qualitative covariates with F factors are converted in F-1 binary covariates.

<a name="gtf">

4. Preparing the GTF (only going to look at protein coding genes and linc)

I used the previously filtered GTF created as an intermediate product for the [bam QC](#bamqc) `  gencode.v24.primary_assembly.annotation.lincNprot_fix.gtf`. I had decided to use the gencode v24 gtf because of the chr issue and its what the creators of QTLtools use, so there should not be an error when I try to quantify the reads. Now I just need to gz the file.

```bash
gzip gencode.v24.primary_assembly.annotation.lincNprot_fix.gtf
```

The file is found here: `/home/greally-lab/indexes/hg38/GenCode/gencode.v24.primary_assembly.annotation.lincNprot_fix.gtf.gz`

<a name="genoprep">

5. Prepare the genotype data

This was performed in the [Genotyping_QC_calling.md](./Genotyping_QC_calling.md)

The finalized VCF file that I'm using is `/home/greally-lab/Deepa_Andrew/eQTL/Genotype2/SplX_converted_noHH_eQTL_snpfilt_hg38_sort.vcf.gz`.


The major step of matching the RNAseq files to the genotyped individuals was performed in the [Genotyping_QC_calling.md](./Genotyping_QC_calling.md). All of the mismatched bams need to be properly switched and the unused samples (poor quality and related) need to be removed.

With the finalized file I ran PCA on it.

```bash
qsub -S /bin/bash -N PCA_geno -cwd -l h_vmem=10G -j y << EOF
module load QTLtools/1.1/gcc.4.7.4
QTLtools pca --vcf SplX_converted_noHH_eQTL_snpfilt_hg38_sort.vcf.gz --scale --center --maf 0.1 --distance 50000 --out genotypes.pca
EOF
 # I call maf .1 but the variants were previously reduced to this.
 # --distance 50000 to only consider variant sites separated by at least 50kb
 	# how does QTLtools choose which SNPs within 50kb to use?
```

Next I analyzed the PCA data in R (see [PCA_heatmap_analysis2.R](./PCA_heatmap_analysis2.R))

<a name="quant">

6. Quantify gene expression / Preparing the phenotype file

I need to quantify gene expression and create a quantification matrix by using the QTLtools quan mode.

Using aforementioned gencode gtf, I run the quan command as follows in `/home/greally-lab/Deepa_Andrew/eQTL/Bams/filt_bams`:

```bash
for f1 in *.bam;
do
SAMPLE="$(echo ${f1} | cut -c1-3)"
qsub -S /bin/bash -N q2_${SAMPLE} -cwd -l h_vmem=10G -j y << EOF
module load QTLtools/1.1/gcc.4.7.4
QTLtools quan --bam $f1 --gtf /home/greally-lab/indexes/hg38/GenCode/gencode.v24.primary_assembly.annotation.lincNprot_fix.gtf.gz --samples $SAMPLE --out-prefix ../quants2/$SAMPLE --filter-mapping-quality 0 --filter-mismatch 5 --rpkm --no-merge
EOF
done
```

Next I created the quantification matrix. First, I extracted the annotation columns using cut -f1-6 on one of the files

```bash
cut -f1-6  A01.exon.rpkm.bed > exon_rpkm_annotation.txt
cut -f1-6  A01.gene.rpkm.bed > gene_rpkm_annotation.txt
```

Then I moved all of the unused bam files into a separate folder `unused_bams`. Here are the unused bams placed in `unused_bams.txt`:
A03
A07
A11
A14
A18
A19
A31
A41
A47
A52
A58
B03
B04
B11
B14
B17
B23
B25
B33
B48
B52

```bash
mkdir unused_bams

 # move all of the unused bams
while IFS='' read -r sample; do
echo $sample
mv ${sample}* unused_bams/
 # mv *${sample} unused_bams/
done < unused_bams.txt
```

Extract the per sample quantifications using cut -f7 for each file separately

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

Build the entire quantification matrix using the paste to merge the annotations with all quantifications.

```bash
paste exon_rpkm_annotation.txt sample_exon_rpkm_* > matrix_exon_rpkm.txt
paste gene_rpkm_annotation.txt sample_gene_rpkm_* > matrix_gene_rpkm.txt

 # Sanity check - answer should correspond to the number of samples +6
 awk '{print NF}' matrix_exon_rpkm.txt | sort -nu | tail -n 1 # 97 = 91+6
 awk '{print NF}' matrix_gene_rpkm.txt | sort -nu | tail -n 1 # 97 = 91+6

```

Filter all poorly quantified exons/genes using awk for example

```R
 # filtering the genes that are not well expressed in majority of samples
options(stringsAsFactors = F)
options(scipen = 999)
 # load libraries
library(data.table)
library(genefilter)
 # read in quant matrix
mat <- fread("matrix_gene_rpkm.txt", header=T, sep="\t")
mat <- as.data.frame(mat)
dim(mat) # 27518    97
head(mat,1)

 # Table S3 in Geavardis paper says the following:
 # "All eQTL counts are for autosomal genes, with a filter of
 # quantification in >90% of samples for genes, exons, transcripts"
table(mat$`#chr`)
non_auto_chr <- c("M", "X", "Y", "GL", "KI")
idx_nonAuto <- grep(pattern = paste(non_auto_chr, collapse = "|"),
                    x = mat$`#chr`)
length(idx_nonAuto) # 1138
mat_auto <- mat[-idx_nonAuto,]
table(mat_auto$`#chr`)
dim(mat_auto) # 26380    97

mat3 <- NULL
for (i in 1:nrow(mat_auto)){
  if (mean(mat_auto[i,-c(1:6)]==0)<.1){
    mat3 <- rbind(mat3, mat_auto[i,])
  }
}
dim(mat3) # 12960    97
head(mat3,1)
 # this number of genes is similar to the 13,703 genes used in the Geavardis
 # paper for eQTL analysis (see Table S3).
hist(rowMeans(log(mat3[,-c(1:6)]+2)))
means_row <- rowMeans(mat3[,-c(1:6)])
summary(means_row)

 # Now I'll make the sample swaps with the IDs
 # reading in key generated for BAM QC
sample_key <- read.csv("../../Bams/sample_prot_coding.csv")
nrow(sample_key)
head(sample_key)
df_colnames <- data.frame(bam_name=colnames(mat3)[-c(1:6)])
df_colnames_key <- merge(x = df_colnames, y = sample_key, by.x = "bam_name", by.y = "Bam_sample")
head(df_colnames_key)
nrow(df_colnames_key)
colnames(mat3)[-c(1:6)] <- df_colnames_key$Act_sample

 # Now I need to make sure that I am calling my samples the same name and get
 # them in the same order
geno_pca <- read.table("../../Genotype2/genotypes.pca.pca", header=T)
head(geno_pca)
geno_IDs <- substring(colnames(geno_pca)[-1], 2, 8)

idx_match <- match(colnames(mat3[,-c(1:6)]), do.call(rbind, strsplit(geno_IDs, "_"))[,2])
same_ids <- which(colnames(mat3) %in% do.call(rbind, strsplit(geno_IDs, "_"))[,2])
length(same_ids)
colnames(mat3)[-same_ids]

idx_match2 <- match(do.call(rbind, strsplit(geno_IDs, "_"))[,2],colnames(mat3[,-c(1:6)]))

 # add the correct genotype names to expression matrix
mat4 <- mat3
colnames(mat4)[-c(1:6)] <- geno_IDs[idx_match]
order(colnames(mat4)[-c(1:6)])
 # reorder the expression columns to match the order in the genotype pca file.
mat4_lite <- mat4[,-c(1:6)]
mat4_lite <- mat4_lite[,idx_match2]
head(mat4_lite)

 # combine the annotation of genes with the reordered sample columns.
mat5 <- cbind(mat4[,c(1:6)], mat4_lite)
dim(mat5) # [1] 12960    97
head(mat5,1)

 # write out the expression table
write.table(mat5, "Filtered_matrix_gene_rpkm.bed", col.names = T, row.names = F,
            append = F, quote = F, sep = "\t")

 # make a file for cell proportion calling
cell_prop <- read.table("Filtered_matrix_gene_rpkm.bed", header=T)
cell_prop <- cell_prop[,-c(1:3,5:6)]
colnames(cell_prop)
colnames(cell_prop)[-1] <- do.call(rbind, strsplit(colnames(cell_prop)[-1], "_"))[,2]
cell_prop[,1] <- do.call(rbind, strsplit(cell_prop[,1], "[.]"))[,1]

write.table(cell_prop, "Filtered_matrix_gene_rpkm_cellProp.txt", col.names = T, row.names = F,
            append = F, quote = F, sep = "\t")
```

Make into proper phenotype file
```bash
 # mv Filtered_matrix_gene_rpkm.txt Filtered_matrix_gene_rpkm.bed
 # if needed:
 # vi Filtered_matrix_gene_rpkm.bed #  add a "#" in front of the first line aka comment it out
 # awk '{if($0 !~ /^#/) print "chr"$0; else print $0}' Filtered_matrix_gene_rpkm.bed > Filtered_matrix_gene_rpkm_chr.bed

 # bgzip and tabix
module load htslib/1.2.1/gcc.4.4.7
bgzip Filtered_matrix_gene_rpkm.bed && tabix -p bed Filtered_matrix_gene_rpkm.bed.gz
```

Finally, perform PCA analysis on the expression data
```bash
 # Gene Expression PCA
qsub -S /bin/bash -N PCA_gene_exp -cwd -l h_vmem=10G -j y << EOF
module load QTLtools/1.1/gcc.4.7.4
QTLtools pca --bed Filtered_matrix_gene_rpkm.bed.gz --scale --center --out genes_90percent
EOF
```

<a name="PCA">

7. Explore Expression PCA

Using the sample metadata, observe which covariates are contributing the most to the various expression PCs. Visualizing using PCA heatmaps.


Expression PC: 1-10?,+13 w/o PC4?
Geno PC: 1,2
cell prop?


8. Nominal eQTL passes

I needed to try out different combinations of covariates based on what I saw in the PCA heatmaps to make educated guesses on which PCs to include to optimize the eQTL calling.

First I selected:
Expression PC: 1-10,13
Geno PC: 1,2

```bash
cd /home/greally-lab/Deepa_Andrew/eQTL/Genotype2
awk 'NR<4' ../Genotype2/genotypes.pca.pca > Cov_geno_2PCs.txt

cd /home/greally-lab/Deepa_Andrew/eQTL/Nominal2

 # creating cov file for first 10 gene expression PCs plus the first 2 geno PCs
awk 'NR>1 && NR<12' /home/greally-lab/Deepa_Andrew/eQTL/Bams/quants2/genes_90percent.pca | cat /home/greally-lab/Deepa_Andrew/eQTL/Genotype2/Cov_geno_2PCs.txt - > Cov_2PC_geno_10PC_gexprs.txt

 # Adding the 13th PC
awk 'NR==14' /home/greally-lab/Deepa_Andrew/eQTL/Bams/quants2/genes_90percent.pca | cat Cov_2PC_geno_10PC_gexprs.txt - > Cov_2PC_geno_10_13PC_gexprs.txt

 # the files need to be gzip'd
bgzip Cov_2PC_geno_10_13PC_gexprs.txt
bgzip Cov_2PC_geno_10PC_gexprs.txt

cd nominal_10_13
 # finding nominal
for j in $(seq 1 200);
do
qsub -S /bin/bash -N eQTL_${j}_chunk_nom -cwd -l h_vmem=5.6G -j y << EOF
module load QTLtools/1.1/gcc.4.7.4
QTLtools cis --vcf ../../Genotype2/SplX_converted_noHH_eQTL_snpfilt_hg38_sort.vcf.gz --bed ../../Bams/quants2/Filtered_matrix_gene_rpkm.bed.gz --cov ../Cov_2PC_geno_10_13PC_gexprs.txt.gz --nominal 0.01  --chunk $j 200 --out chunk_$j\_200.txt --normal --window 1000000
EOF
done

cat chunk*_200.txt > nominals_full.txt
wc -l nominals_full.txt
 # 54,216 eQTLs (pval <.01) found in nominals_full.txt

 # top eQTL for gene
awk '$14==1'  nominals_full.txt | wc -l
 # 10,714 top eQTLs for gene; so only 10828 genes had a variant associated with it.
```

Now let's run the nominal pass without the 13th PC

```bash
cd ../
mkdir nominal_10
cd nominal_10

 # running the nominal
for j in $(seq 1 200);
do
qsub -S /bin/bash -N eQTL_${j}_chunk_nom -cwd -l h_vmem=5.6G -j y << EOF
module load QTLtools/1.1/gcc.4.7.4
QTLtools cis --vcf ../../Genotype2/SplX_converted_noHH_eQTL_snpfilt_hg38_sort.vcf.gz --bed ../../Bams/quants2/Filtered_matrix_gene_rpkm.bed.gz --cov ../Cov_2PC_geno_10PC_gexprs.txt.gz --nominal 0.01  --chunk $j 200 --out chunk_$j\_200.txt --normal --window 1000000
EOF
done

cat chunk*_200.txt > nominals_full.txt
wc -l nominals_full.txt
 # 54,448 eQTLs (pval <.01) found in nominals_full.txt

 # top eQTL for gene
awk '$14==1'  nominals_full.txt | wc -l
 # 10,764 top eQTLs for gene; so only 10764 genes had a variant associated with it.
```

I'm worried that I should be testing using permutation tests so let's see if there's a similar difference.

9. Permutation Passes

Let's examine the number of eQTLs called when using permutation.

As a side note, I am calling using the `--normal` argument since the software uses linear modeling which assumes normally distributed phenotypes. Therefore, the argument is used to force the phenotype to match a normal distribution.

Using first 2 PCs for geno and 10 for gene expression.
```bash
cd /home/greally-lab/Deepa_Andrew/eQTL/Permutation2
mkdir perm_pc10
cd perm_pc10

 # call the permutations
for j in $(seq 1 200);
do
qsub -S /bin/bash -N eQTL_${j}_chunk -cwd -R y -l h_vmem=5.6G -j y << EOF
module load QTLtools/1.1/gcc.4.7.4
QTLtools cis --vcf ../../Genotype2/SplX_converted_noHH_eQTL_snpfilt_hg38_sort.vcf.gz --bed ../../Bams/quants2/Filtered_matrix_gene_rpkm.bed.gz --cov ../../Nominal2/Cov_2PC_geno_10PC_gexprs.txt.gz --permute 1000  --chunk $j 200 --out perms_$j\_200.txt --normal --window 1000000 --seed 123456
EOF
done

 # combine the permutation files together
cat perms_*200.txt > perms_full.txt
wc -l perms_full.txt
 # 12960 perms_full.txt

module load R/3.4.0/gcc.4.7.4
Rscript plot_beta.R perms_full.txt G_12_E_10
 # Found 452 significant QTls.
```

Using first 2 PCs for geno and 10 plus PC 13 for gene expression.
```bash
cd /home/greally-lab/Deepa_Andrew/eQTL/Permutation2
mkdir perm_pc10_13
cd perm_pc10_13

 # call the permutations
for j in $(seq 1 200);
do
qsub -S /bin/bash -N eQTL_${j}_chunk -cwd -R y -l h_vmem=5.6G -j y << EOF
module load QTLtools/1.1/gcc.4.7.4
QTLtools cis --vcf ../../Genotype2/SplX_converted_noHH_eQTL_snpfilt_hg38_sort.vcf.gz --bed ../../Bams/quants2/Filtered_matrix_gene_rpkm.bed.gz --cov ../../Nominal2/Cov_2PC_geno_10_13PC_gexprs.txt.gz --permute 1000  --chunk $j 200 --out perms_$j\_200.txt --normal --window 1000000 --seed 123456
EOF
done

 # combine the permutation files together
cat perms_*200.txt > perms_full.txt
wc -l perms_full.txt
 # 12960 perms_full.txt

cp ../perm_pc10/plot_beta.R .
 # determine the number of significant genes
module load R/3.4.0/gcc.4.7.4
Rscript plot_beta.R perms_full.txt G_12_E_10_13
 # Found 446 significant QTls.
```

Using first 2 PCs for geno and 10 minus PC 4 for gene expression.
```bash
cd /home/greally-lab/Deepa_Andrew/eQTL/Permutation2
cp ../Nominal2/Cov_2PC_geno_10PC_gexprs.txt.gz .
gunzip -d Cov_2PC_geno_10PC_gexprs.txt.gz
awk 'NR!=7' Cov_2PC_geno_10PC_gexprs.txt > Cov_2PC_geno_10_no4_PC_gexprs.txt
bgzip Cov_2PC_geno_10_no4_PC_gexprs.txt

mkdir perm_pc10_no4
cd perm_pc10_no4

 # call the permutations
for j in $(seq 1 200);
do
qsub -S /bin/bash -N eQTL_${j}_chunk -cwd -R y -l h_vmem=5.6G -j y << EOF
module load QTLtools/1.1/gcc.4.7.4
QTLtools cis --vcf ../../Genotype2/SplX_converted_noHH_eQTL_snpfilt_hg38_sort.vcf.gz --bed ../../Bams/quants2/Filtered_matrix_gene_rpkm.bed.gz --cov ../Cov_2PC_geno_10_no4_PC_gexprs.txt.gz --permute 1000  --chunk $j 200 --out perms_$j\_200.txt --normal --window 1000000 --seed 123456
EOF
done

 # combine the permutation files together
cat perms_*200.txt > perms_full.txt
wc -l perms_full.txt
 # 12960 perms_full.txt

cp ../perm_pc10/plot_beta.R .
 # determine the number of significant genes
module load R/3.4.0/gcc.4.7.4
Rscript plot_beta.R perms_full.txt G_12_E_10_13
 # Found 417 significant QTls.
```

I'll generate a covariate file containing 1-5 PCs + additional

```bash
cd /home/greally-lab/Deepa_Andrew/eQTL/Permutation2
module load htslib/1.2.1/gcc.4.4.7

awk 'NR<9' Cov_2PC_geno_10PC_gexprs.txt > Cov_2PC_geno_1-5_PC_gexprs.txt
bgzip Cov_2PC_geno_1-5_PC_gexprs.txt
awk 'NR<10' Cov_2PC_geno_10PC_gexprs.txt > Cov_2PC_geno_1-6_PC_gexprs.txt
bgzip Cov_2PC_geno_1-6_PC_gexprs.txt
awk 'NR<11' Cov_2PC_geno_10PC_gexprs.txt > Cov_2PC_geno_1-7_PC_gexprs.txt
bgzip Cov_2PC_geno_1-7_PC_gexprs.txt
awk 'NR<12' Cov_2PC_geno_10PC_gexprs.txt > Cov_2PC_geno_1-8_PC_gexprs.txt
bgzip Cov_2PC_geno_1-8_PC_gexprs.txt
awk 'NR<13' Cov_2PC_geno_10PC_gexprs.txt > Cov_2PC_geno_1-9_PC_gexprs.txt
bgzip Cov_2PC_geno_1-9_PC_gexprs.txt
```

Now to call QTLs with the various PCs

Using first 2 geno and 1-5 exprs PCs
```bash
mkdir perm_pc5
cd perm_pc5

 # call the permutations
for j in $(seq 1 200);
do
qsub -S /bin/bash -N eQTL_${j}_chunk -cwd -R y -l h_vmem=5.6G -j y << EOF
module load QTLtools/1.1/gcc.4.7.4
QTLtools cis --vcf ../../Genotype2/SplX_converted_noHH_eQTL_snpfilt_hg38_sort.vcf.gz --bed ../../Bams/quants2/Filtered_matrix_gene_rpkm.bed.gz --cov ../Cov_2PC_geno_1-5_PC_gexprs.txt.gz --permute 1000  --chunk $j 200 --out perms_$j\_200.txt --normal --window 1000000 --seed 123456
EOF
done

 # combine the permutation files together
cat perms_*200.txt > perms_full.txt
wc -l perms_full.txt
 # 12960 perms_full.txt

cp ../perm_pc10/plot_beta.R .
 # determine the number of significant genes
module load R/3.4.0/gcc.4.7.4
Rscript plot_beta.R perms_full.txt G_12_E_10_13
 # Found 418 significant QTls.
```

```bash
cd /home/greally-lab/Deepa_Andrew/eQTL/Permutation2
mkdir perm_pc6
cd perm_pc6

 # call the permutations
for j in $(seq 1 200);
do
qsub -S /bin/bash -N eQTL_${j}_chunk -cwd -R y -l h_vmem=5.6G -j y << EOF
module load QTLtools/1.1/gcc.4.7.4
QTLtools cis --vcf ../../Genotype2/SplX_converted_noHH_eQTL_snpfilt_hg38_sort.vcf.gz --bed ../../Bams/quants2/Filtered_matrix_gene_rpkm.bed.gz --cov ../Cov_2PC_geno_1-6_PC_gexprs.txt.gz --permute 1000  --chunk $j 200 --out perms_$j\_200.txt --normal --window 1000000 --seed 123456
EOF
done

 # combine the permutation files together
cat perms_*200.txt > perms_full.txt
wc -l perms_full.txt
 # 12960 perms_full.txt

cp ../perm_pc10/plot_beta.R .
 # determine the number of significant genes
module load R/3.4.0/gcc.4.7.4
Rscript plot_beta.R perms_full.txt G_12_E_123456
 # Found 418 significant QTls.
```

```bash
cd /home/greally-lab/Deepa_Andrew/eQTL/Permutation2
mkdir perm_pc7
cd perm_pc7

 # call the permutations
for j in $(seq 1 200);
do
qsub -S /bin/bash -N eQTL_${j}_chunk -cwd -R y -l h_vmem=5.6G -j y << EOF
module load QTLtools/1.1/gcc.4.7.4
QTLtools cis --vcf ../../Genotype2/SplX_converted_noHH_eQTL_snpfilt_hg38_sort.vcf.gz --bed ../../Bams/quants2/Filtered_matrix_gene_rpkm.bed.gz --cov ../Cov_2PC_geno_1-7_PC_gexprs.txt.gz --permute 1000  --chunk $j 200 --out perms_$j\_200.txt --normal --window 1000000 --seed 123456
EOF
done

 # combine the permutation files together
cat perms_*200.txt > perms_full.txt
wc -l perms_full.txt
 # 12960 perms_full.txt

cp ../perm_pc10/plot_beta.R .
 # determine the number of significant genes
module load R/3.4.0/gcc.4.7.4
Rscript plot_beta.R perms_full.txt G_12_E_1234567
 # Found 439 significant QTls.
```

```bash
cd /home/greally-lab/Deepa_Andrew/eQTL/Permutation2
mkdir perm_pc8
cd perm_pc8

 # call the permutations
for j in $(seq 1 200);
do
qsub -S /bin/bash -N eQTL_${j}_chunk -cwd -R y -l h_vmem=5.6G -j y << EOF
module load QTLtools/1.1/gcc.4.7.4
QTLtools cis --vcf ../../Genotype2/SplX_converted_noHH_eQTL_snpfilt_hg38_sort.vcf.gz --bed ../../Bams/quants2/Filtered_matrix_gene_rpkm.bed.gz --cov ../Cov_2PC_geno_1-8_PC_gexprs.txt.gz --permute 1000  --chunk $j 200 --out perms_$j\_200.txt --normal --window 1000000 --seed 123456
EOF
done

 # combine the permutation files together
cat perms_*200.txt > perms_full.txt
wc -l perms_full.txt
 # 12960 perms_full.txt

cp ../perm_pc10/plot_beta.R .
 # determine the number of significant genes
module load R/3.4.0/gcc.4.7.4
Rscript plot_beta.R perms_full.txt G_12_E_12345678
 # Found 446 significant QTls.
```

```bash
cd /home/greally-lab/Deepa_Andrew/eQTL/Permutation2
mkdir perm_pc9
cd perm_pc9

 # call the permutations
for j in $(seq 1 200);
do
qsub -S /bin/bash -N eQTL_${j}_chunk -cwd -R y -l h_vmem=5.6G -j y << EOF
module load QTLtools/1.1/gcc.4.7.4
QTLtools cis --vcf ../../Genotype2/SplX_converted_noHH_eQTL_snpfilt_hg38_sort.vcf.gz --bed ../../Bams/quants2/Filtered_matrix_gene_rpkm.bed.gz --cov ../Cov_2PC_geno_1-9_PC_gexprs.txt.gz --permute 1000  --chunk $j 200 --out perms_$j\_200.txt --normal --window 1000000 --seed 123456
EOF
done

 # combine the permutation files together
cat perms_*200.txt > perms_full.txt
wc -l perms_full.txt
 # 12960 perms_full.txt

cp ../perm_pc10/plot_beta.R .
 # determine the number of significant genes
module load R/3.4.0/gcc.4.7.4
Rscript plot_beta.R perms_full.txt G_12_E_123456789
 # Found 452 significant QTls.
```

maybe I can get rid of PC 6 and 10? what about adding 11?

```bash
cd /home/greally-lab/Deepa_Andrew/eQTL/Permutation2
module load htslib/1.2.1/gcc.4.4.7
gunzip -d Cov_2PC_geno_1-9_PC_gexprs.txt.gz
awk 'NR!=9' Cov_2PC_geno_1-9_PC_gexprs.txt > Cov_2PC_geno_1-9_no6_PC_gexprs.txt
bgzip Cov_2PC_geno_1-9_no6_PC_gexprs.txt

 # Adding the 13th PC
awk 'NR==12' /home/greally-lab/Deepa_Andrew/eQTL/Bams/quants2/genes_90percent.pca | cat Cov_2PC_geno_10PC_gexprs.txt - > Cov_2PC_geno_11PC_gexprs.txt
bgzip Cov_2PC_geno_11PC_gexprs.txt
```


```bash
cd /home/greally-lab/Deepa_Andrew/eQTL/Permutation2
mkdir perm_pc9_no6
cd perm_pc9_no6

 # call the permutations
for j in $(seq 1 200);
do
qsub -S /bin/bash -N eQTL_${j}_chunk -cwd -R y -l h_vmem=5.6G -j y << EOF
module load QTLtools/1.1/gcc.4.7.4
QTLtools cis --vcf ../../Genotype2/SplX_converted_noHH_eQTL_snpfilt_hg38_sort.vcf.gz --bed ../../Bams/quants2/Filtered_matrix_gene_rpkm.bed.gz --cov ../Cov_2PC_geno_1-9_no6_PC_gexprs.txt.gz --permute 1000  --chunk $j 200 --out perms_$j\_200.txt --normal --window 1000000 --seed 123456
EOF
done

 # combine the permutation files together
cat perms_*200.txt > perms_full.txt
wc -l perms_full.txt
 # 12960 perms_full.txt

cp ../perm_pc10/plot_beta.R .
 # determine the number of significant genes
module load R/3.4.0/gcc.4.7.4
Rscript plot_beta.R perms_full.txt G_12_E_12345789
 # Found 447 significant QTls.
```

```bash
cd /home/greally-lab/Deepa_Andrew/eQTL/Permutation2
mkdir perm_pc11
cd perm_pc11

 # call the permutations
for j in $(seq 1 200);
do
qsub -S /bin/bash -N eQTL_${j}_chunk -cwd -R y -l h_vmem=5.6G -j y << EOF
module load QTLtools/1.1/gcc.4.7.4
QTLtools cis --vcf ../../Genotype2/SplX_converted_noHH_eQTL_snpfilt_hg38_sort.vcf.gz --bed ../../Bams/quants2/Filtered_matrix_gene_rpkm.bed.gz --cov ../Cov_2PC_geno_11PC_gexprs.txt.gz --permute 1000  --chunk $j 200 --out perms_$j\_200.txt --normal --window 1000000 --seed 123456
EOF
done

 # combine the permutation files together
cat perms_*200.txt > perms_full.txt
wc -l perms_full.txt
 # 12960 perms_full.txt

cp ../perm_pc10/plot_beta.R .
 # determine the number of significant genes
module load R/3.4.0/gcc.4.7.4
Rscript plot_beta.R perms_full.txt G_12_E_123457801011
 # Found 446 significant QTls.
```

what about 3 Geno PCs?

```bash
cd /home/greally-lab/Deepa_Andrew/eQTL/Permutation2
awk 'NR<5' /home/greally-lab/Deepa_Andrew/eQTL/Genotype2/genotypes.pca.pca > Cov_geno_3PCs.txt

 # creating cov file for first 10 gene expression PCs plus the first 3 geno PCs
awk 'NR>1 && NR<12' /home/greally-lab/Deepa_Andrew/eQTL/Bams/quants2/genes_90percent.pca | cat Cov_geno_3PCs.txt - > Cov_3PC_geno_10PC_gexprs.txt

awk 'NR<14' Cov_3PC_geno_10PC_gexprs.txt > Cov_3PC_geno_9PC_gexprs.txt
bgzip Cov_3PC_geno_9PC_gexprs.txt
```

```bash
cd /home/greally-lab/Deepa_Andrew/eQTL/Permutation2
mkdir perm_pc9_3g
cd perm_pc9_3g

 # call the permutations
for j in $(seq 1 200);
do
qsub -S /bin/bash -N eQTL_${j}_chunk -cwd -R y -l h_vmem=5.6G -j y << EOF
module load QTLtools/1.1/gcc.4.7.4
QTLtools cis --vcf ../../Genotype2/SplX_converted_noHH_eQTL_snpfilt_hg38_sort.vcf.gz --bed ../../Bams/quants2/Filtered_matrix_gene_rpkm.bed.gz --cov ../Cov_3PC_geno_9PC_gexprs.txt.gz --permute 1000  --chunk $j 200 --out perms_$j\_200.txt --normal --window 1000000 --seed 123456
EOF
done

 # combine the permutation files together
cat perms_*200.txt > perms_full.txt
wc -l perms_full.txt
 # 12960 perms_full.txt

cp ../perm_pc10/plot_beta.R .
 # determine the number of significant genes
module load R/3.4.0/gcc.4.7.4
Rscript plot_beta.R perms_full.txt G_123_E_12345789
 # Found 430 significant QTls.
```

It seems that using the first 10 expression PCs and the first 2 genotype PCs gives the most eQTLs in the result. So I ran QTL analysis through 10,000 permutations to generate a "final" eQTL list.

```bash
cd /home/greally-lab/Deepa_Andrew/eQTL/Permutation2
mkdir perm_pc10_manyperm
cd perm_pc10_manyperm

 # call the permutations
for j in $(seq 1 200);
do
qsub -S /bin/bash -N eQTL_${j}_chunk -cwd -R y -l h_vmem=5.6G -j y << EOF
module load QTLtools/1.1/gcc.4.7.4
QTLtools cis --vcf ../../Genotype2/SplX_converted_noHH_eQTL_snpfilt_hg38_sort.vcf.gz --bed ../../Bams/quants2/Filtered_matrix_gene_rpkm.bed.gz --cov ../../Nominal2/Cov_2PC_geno_10PC_gexprs.txt.gz --permute 10000  --chunk $j 200 --out perms_$j\_200.txt --normal --window 1000000 --seed 123456
EOF
done

 # combine the permutation files together
cat perms_*200.txt > eQTL_final1.txt
wc -l eQTL_final1.txt
 # 12960 eQTL_final1.txt

module load R/3.4.0/gcc.4.7.4
Rscript plot_beta.R eQTL_final1.txt G_12_E_10_manyperm
 # Found 443 significant QTls.

gzip eQTL_final1.txt
Rscript RunFDR_cis.R eQTL_final1.txt.gz 0.05 G_12_E_10_final

wc -l G_12_E_10_final.significant.txt
 #443 G_12_E_10_final.significant.txt
```

Now we have a list of 443 variants (top variants of genes), which pass the significance threshold. Now we'll make a conditional pass to find other significant varaints, which are independent of the top variant

```bash
cd ../
mkdir perm_pc10_cond
cd perm_pc10_cond

for j in $(seq 1 200);
do
qsub -S /bin/bash -N eQTL_${j}_chunk -cwd -R y -l h_vmem=5.6G -j y << EOF
module load QTLtools/1.1/gcc.4.7.4
QTLtools cis --vcf ../../Genotype2/SplX_converted_noHH_eQTL_snpfilt_hg38_sort.vcf.gz --bed ../../Bams/quants2/Filtered_matrix_gene_rpkm.bed.gz --cov ../../Nominal2/Cov_2PC_geno_10PC_gexprs.txt.gz --mapping ../perm_pc10_manyperm/G_12_E_10_final.thresholds.txt --normal --window 1000000 --chunk $j 200 --out cond_$j\_200.txt
EOF
done

cat cond_*200.txt > eQTL_final2.txt
wc -l eQTL_final2.txt
 # 1582 eQTL_final2.txt
 # take only those that pass backward significance
awk '$20==1' eQTL_final2.txt > eQTL_final3.txt

 #how many SNPs
wc -l eQTL_final3.txt
 # 1574 eQTL_final3.txt

 # how many genes?
awk '$19==1' eQTL_final3.txt | wc -l #449
awk '{print $1}' eQTL_final3.txt | uniq | wc -l #443, which means there's a tie somewhere?
```

So let's examine these 1574 eQTLs and 443 eGenes.

Sanity check that the perm genes are in the conditional pass. go back to perm directory, `cd ../perm_pc10_manyperm`


```R
fpass <- read.table("G_12_E_10_final.significant.txt")
spass <- read.table("../perm_pc10_cond/eQTL_final3.txt")

sum(fpass$V1 %in% spass$V1) #443
```

10. Exploring eQTLs

First let's jsut see if any of the top 11 genes came up.
VAV2	ENSG00000160293
DOCK5	ENSG00000147459
PAK3	ENSG00000077264 # not in the filtered genes b/c on X
PLD1	ENSG00000075651
CDC42EP4	ENSG00000179604
CDC42PBB	ENSG00000070831 # couldn't find an ENS ID specific to it.
MLK3	ENSG00000173327
CDC42	ENSG00000070831
MEN1	ENSG00000133895
PIK3CA	 ENSG00000121879
PIK3CD	ENSG00000171608
PIK3CB	ENSG00000051382
PIK3R2	ENSG00000105647
PIK3R1	ENSG00000145675

```bash
gunzip -d eQTL_final1.txt.gz

 # which genes were tested? (all above not commented on)
for f1 in ENSG00000160293 ENSG00000147459 ENSG00000077264 ENSG00000075651 ENSG00000179604 ENSG00000070831 ENSG00000173327 ENSG00000133895 ENSG00000121879 ENSG00000171608 ENSG00000051382 ENSG00000105647 ENSG00000145675;
do
grep $f1 eQTL_final1.txt
done

for f1 in ENSG00000160293 ENSG00000147459 ENSG00000077264 ENSG00000075651 ENSG00000179604 ENSG00000070831 ENSG00000173327 ENSG00000133895 ENSG00000121879 ENSG00000171608 ENSG00000051382 ENSG00000105647 ENSG00000145675;
do
grep $f1 SigG_12_E_10_manyperm.txt
done
 # MLK3 is sig
ENSG00000173327.7 11 65615382 65615382 - 257 -7632 rs1151490 11 65623014 65623014 89 74.1029 1.06373 63.9991 2.95671e-10 0.872824 9.999e-05 2.27281e-07 2.49899529207714e-05
```

let's grab all of the ENSIDs from the eGenes and pass through string.

```bash
 #grabbing IDs
awk '{print $1}' SigG_12_E_10_manyperm.txt | cut -d '.' -f1 - > SigG_12_E_10_manyperm_ensIDs.txt
```

Some genes were not found by string:
ENSG00000263753
ENSG00000204110
ENSG00000227920
ENSG00000278705
ENSG00000280670
ENSG00000225783
ENSG00000262370
ENSG00000184274
ENSG00000275395
ENSG00000278129
ENSG00000268362
ENSG00000258504
ENSG00000196912
ENSG00000250312
ENSG00000260274
ENSG00000275342
ENSG00000237753
ENSG00000232229
ENSG00000282889
ENSG00000251136
ENSG00000258181
ENSG00000274245
ENSG00000273841
ENSG00000245954
ENSG00000259715
25 genes not found

ENSP00000353157 (SEPT2 - physically interacts with CDC42)


```bash
 # Pak1 (ENSG00000149269) shows up which interacts with cdc42
grep "ENSG00000149269" SigG_12_E_10_manyperm_ensIDs.txt

for f1 in ENSG00000149269;
do
grep $f1 SigG_12_E_10_manyperm.txt
done
 # ENSG00000149269.9 11 77474635 77474635 - 292 6982 rs2154754 11 77467653 77467653 89 76.7451 1.06402 101.331 3.71897e-09 0.773758 9.999e-05 1.98301e-06 0.000161672581847711
```

rs2154754 is the SNP name

Let's make some analysis plots!

ENSG00000172057 is ORMDL3, which is the number one gene associated with atopic asthma.

```bash
 # first I want to make a file with all of the information for the genes of interest
 # these are the genes that Deepa identified previously (22/24 are still present)
 #  (ENSG00000105793, ENSG00000106992) are no longer present
for f1 in ENSG00000064012 ENSG00000105793 ENSG00000106992 ENSG00000111450 ENSG00000113319 ENSG00000113658 ENSG00000121742 ENSG00000125245 ENSG00000149269 ENSG00000166263 ENSG00000169047 ENSG00000171033 ENSG00000171105 ENSG00000172057 ENSG00000173327 ENSG00000173890 ENSG00000176049 ENSG00000177426 ENSG00000179344 ENSG00000196126 ENSG00000196526 ENSG00000196735 ENSG00000198502 ENSG00000204138;
do
grep $f1 SigG_12_E_10_manyperm.txt >> SigG_12_E_10_manyperm_selected_old.txt
done

wc -l SigG_12_E_10_manyperm_selected_old.txt
 # 22
```

We want to check to see if genotype is contributing to the differential expression of any of the genes found within the original 89 genes identifeid via RNAseq. These genes can be found in the following file: `/home/greally-lab/Deepa_Andrew/04-2018_Helptagging_cellprop_analysis/04_2018_hcgenes_rnaseqpaper_fix.txt`

how many of these genes are also within the eGenes (genes associated with an eQTL)?

```bash
grep -F -f /home/greally-lab/Deepa_Andrew/04-2018_Helptagging_cellprop_analysis/04_2018_hcgenes_rnaseqpaper_fix.txt SigG_12_E_10_manyperm.txt > SigG_12_E_10_manyperm_hcgenes.txt
 # ENSG00000168404.12 (MLKL) and ENSG00000066923.17 (STAG3)
 # only 2 hc_genes found in eQTLs
 # let's include PAK1 too (ENSG00000149269)
```

Let's make some graph for the 2 hc genes found as eGenes

```bash
mkdir ../../Validate2
cd ../../Validate2
mkdir hc_genes
cd hc_genes
cp ../../Permutation2/perm_pc10_manyperm/SigG_12_E_10_manyperm.txt .
cp ../../Bams/quants2/Filtered_matrix_gene_rpkm.bed.gz .
gunzip -d Filtered_matrix_gene_rpkm.bed.gz
cp ../../Bams/PCA2/sample_metadata_eQTL.txt .
cp ../../Genotype2/SplX_converted_noHH_eQTL_snpfilt_hg38_sort.vcf.gz
gunzip -d SplX_converted_noHH_eQTL_snpfilt_hg38_sort.vcf.gz
cp  ../../Permutation2/perm_pc10_cond/eQTL_final3.txt .


module load R/3.4.0/gcc.4.7.4

 # PAK1
Rscript Plot_eQTL.R PAK1 ENSG00000149269 eQTL_final3.txt SplX_converted_noHH_eQTL_snpfilt_hg38_sort.vcf sample_metadata_eQTL.txt Filtered_matrix_gene_rpkm.bed

 # MLKL
Rscript Plot_eQTL.R MLKL ENSG00000168404 eQTL_final3.txt SplX_converted_noHH_eQTL_snpfilt_hg38_sort.vcf sample_metadata_eQTL.txt Filtered_matrix_gene_rpkm.bed

 # STAG3
Rscript Plot_eQTL.R STAG3 ENSG00000066923 eQTL_final3.txt SplX_converted_noHH_eQTL_snpfilt_hg38_sort.vcf sample_metadata_eQTL.txt Filtered_matrix_gene_rpkm.bed
```
for some reason I had to run the script manually instead of as a function because the pdfs were not being successfully generated on the loginnode

We wanted to look at some other genes
```
   hgnc_symbol ensembl_gene_id
1        ACACB ENSG00000076555
2         BDH2 ENSG00000164039
3         CBR3 ENSG00000159231
4         CCT3 ENSG00000163468
5         CCT7 ENSG00000135624
7        EPHA4 ENSG00000116106
8         FASN ENSG00000169710
9          GSR ENSG00000104687
10        INSR ENSG00000171105
11        IRS1 ENSG00000169047
12     MAP3K11 ENSG00000173327
13         ME3 ENSG00000151376
14        PASK ENSG00000115687
```

```bash
while read -r a b c; do
	qsub -S /bin/bash -N plot_${b}_eQTL -cwd -R y -l h_vmem=5.6G -j y << EOF
	module load R/3.4.0/gcc.4.7.4
    Rscript Plot_eQTL.R $b $c eQTL_final3.txt SplX_converted_noHH_eQTL_snpfilt_hg38_sort.vcf sample_metadata_eQTL.txt Filtered_matrix_gene_rpkm.bed
EOF
done < Genes_of_interest.txt
```

I also wanted to plot ORMDL3
```bash
qsub -S /bin/bash -N plot_ORMDL3_eQTL -cwd -R y -l h_vmem=5.6G -j y << EOF
module load R/3.4.0/gcc.4.7.4
Rscript Plot_eQTL.R ORMDL3 ENSG00000172057 eQTL_final3.txt SplX_converted_noHH_eQTL_snpfilt_hg38_sort.vcf sample_metadata_eQTL.txt Filtered_matrix_gene_rpkm.bed
EOF
```

I also wanted to plot SMAD1
```bash
qsub -S /bin/bash -N plot_SMAD1_eQTL -cwd -R y -l h_vmem=5.6G -j y << EOF
module load R/3.4.0/gcc.4.7.4
Rscript Plot_eQTL.R SMAD1 ENSG00000170365 eQTL_final3.txt SplX_converted_noHH_eQTL_snpfilt_hg38_sort.vcf sample_metadata_eQTL.txt Filtered_matrix_gene_rpkm.bed
EOF
```

  I also wanted to plot TG (ENSG00000042832) and RPS27L (ENSG00000185088)
```bash
qsub -S /bin/bash -N plot_TG_eQTL -cwd -R y -l h_vmem=5.6G -j y << EOF
module load R/3.4.0/gcc.4.7.4
Rscript Plot_eQTL.R TG ENSG00000042832 eQTL_final3.txt SplX_converted_noHH_eQTL_snpfilt_hg38_sort.vcf sample_metadata_eQTL.txt Filtered_matrix_gene_rpkm.bed
EOF

qsub -S /bin/bash -N plot_RPS_eQTL -cwd -R y -l h_vmem=5.6G -j y << EOF
module load R/3.4.0/gcc.4.7.4
Rscript Plot_eQTL.R RPS27L ENSG00000185088 eQTL_final3.txt SplX_converted_noHH_eQTL_snpfilt_hg38_sort.vcf sample_metadata_eQTL.txt Filtered_matrix_gene_rpkm.bed
EOF
```
There were several meQTL CpGs that were also identified as differentially methylated, associated with the following genes:
 ACER3 (ENSG00000078124) for promoter, PLA2G6 (ENSG00000184381) for gene body, and SEPT9 (ENSG00000184640) for cis-regulatory site

Are any of these genes also eQTLs? no.
```bash
 grep -E 'ENSG00000078124|ENSG00000184381|ENSG00000184640' eQTL_final3.txt
```


### Annotating the eQTLs

#### TSS Density

 `SigG_12_E_10_manyperm.txt` is our final eQTL list of 326 genes. Now to make a bed file of the QTLs and the TSSs of the significant eQTLs. I'll also take the ATACseq peaks to look at enrichment.  

Looking at the density of eQTLs around TSSs of genes

```bash
 # generate total and eGene beds
cat eQTL_final1.txt | awk '{ print $2, $3-1, $4, $1, $8, $5 }' | tr " " "\t" | sort -k1,1 -k2,2n > eQTL_final1_allGenes.bed
cat SigG_12_E_10_manyperm.txt | awk '{ print $2, $3-1, $4, $1, $8, $5 }' | tr " " "\t" | sort -k1,1 -k2,2n > SigG_12_E_10_manyperm_eGenes.bed


cd ../perm_pc10_cond
cat eQTL_final3.txt | awk '{ print $9, $10-1, $11, $8, $1, $5 }' | tr " " "\t" | sort -k1,1 -k2,2n > eQTL_final3_eQTL.bed

cd ../../Validate2
mkdir TSS_density
cd TSS_density/
cp ../../Permutation2/perm_pc10_cond/eQTL_final3_eQTL.bed .
cp ../../Permutation2/perm_pc10_manyperm/eQTL_final1_allGenes.bed .
cp ../../Permutation2/perm_pc10_manyperm/SigG_12_E_10_manyperm_eGenes.bed .


 # grab the TSSs
 # generate bed file with +/- 100kb from TSS
awk 'OFS="\t" {if ($6=="+") {print $1,$2-100000,$2+100000,$4,$5,$6} if ($6=="-") {print $1,$3-100000,$3+100000,$4,$5,$6}}' eQTL_final1_allGenes.bed | awk 'OFS="\t" {if ($2 < 0) {print $1,1,$3,$4,$5,$6} else {print $0}}' > eQTL_final1_allGenes_100kbTSS.bed

 # grabbing the TSSs from the genes associated with variants
awk 'OFS="\t" {if ($6=="+") {print $1,$2-100000,$2+100000,$4,$5,$6} if ($6=="-") {print $1,$3-100000,$3+100000,$4,$5,$6}}' SigG_12_E_10_manyperm_eGenes.bed | awk 'OFS="\t" {if ($2 < 0) {print $1,1,$3,$4,$5,$6} else {print $0}}' > SigG_12_E_10_manyperm_eGenes_100kbTSS.bed

 # seperate minus and plus stranded genes
awk '$6=="+"' eQTL_final1_allGenes_100kbTSS.bed >  eQTL_final1_allGenes_100kbTSS_plus.bed
awk '$6=="-"' eQTL_final1_allGenes_100kbTSS.bed >  eQTL_final1_allGenes_100kbTSS_minus.bed

awk '$6=="+"' SigG_12_E_10_manyperm_eGenes_100kbTSS.bed >  SigG_12_E_10_manyperm_eGenes_100kbTSS_plus.bed
awk '$6=="-"' SigG_12_E_10_manyperm_eGenes_100kbTSS.bed >  SigG_12_E_10_manyperm_eGenes_100kbTSS_minus.bed

 # creating 1kb windows for the minus and plus
module load bedtools2
bedtools makewindows -b eQTL_final1_allGenes_100kbTSS_plus.bed -w 1000 -i winnum > eQTL_final1_allGenes_100kbTSS_plus_1kbwin.bed
bedtools makewindows -b eQTL_final1_allGenes_100kbTSS_minus.bed -w 1000 -i winnum > eQTL_final1_allGenes_100kbTSS_minus_1kbwin.bed
bedtools makewindows -b SigG_12_E_10_manyperm_eGenes_100kbTSS_plus.bed -w 1000 -i winnum > SigG_12_E_10_manyperm_eGenes_100kbTSS_plus_1kbwin.bed
bedtools makewindows -b SigG_12_E_10_manyperm_eGenes_100kbTSS_minus.bed -w 1000 -i winnum > SigG_12_E_10_manyperm_eGenes_100kbTSS_minus_1kbwin.bed
```

Now I need reverse the minus bins to reflect the directionality of the gene (TSS)

```R
library(data.table)
 # all tested genes: flip the minus window numbers and combine the minus and plus
allG_win_plus <- fread("eQTL_final1_allGenes_100kbTSS_plus_1kbwin.bed")
allG_win_mins <- fread("eQTL_final1_allGenes_100kbTSS_minus_1kbwin.bed")

allG_win_mins$V4 <- (allG_win_mins$V4-201)*(-1)
allG_win_comb <- rbind(allG_win_plus,allG_win_mins)

write.table(allG_win_comb,"eQTL_final1_allGenes_100kbTSS_1kbwin.bed",col.names=F,row.names=F,quote=F,sep="\t",append=F)

 # eGenes: flip the minus window numbers and combine the minus and plus
eGene_win_plus <- fread("SigG_12_E_10_manyperm_eGenes_100kbTSS_plus_1kbwin.bed")
eGene_win_mins <- fread("SigG_12_E_10_manyperm_eGenes_100kbTSS_minus_1kbwin.bed")

eGene_win_mins$V4 <- (eGene_win_mins$V4-201)*(-1)
eGene_win_comb <- rbind(eGene_win_plus,eGene_win_mins)

write.table(eGene_win_comb,"SigG_12_E_10_manyperm_eGenes_100kbTSS_1kbwin.bed",col.names=F,row.names=F,quote=F,sep="\t",append=F)
```

Now intersect the bins with the all tested variants / eQTL variants

```bash
 # intersect all tested variants with all TSSs
 # need to make a bed file of all of the variants
cp ../hc_genes/SplX_converted_noHH_eQTL_snpfilt_hg38_sort.vcf .
awk 'NR>28 {OFS="\t"; print $1,$2,$2+1,$3}' SplX_converted_noHH_eQTL_snpfilt_hg38_sort.vcf > SplX_converted_noHH_eQTL_snpfilt_hg38_sort_allVar.bed
bedtools intersect -a eQTL_final1_allGenes_100kbTSS_1kbwin.bed -b SplX_converted_noHH_eQTL_snpfilt_hg38_sort_allVar.bed -wo > allVar_allGene_100kbTSS_1kbwin.txt

 # note that some variants are not within 100kb of TSSs
 # OR different chromosomes
wc -l SplX_converted_noHH_eQTL_snpfilt_hg38_sort_allVar.bed
 # 497036 SplX_converted_noHH_eQTL_snpfilt_hg38_sort_allVar.bed
wc -l allVar_allGene_100kbTSS_1kbwin.txt
 # 429183 allVar_allGene_100kbTSS_1kbwin.txt

 # intersect the eQTLs with sig genes
bedtools intersect -a SigG_12_E_10_manyperm_eGenes_100kbTSS_1kbwin.bed -b eQTL_final3_eQTL.bed -wo > eQTL_final3_eQTL_eGene_100kbTSS_1kbwin.txt
```
see TSS_plots.R and the resulting plots [TSS_allvar_allG.pdf](TSS_allvar_allG.pdf) and
[TSS_eQTL_eGenes.pdf](TSS_eQTL_eGenes.pdf)


Using QTLtools to get the density of TSSs around eQTLs...
```bash
 #all genes
qsub -S /bin/bash -N eQTL_TSS_dens -cwd -R y -l h_vmem=5.6G -j y << EOF
module load QTLtools/1.1/gcc.4.7.4
QTLtools fdensity --qtl eQTL_final3_eQTL.bed --bed eQTL_final1_allGenes.bed --out dens_TSS_eQTL.txt
EOF
 # eGenes
qsub -S /bin/bash -N eQTL_TSS_dens -cwd -R y -l h_vmem=5.6G -j y << EOF
module load QTLtools/1.1/gcc.4.7.4
QTLtools fdensity --qtl eQTL_final3_eQTL.bed --bed SigG_12_E_10_manyperm_eGenes.bed --out dens_TSS_eGene_eQTL.txt
EOF
```

Now plot the density using R
```R
 # What about density from QTL tools - all TSSs
D <- read.table("dens_TSS_eQTL.txt", head=FALSE, stringsAsFactors=FALSE)

pdf("allTSS_dens_eQTL.pdf", width=6, height = 4, family="ArialMT")
plot((D$V1+D$V2)/2, D$V3, type="l", xlab="Distance to QTLs", ylab="#annotations/kb",
     main="all gene TSSs around eQTLs")
dev.off()

 # What about density from QTL tools - only eGene TSSs
D_Egene <- read.table("dens_TSS_eGene_eQTL.txt", head=FALSE, stringsAsFactors=FALSE)

pdf("eGeneTSS_dens_eQTL.pdf", width=6, height = 4, family="ArialMT")
plot((D_Egene$V1+D_Egene$V2)/2, D_Egene$V3, type="l", xlab="Distance to QTLs",
     ylab="#annotations/kb", main = "Egene TSSs around eQTLs")
dev.off()
```

see the plots `eGeneTSS_dens_eQTL.pdf` and `allTSS_dens_eQTL.pdf`. both show clear enrichment of eQTLs around TSSs.

#### enhancer enrichment

Let's see if QTLs are enriched in enhancer annotations.

I have enhancer enrichments in `ATAC_ens_enh_hg38_no10kb_merged.bed` that were made for annotating the CpGs in differential methylation analysis. See `Compare_annotations.html` and `Annotate_enhancers.html` for information on the enahncer annotations.

```bash
cd /home/greally-lab/Deepa_Andrew/Deepa-helptagging/CpG_Annotations
mkdir ../../eQTL/Validate2/enh_enrich/
cp ATAC_ens_enh_hg38_no10kb_merged.bed ../../eQTL/Validate2/enh_enrich/
cd ../../eQTL/Validate2/enh_enrich/
cp ../TSS_density/eQTL_final3_eQTL.bed .

wc -l ATAC_ens_enh_hg38_no10kb_merged.bed # 65510

 # removing chr prefix on chromosomes
cat ATAC_ens_enh_hg38_no10kb_merged.bed | sed 's/^chr//' > ATAC_ens_enh_hg38_no10kb_merged_nochr.bed

qsub -S /bin/bash -N eQTL_enh_enrich -cwd -R y -l h_vmem=5.6G -j y << EOF
module load QTLtools/1.1/gcc.4.7.4
QTLtools fenrich --qtl eQTL_final3_eQTL.bed --tss eQTL_final1_allGenes.bed --bed ATAC_ens_enh_hg38_no10kb_merged_nochr.bed --permute 1000 --out enrich_enh_eQTL.txt
EOF
```

Here are the results. Enhancers are enriched
```
Enrichment analysis
  * #observed overlaps = 72 / 1574 (4.57%)
  * #null overlaps = 40.21 +/- 6.32 (2.55% +/- 0.40%)
  * empirical p-value = 0.000999001
  * Odd ratio = 1.8383 [2.6468,1.3757]
```

#### other enrichments (Ensembl Anno, ATAC)


First I need to get all of the other enrichments. I can take the various annotations in the ensembl regulatory build.

downloaded the the regulatory build from ftp://ftp.ensembl.org/pub/release-90/regulation/homo_sapiens/RegulatoryFeatureActivity/. I got the CD4 TCell Venous blood build.

```bash
cp ../../../Deepa-helptagging/CpG_Annotations/homo_sapiens.GRCh38.CD4_ab_T_cell_VB.Regulatory_Build.regulatory_activity.20161111.gff .

 ### what are the different annotations in the build and how many?
awk '{print $3}' homo_sapiens.GRCh38.CD4_ab_T_cell_VB.Regulatory_Build.regulatory_activity.20161111.gff | sort | uniq -c
 ####  92171 CTCF_binding_site
 ####  36452 enhancer
 #### 129294 open_chromatin_region
 ####  16384 promoter
 ####  45706 promoter_flanking_region
 ####  21922 TF_binding_site

 # grabbing all of the different annotations
awk '$3 == "enhancer" {print $1,$4,$5,"ens_enh_",".","."}' F="\t" OFS="\t" homo_sapiens.GRCh38.CD4_ab_T_cell_VB.Regulatory_Build.regulatory_activity.20161111.gff | sort -k 1,1 -k2,2n | awk '{print $1,$2,$3,$4NR,$5,$6}' OFS="\t" > ens_enhancer_hg38.bed

awk '$3 == "TF_binding_site" {print $1,$4,$5,"ens_enh_",".","."}' F="\t" OFS="\t" homo_sapiens.GRCh38.CD4_ab_T_cell_VB.Regulatory_Build.regulatory_activity.20161111.gff | sort -k 1,1 -k2,2n | awk '{print $1,$2,$3,$4NR,$5,$6}' OFS="\t" > ens_TF_binding_hg38.bed

awk '$3 == "CTCF_binding_site" {print $1,$4,$5,"ens_enh_",".","."}' F="\t" OFS="\t" homo_sapiens.GRCh38.CD4_ab_T_cell_VB.Regulatory_Build.regulatory_activity.20161111.gff | sort -k 1,1 -k2,2n | awk '{print $1,$2,$3,$4NR,$5,$6}' OFS="\t" > ens_CTCF_binding_site_hg38.bed

awk '$3 == "promoter" {print $1,$4,$5,"ens_enh_",".","."}' F="\t" OFS="\t" homo_sapiens.GRCh38.CD4_ab_T_cell_VB.Regulatory_Build.regulatory_activity.20161111.gff | sort -k 1,1 -k2,2n | awk '{print $1,$2,$3,$4NR,$5,$6}' OFS="\t" > ens_promoter_hg38.bed

awk '$3 == "promoter_flanking_region" {print $1,$4,$5,"ens_enh_",".","."}' F="\t" OFS="\t" homo_sapiens.GRCh38.CD4_ab_T_cell_VB.Regulatory_Build.regulatory_activity.20161111.gff | sort -k 1,1 -k2,2n | awk '{print $1,$2,$3,$4NR,$5,$6}' OFS="\t" > ens_promoter_flanking_region_hg38.bed

awk '$3 == "open_chromatin_region" {print $1,$4,$5,"ens_enh_",".","."}' F="\t" OFS="\t" homo_sapiens.GRCh38.CD4_ab_T_cell_VB.Regulatory_Build.regulatory_activity.20161111.gff | sort -k 1,1 -k2,2n | awk '{print $1,$2,$3,$4NR,$5,$6}' OFS="\t" > ens_open_chromatin_region_hg38.bed

qsub -S /bin/bash -N eQTL_ens_enh_enrich -cwd -R y -l h_vmem=5.6G -j y << EOF
module load QTLtools/1.1/gcc.4.7.4
QTLtools fenrich --qtl eQTL_final3_eQTL.bed --tss eQTL_final1_allGenes.bed --bed ens_enhancer_hg38.bed --permute 1000 --out enrich_ens_enh_eQTL.txt
EOF

qsub -S /bin/bash -N eQTL_ens_TF_enrich -cwd -R y -l h_vmem=5.6G -j y << EOF
module load QTLtools/1.1/gcc.4.7.4
QTLtools fenrich --qtl eQTL_final3_eQTL.bed --tss eQTL_final1_allGenes.bed --bed ens_TF_binding_hg38.bed --permute 1000 --out enrich_ens_TF_eQTL.txt
EOF

qsub -S /bin/bash -N eQTL_ens_CTCF_enrich -cwd -R y -l h_vmem=5.6G -j y << EOF
module load QTLtools/1.1/gcc.4.7.4
QTLtools fenrich --qtl eQTL_final3_eQTL.bed --tss eQTL_final1_allGenes.bed --bed ens_CTCF_binding_site_hg38.bed --permute 1000 --out enrich_ens_CTCF_eQTL.txt
EOF

qsub -S /bin/bash -N eQTL_ens_prom_enrich -cwd -R y -l h_vmem=5.6G -j y << EOF
module load QTLtools/1.1/gcc.4.7.4
QTLtools fenrich --qtl eQTL_final3_eQTL.bed --tss eQTL_final1_allGenes.bed --bed ens_promoter_hg38.bed --permute 1000 --out enrich_ens_prom_eQTL.txt
EOF

qsub -S /bin/bash -N eQTL_ens_open_chrom_enrich -cwd -R y -l h_vmem=5.6G -j y << EOF
module load QTLtools/1.1/gcc.4.7.4
QTLtools fenrich --qtl eQTL_final3_eQTL.bed --tss eQTL_final1_allGenes.bed --bed ens_open_chromatin_region_hg38.bed --permute 1000 --out enrich_ens_open_chrom_eQTL.txt
EOF

qsub -S /bin/bash -N eQTL_ens_prom_flank_enrich -cwd -R y -l h_vmem=5.6G -j y << EOF
module load QTLtools/1.1/gcc.4.7.4
QTLtools fenrich --qtl eQTL_final3_eQTL.bed --tss eQTL_final1_allGenes.bed --bed ens_promoter_flanking_region_hg38.bed --permute 1000 --out enrich_ens_prom_flank_chrom_eQTL.txt
EOF

```


 enhancer Odds:
```
  * #observed overlaps = 17 / 1574 (1.08%)
  * #null overlaps = 18.22 +/- 4.21 (1.16% +/- 0.27%)
  * empirical p-value = 0.868132
  * Odd ratio = 0.9438 [1.7076,0.6256]
```

 TF binding Odds:
```
  * #observed overlaps = 8 / 1574 (0.51%)
  * #null overlaps = 12.41 +/- 3.42 (0.79% +/- 0.22%)
  * empirical p-value = 0.230769
  * Odd ratio = 0.6650 [1.3350,0.3969]
```

 promoter Odds:
```
  * #observed overlaps = 92 / 1574 (5.84%)
  * #null overlaps = 85.83 +/- 7.25 (5.45% +/- 0.46%)
  * empirical p-value = 0.422577
  * Odd ratio = 1.0741 [1.2950,0.9150]
```

 promoter flank Odds:
```
  * #observed overlaps = 79 / 1574 (5.02%)
  * #null overlaps = 72.24 +/- 8.27 (4.59% +/- 0.53%)
  * empirical p-value = 0.438561
  * Odd ratio = 1.1024 [1.4324,0.8923]
```

 open chrom Odds:
```
  * #observed overlaps = 46 / 1574 (2.92%)
  * #null overlaps = 40.14 +/- 6.44 (2.55% +/- 0.41%)
  * empirical p-value = 0.410589
  * Odd ratio = 1.1545 [1.6622,0.8639]
```

 CTCF Odds:
```
  * #observed overlaps = 34 / 1574 (2.16%)
  * #null overlaps = 36.38 +/- 5.99 (2.31% +/- 0.38%)
  * empirical p-value = 0.782218
  * Odd ratio = 0.9432 [1.3679,0.7019]
```

Let's also look at ATAC peaks (not within promoter)
```
cp ../../../Deepa-helptagging/CpG_Annotations/CD4_ATAC_peaks_hg38_noProm_sort.bed .
wc -l CD4_ATAC_peaks_hg38_noProm_sort.bed
 # 31709 CD4_ATAC_peaks_hg38_noProm_sort.bed

 # need to remove chr
cat CD4_ATAC_peaks_hg38_noProm_sort.bed | sed 's/^chr//' > CD4_ATAC_peaks_hg38_noProm_sort_nochr.bed

qsub -S /bin/bash -N eQTL_ATAC_enrich -cwd -R y -l h_vmem=5.6G -j y -q highmem.q << EOF
module load QTLtools/1.1/gcc.4.7.4
QTLtools fenrich --qtl eQTL_final3_eQTL.bed --tss eQTL_final1_allGenes.bed --bed CD4_ATAC_peaks_hg38_noProm_sort_nochr.bed --permute 1000 --out enrich_ATAC_eQTL.txt
EOF

 A very large enrichment in ATAC peaks
# * observed overlaps = 59 / 1574 (3.75%)
# * null overlaps = 22.93 +/- 4.82 (1.46% +/- 0.31%)
# * empirical p-value = 0.000999001
# * Odd ratio = 2.6262 [4.3395,1.8766]
```

#### IDEAS enrichment

I wanted to have a Segway/ChromHMM track to look for enrichment. I stumbled across the roadmap consortiums "chromHMM" paper where they use a slightly modified algorithm IDEAS to annotate the genome. They show that their algorithm better defines active annotations that are enriched in eQTLs.

1. Zhang, Y., & Hardison, R. C. (2017). Accurate and reproducible functional maps in 127 human cell types via 2D genome segmentation. Nucleic Acids Research, 45(17), 9823–9836. http://doi.org/10.1093/nar/gkx659 (IDEAS and eQTLs)
2. Zhang, Y., An, L., Yue, F., & Hardison, R. C. (2016). Jointly characterizing epigenetic dynamics across multiple human cell types. Nucleic Acids Research, 44(14), 6721–6731. http://doi.org/10.1093/nar/gkw278 (Initial IDEAS paper)


I went to the penn state genome browser: http://main.genome-browser.bx.psu.edu/cgi-bin/hgTables?hgsid=322193_XWMysfws5Hx8lFaToIhaCY8RfAcl from http://personal.psu.edu/yzz2/IDEAS/

I downloaded the IDEAS track for E038_BLD.CD4.NPC (E038_Primary T helper naive cells from peripheral blood; `CD4_IDEAS.gtf`), E039_BLD.CD4.CD25M.CD45RA.NPC ( E039_Primary T helper naive cells from peripheral blood; `CD4_naive_IDEAS.gtf`), EO43_BLD.CD4.CD25M.TPC (E043_Primary T helper cells from peripheral blood; `CD4_TPC_IDEAS.gtf`)


Let's look at the various states
```bash
awk '{print $10}' CD4_IDEAS.gtf | tr -d '";' | sort | uniq -c
```
```
 704660 0_Quies
  12380 10_TssA
  20197 11_EnhBiv
  26284 12_Het/ReprPC
   6864 13_ReprPC
  11535 14_TssWk
   4635 15_TssBiv
   8821 16_TxRepr
  25667 17_EnhGA (must be active? genic enhancer, they feel stonger effect size)
   6146 18_Enh/Het
   1966 19_Enh/ReprPC
 253055 1_ReprPCWk
 157923 2_TxWk
 259725 3_HetWk
 106153 4_Enh
  58358 5_Tx
  53047 6_EnhG (Genic Enhancers, within gene)
  15641 7_ZNF/Rpts
  29799 8_TssAFlnk
   6473 9_Het
```

Let's segment and find the enrichment for the various states

First up is CD4_IDEAS.gtf

```bash

 # need to remove chr
cat CD4_IDEAS.gtf | sed 's/^chr//' > CD4_IDEAS_nochr.gtf

 # now to seperate by states
mkdir CD4_IDEAS_states

awk '{print $10}' CD4_IDEAS.gtf | tr -d '";' | sort -u > IDEAS_states.txt

 # need to edit the /'s in names in IDEAS_states.txt

while read state;
do
echo "$state"
awk -v i="$state" '$10~i {OFS="\t"; print $1,$4,$5,$12,$6,$7}' CD4_IDEAS_nochr.gtf | tr -d '";' > CD4_IDEAS_states/CD4_IDEAS_nochr_${state}.bed
done < IDEAS_states.txt

while read state;
do
echo "$state"
qsub -S /bin/bash -N eQTL_CD4_${state}_enrich -cwd -R y -l h_vmem=5.6G -j y -q highmem.q << EOF
module load QTLtools/1.1/gcc.4.7.4
QTLtools fenrich --qtl eQTL_final3_eQTL.bed --tss eQTL_final1_allGenes.bed --bed CD4_IDEAS_states/CD4_IDEAS_nochr_${state}.bed --permute 1000 --out CD4_IDEAS_states/enrich_${state}_eQTL.txt
EOF
done < IDEAS_states.txt
```
results:
```
0 Quies
Enrichment analysis
  * #observed overlaps = 743 / 1574 (47.20%)
  * #null overlaps = 774.62 +/- 19.48 (49.21% +/- 1.24%)
  * empirical p-value = 0.104895
  * Odd ratio = 0.9218 [1.0206,0.8390]

1 ReprPCwk (?)
Enrichment analysis
  * #observed overlaps = 124 / 1574 (7.88%)
  * #null overlaps = 144.27 +/- 11.52 (9.17% +/- 0.73%)
  * empirical p-value = 0.0889111
  * Odd ratio = 0.8492 [1.0088,0.7157]

2 Tx Weak
Enrichment analysis
  * #observed overlaps = 175 / 1574 (11.12%)
  * #null overlaps = 187.53 +/- 12.37 (11.91% +/- 0.79%)
  * empirical p-value = 0.334665
  * Odd ratio = 0.9222 [1.0755,0.7993]

3 Het Weak
Enrichment analysis
  * #observed overlaps = 50 / 1574 (3.18%)
  * #null overlaps = 42.66 +/- 6.27 (2.71% +/- 0.40%)
  * empirical p-value = 0.260739
  * Odd ratio = 1.1681 [1.6330,0.8893]

4 Enh
Enrichment analysis
  * #observed overlaps = 108 / 1574 (6.86%)
  * #null overlaps = 85.21 +/- 8.79 (5.41% +/- 0.56%)
  * empirical p-value = 0.014985
  * Odd ratio = 1.2905 [1.6316,1.0632]

5 Tx
Enrichment analysis
  * #observed overlaps = 123 / 1574 (7.81%)
  * #null overlaps = 106.37 +/- 9.48 (6.76% +/- 0.60%)
  * empirical p-value = 0.0949051
  * Odd ratio = 1.1740 [1.4314,0.9826]

6 EnhG
Enrichment analysis
  * #observed overlaps = 44 / 1574 (2.80%)
  * #null overlaps = 50.56 +/- 6.80 (3.21% +/- 0.43%)
  * empirical p-value = 0.374625
  * Odd ratio = 0.8765 [1.1624,0.6785]

7 ZNF/Rpts
Enrichment analysis
  * #observed overlaps = 32 / 1574 (2.03%)
  * #null overlaps = 12.42 +/- 3.50 (0.79% +/- 0.22%)
  * empirical p-value = 0.000999001
  * Odd ratio = 2.7013 [5.4233,1.6125]

8 TSSA flank
Enrichment analysis
  * #observed overlaps = 25 / 1574 (1.59%)
  * #null overlaps = 25.88 +/- 5.01 (1.64% +/- 0.32%)
  * empirical p-value = 1
  * Odd ratio = 1.0000 [1.5716,0.6895]

9 Het
Enrichment analysis
  * #observed overlaps = 6 / 1574 (0.38%)
  * #null overlaps = 3.09 +/- 1.74 (0.20% +/- 0.11%)
  * empirical p-value = 0.190809
  * Odd ratio = 2.0038 [inf,0.8566]

10 TSS
Enrichment analysis
  * #observed overlaps = 15 / 1574 (0.95%)
  * #null overlaps = 15.79 +/- 4.07 (1.00% +/- 0.26%)
  * empirical p-value = 0.996004
  * Odd ratio = 0.9369 [1.6731,0.6214]

11 Enh Biv
Enrichment analysis
  * #observed overlaps = 12 / 1574 (0.76%)
  * #null overlaps = 11.24 +/- 3.50 (0.71% +/- 0.22%)
  * empirical p-value = 0.886114
  * Odd ratio = 1.0916 [2.4108,0.6287]

12 Het/ReprPC
Enrichment analysis
  * #observed overlaps = 1 / 1574 (0.06%)
  * #null overlaps = 6.08 +/- 2.51 (0.39% +/- 0.16%)
  * empirical p-value = 0.040959
  * Odd ratio = 0.1661 [0.4997,0.0903]

13 RepR
Enrichment analysis
  * #observed overlaps = 1 / 1574 (0.06%)
  * #null overlaps = 3.69 +/- 1.93 (0.23% +/- 0.12%)
  * empirical p-value = 0.246753
  * Odd ratio = 0.2495 [inf,0.1244]

14 TSS weak
Enrichment analysis
  * #observed overlaps = 2 / 1574 (0.13%)
  * #null overlaps = 5.12 +/- 2.22 (0.33% +/- 0.14%)
  * empirical p-value = 0.234765
  * Odd ratio = 0.3992 [2.0013,0.1990]

15 TSS Biv
Enrichment analysis
  * #observed overlaps = 0 / 1574 (0.00%)
  * #null overlaps = 2.93 +/- 1.67 (0.19% +/- 0.11%)
  * empirical p-value = 0.104895
  * Odd ratio = 0.0000 [-nan,0.0000]

16 TxRepR
Enrichment analysis
  * #observed overlaps = 16 / 1574 (1.02%)
  * #null overlaps = 7.71 +/- 2.81 (0.49% +/- 0.18%)
  * empirical p-value = 0.016983
  * Odd ratio = 2.0103 [5.3778,1.2331]

17 EnhGA
Enrichment analysis
  * #observed overlaps = 47 / 1574 (2.99%)
  * #null overlaps = 29.66 +/- 5.34 (1.88% +/- 0.34%)
  * empirical p-value = 0.004995
  * Odd ratio = 1.5841 [2.3916,1.1508]

18 Enh/Het
Enrichment analysis
  * #observed overlaps = 1 / 1574 (0.06%)
  * #null overlaps = 2.77 +/- 1.66 (0.18% +/- 0.11%)
  * empirical p-value = 0.480519
  * Odd ratio = 0.3329 [inf,0.1661]

19 Enh/ReprPC (PC = polycomb)
Enrichment analysis
  * #observed overlaps = 0 / 1574 (0.00%)
  * #null overlaps = 0.56 +/- 0.76 (0.04% +/- 0.05%)
  * empirical p-value = 1
  * Odd ratio = -nan [-nan,0.0000]
```

eQTLs were enriched in the following:
4 Enh
7 ZNF/Rpts
16 TxRepR
17 EnhGA


Next is CD4_naive_IDEAS.gtf

```bash

 # need to remove chr
cat CD4_naive_IDEAS.gtf | sed 's/^chr//' > CD4_naive_IDEAS_nochr.gtf

 # now to seperate by states
mkdir CD4_naive_IDEAS

while read state;
do
echo "$state"
awk -v i="$state" '$10~i {OFS="\t"; print $1,$4,$5,$12,$6,$7}' CD4_naive_IDEAS_nochr.gtf | tr -d '";' > CD4_naive_IDEAS/CD4_naive_IDEAS_nochr_${state}.bed
qsub -S /bin/bash -N eQTL_CD4_naive_${state}_enrich -cwd -R y -l h_vmem=5.6G -j y -q highmem.q << EOF
module load QTLtools/1.1/gcc.4.7.4
QTLtools fenrich --qtl eQTL_final3_eQTL.bed --tss eQTL_final1_allGenes.bed --bed CD4_naive_IDEAS/CD4_naive_IDEAS_nochr_${state}.bed --permute 1000 --out CD4_naive_IDEAS/enrich_${state}_eQTL.txt
EOF
done < IDEAS_states.txt
```

results:

```
0 Quies
Enrichment analysis
  * #observed overlaps = 667 / 1574 (42.38%)
  * #null overlaps = 670.46 +/- 19.40 (42.60% +/- 1.23%)
  * empirical p-value = 0.856144
  * Odd ratio = 0.9897 [1.0903,0.8972]

1 ReprPCwk (?)
Enrichment analysis
  * #observed overlaps = 154 / 1574 (9.78%)
  * #null overlaps = 192.57 +/- 13.18 (12.23% +/- 0.84%)
  * empirical p-value = 0.002997
  * Odd ratio = 0.7760 [0.9137,0.6710]

2 Tx Weak
Enrichment analysis
  * #observed overlaps = 149 / 1574 (9.47%)
  * #null overlaps = 165.60 +/- 12.46 (10.52% +/- 0.79%)
  * empirical p-value = 0.206793
  * Odd ratio = 0.8869 [1.0545,0.7616]

3 Het Weak
Enrichment analysis
  * #observed overlaps = 63 / 1574 (4.00%)
  * #null overlaps = 56.32 +/- 7.45 (3.58% +/- 0.47%)
  * empirical p-value = 0.406593
  * Odd ratio = 1.1302 [1.5208,0.8698]

4 Enh
Enrichment analysis
  * #observed overlaps = 121 / 1574 (7.69%)
  * #null overlaps = 85.32 +/- 8.61 (5.42% +/- 0.55%)
  * empirical p-value = 0.000999001
  * Odd ratio = 1.4588 [1.8443,1.2018]

5 Tx
Enrichment analysis
  * #observed overlaps = 161 / 1574 (10.23%)
  * #null overlaps = 142.92 +/- 10.91 (9.08% +/- 0.69%)
  * empirical p-value = 0.104895
  * Odd ratio = 1.1402 [1.3561,0.9730]

6 EnhG
Enrichment analysis
  * #observed overlaps = 43 / 1574 (2.73%)
  * #null overlaps = 56.65 +/- 7.37 (3.60% +/- 0.47%)
  * empirical p-value = 0.0609391
  * Odd ratio = 0.7475 [1.0000,0.5946]

7 ZNF/Rpts
Enrichment analysis
  * #observed overlaps = 34 / 1574 (2.16%)
  * #null overlaps = 16.62 +/- 4.15 (1.06% +/- 0.26%)
  * empirical p-value = 0.000999001
  * Odd ratio = 2.0221 [3.8391,1.3679]

8 TSSA flank
Enrichment analysis
  * #observed overlaps = 25 / 1574 (1.59%)
  * #null overlaps = 28.40 +/- 5.22 (1.80% +/- 0.33%)
  * empirical p-value = 0.594406
  * Odd ratio = 0.8911 [1.3209,0.6352]

9 Het
Enrichment analysis
  * #observed overlaps = 13 / 1574 (0.83%)
  * #null overlaps = 6.75 +/- 2.60 (0.43% +/- 0.17%)
  * empirical p-value = 0.048951
  * Odd ratio = 1.8643 [6.5458,1.0840]

10 TSS
Enrichment analysis
  * #observed overlaps = 17 / 1574 (1.08%)
  * #null overlaps = 17.66 +/- 4.31 (1.12% +/- 0.27%)
  * empirical p-value = 0.998002
  * Odd ratio = 0.9438 [1.7076,0.6256]

11 Enh Biv
Enrichment analysis
  * #observed overlaps = 11 / 1574 (0.70%)
  * #null overlaps = 9.65 +/- 3.23 (0.61% +/- 0.20%)
  * empirical p-value = 0.734266
  * Odd ratio = 1.2238 [2.7623,0.6446]

12 Het/ReprPC
Enrichment analysis
  * #observed overlaps = 5 / 1574 (0.32%)
  * #null overlaps = 10.75 +/- 3.37 (0.68% +/- 0.21%)
  * empirical p-value = 0.112887
  * Odd ratio = 0.4528 [1.2508,0.2755]

13 RepR
Enrichment analysis
  * #observed overlaps = 1 / 1574 (0.06%)
  * #null overlaps = 8.00 +/- 2.82 (0.51% +/- 0.18%)
  * empirical p-value = 0.000999001
  * Odd ratio = 0.1244 [0.3329,0.0708]

14 TSS weak
Enrichment analysis
  * #observed overlaps = 3 / 1574 (0.19%)
  * #null overlaps = 3.15 +/- 1.70 (0.20% +/- 0.11%)
  * empirical p-value = 1
  * Odd ratio = 1.0000 [inf,0.4275]

15 TSS Biv
Enrichment analysis
  * #observed overlaps = 4 / 1574 (0.25%)
  * #null overlaps = 6.02 +/- 2.51 (0.38% +/- 0.16%)
  * empirical p-value = 0.596404
  * Odd ratio = 0.6658 [2.0025,0.3620]

16 TxRepR
Enrichment analysis
  * #observed overlaps = 1 / 1574 (0.06%)
  * #null overlaps = 2.44 +/- 1.58 (0.15% +/- 0.10%)
  * empirical p-value = 0.6004
  * Odd ratio = 0.4997 [inf,0.1661]

17 EnhGA
Enrichment analysis
  * #observed overlaps = 52 / 1574 (3.30%)
  * #null overlaps = 34.96 +/- 6.11 (2.22% +/- 0.39%)
  * empirical p-value = 0.010989
  * Odd ratio = 1.5023 [2.2065,1.1100]

18 Enh/Het
Enrichment analysis
  * #observed overlaps = 1 / 1574 (0.06%)
  * #null overlaps = 3.48 +/- 1.94 (0.22% +/- 0.12%)
  * empirical p-value = 0.298701
  * Odd ratio = 0.3329 [inf,0.1244]

19 Enh/ReprPC
Enrichment analysis
  * #observed overlaps = 0 / 1574 (0.00%)
  * #null overlaps = 0.43 +/- 0.66 (0.03% +/- 0.04%)
  * empirical p-value = 1
  * Odd ratio = -nan [-nan,0.0000]
```

eQTLs were enriched in the following:
4 Enh
7 ZNF/Rpts
9 Het is sig (but fairly low in terms of numbers only 13 eQTLs)
17 EnhGA

Last is CD4_TPC_IDEAS.gtf

```bash

 # need to remove chr
cat CD4_TPC_IDEAS.gtf | sed 's/^chr//' > CD4_TPC_IDEAS_nochr.gtf

 # now to seperate by states
mkdir CD4_TPC_IDEAS

while read state;
do
echo "$state"
awk -v i="$state" '$10~i {OFS="\t"; print $1,$4,$5,$12,$6,$7}' CD4_TPC_IDEAS_nochr.gtf | tr -d '";' > CD4_TPC_IDEAS/CD4_TPC_IDEAS_nochr_${state}.bed
qsub -S /bin/bash -N eQTL_CD4_TPC_${state}_enrich -cwd -R y -l h_vmem=5.6G -j y -q highmem.q << EOF
module load QTLtools/1.1/gcc.4.7.4
QTLtools fenrich --qtl eQTL_final3_eQTL.bed --tss eQTL_final1_allGenes.bed --bed CD4_TPC_IDEAS/CD4_TPC_IDEAS_nochr_${state}.bed --permute 1000 --out CD4_TPC_IDEAS/enrich_${state}_eQTL.txt
EOF
done < IDEAS_states.txt
```

results:

```
0 Quies
Enrichment analysis
  * #observed overlaps = 565 / 1574 (35.90%)
  * #null overlaps = 585.42 +/- 18.67 (37.19% +/- 1.19%)
  * empirical p-value = 0.296703
  * Odd ratio = 0.9467 [1.0455,0.8593]

1 ReprPCwk (?)
Enrichment analysis
  * #observed overlaps = 229 / 1574 (14.55%)
  * #null overlaps = 241.49 +/- 14.59 (15.34% +/- 0.93%)
  * empirical p-value = 0.408591
  * Odd ratio = 0.9417 [1.0938,0.8260]

2 Tx Weak
Enrichment analysis
  * #observed overlaps = 126 / 1574 (8.01%)
  * #null overlaps = 145.77 +/- 11.76 (9.26% +/- 0.75%)
  * empirical p-value = 0.106893
  * Odd ratio = 0.8511 [1.0265,0.7187]

3 Het Weak
Enrichment analysis
  * #observed overlaps = 70 / 1574 (4.45%)
  * #null overlaps = 62.63 +/- 7.45 (3.98% +/- 0.47%)
  * empirical p-value = 0.348651
  * Odd ratio = 1.1163 [1.4797,0.9049]

4 Enh
Enrichment analysis
  * #observed overlaps = 117 / 1574 (7.43%)
  * #null overlaps = 89.22 +/- 9.34 (5.67% +/- 0.59%)
  * empirical p-value = 0.004995
  * Odd ratio = 1.3399 [1.6752,1.1010]

5 Tx
Enrichment analysis
  * #observed overlaps = 182 / 1574 (11.56%)
  * #null overlaps = 166.59 +/- 11.78 (10.58% +/- 0.75%)
  * empirical p-value = 0.206793
  * Odd ratio = 1.1016 [1.2984,0.9467]

6 EnhG
Enrichment analysis
  * #observed overlaps = 56 / 1574 (3.56%)
  * #null overlaps = 64.83 +/- 7.84 (4.12% +/- 0.50%)
  * empirical p-value = 0.286713
  * Odd ratio = 0.8704 [1.1244,0.6712]

7 ZNF/Rpts
Enrichment analysis
  * #observed overlaps = 32 / 1574 (2.03%)
  * #null overlaps = 16.22 +/- 3.94 (1.03% +/- 0.25%)
  * empirical p-value = 0.000999001
  * Odd ratio = 2.0208 [3.6086,1.3403]

8 TSSA flank
Enrichment analysis
  * #observed overlaps = 30 / 1574 (1.91%)
  * #null overlaps = 32.36 +/- 5.56 (2.06% +/- 0.35%)
  * empirical p-value = 0.768232
  * Odd ratio = 0.9363 [1.3707,0.6918]

9 Het
Enrichment analysis
  * #observed overlaps = 14 / 1574 (0.89%)
  * #null overlaps = 8.13 +/- 2.77 (0.52% +/- 0.18%)
  * empirical p-value = 0.0529471
  * Odd ratio = 1.7567 [4.6996,1.0000]

10 TSSA
Enrichment analysis
  * #observed overlaps = 15 / 1574 (0.95%)
  * #null overlaps = 14.19 +/- 3.81 (0.90% +/- 0.24%)
  * empirical p-value = 0.91009
  * Odd ratio = 1.0721 [2.1539,0.6788]

11 Enh Biv
Enrichment analysis
  * #observed overlaps = 13 / 1574 (0.83%)
  * #null overlaps = 12.34 +/- 3.75 (0.78% +/- 0.24%)
  * empirical p-value = 0.912088
  * Odd ratio = 1.0840 [2.1764,0.6471]

12 Het/ReprPC
  * #observed overlaps = 9 / 1574 (0.57%)
  * #null overlaps = 14.75 +/- 3.90 (0.94% +/- 0.25%)
  * empirical p-value = 0.172827
  * Odd ratio = 0.5977 [1.1257,0.3878]

13 RepR
Enrichment analysis
  * #observed overlaps = 3 / 1574 (0.19%)
  * #null overlaps = 10.63 +/- 3.21 (0.68% +/- 0.20%)
  * empirical p-value = 0.014985
  * Odd ratio = 0.2987 [0.5992,0.1651]

14 TSS weak
Enrichment analysis
  * #observed overlaps = 1 / 1574 (0.06%)
  * #null overlaps = 2.74 +/- 1.66 (0.17% +/- 0.11%)
  * empirical p-value = 0.470529
  * Odd ratio = 0.3329 [inf,0.1423]

15 TSS Biv
Enrichment analysis
  * #observed overlaps = 6 / 1574 (0.38%)
  * #null overlaps = 7.24 +/- 2.75 (0.46% +/- 0.17%)
  * empirical p-value = 0.846154
  * Odd ratio = 0.8566 [3.0077,0.4595]

16 TxRepR
Enrichment analysis
  * #observed overlaps = 2 / 1574 (0.13%)
  * #null overlaps = 0.66 +/- 0.79 (0.04% +/- 0.05%)
  * empirical p-value = 0.276723
  * Odd ratio = inf [inf,1.0000]

17 EnhGA
Enrichment analysis
  * #observed overlaps = 54 / 1574 (3.43%)
  * #null overlaps = 39.57 +/- 6.26 (2.51% +/- 0.40%)
  * empirical p-value = 0.036963
  * Odd ratio = 1.3983 [1.9616,1.0195]

18 Enh/Het
Enrichment analysis
  * #observed overlaps = 1 / 1574 (0.06%)
  * #null overlaps = 3.04 +/- 1.73 (0.19% +/- 0.11%)
  * empirical p-value = 0.400599
  * Odd ratio = 0.3329 [inf,0.1661]

19 Enh/ReprPC
Enrichment analysis
  * #observed overlaps = 0 / 1574 (0.00%)
  * #null overlaps = 0.32 +/- 0.56 (0.02% +/- 0.04%)
  * empirical p-value = 1
  * Odd ratio = -nan [-nan,0.0000]
```

eQTLs were enriched in the following:
4 Enh
7 ZNF/Rpts
17 EnhGA


12. Overlap with Previous Studies


A study looking at autoimmune diseases.

>"The cohort specifically consists of 162 African American subjects of European and African ancestry (AA), 155 East Asian subjects of Chinese, Japanese, or Korean ancestry (EA), and 377 Caucasian subjects of European ancestry (EU). The median age of the subjects is 24. (After genotype and expression QC, total number individuals analyzed are 461, of which 401 have monocyte data and 407 have CD4+ T cell data, respectively)"


We wanted to see how many SNPs were found as eQTLs in both studies
```bash
cd /home/greally-lab/Deepa_Andrew/eQTL/Validate2/eQTL_datasets/
cp ../hc_genes/eQTL_final3.txt  # SNP is in 8th column

cd Boston
 # SNP name is in the first column
wc -l tableS6_aa_cd4T_cis_fdr05.tsv
 # 38,238 tableS6_aa_cd4T_cis_fdr05.tsv snps in the study
```

Now to look at the overlap

grabbing the total number of variants
```bash
cd /home/greally-lab/Deepa_Andrew/eQTL/Genotype2
awk 'NR>27 {OFS="\t"; print $1,$2,$3}' SplX_converted_noHH_eQTL_snpfilt_hg38_sort.vcf > SplX_converted_noHH_eQTL_snpfilt_hg38_sort_probes
wc -l SplX_converted_noHH_eQTL_snpfilt_hg38_sort_probes
 # 497,037 SplX_converted_noHH_eQTL_snpfilt_hg38_sort_probes

```
```R
options(StringsAsFactors=F)
mine <- read.table("../eQTL_final3.txt")
mine$V8 <- as.character(mine$V8)

boston_aa <- read.table("tableS6_aa_cd4T_cis_fdr05.tsv",header=T)
boston_aa$SNP <- as.character(boston_aa$SNP)

sum(mine$V8 %in% boston_aa$SNP)  # 337
 #  0.2141042

boston_eu <- read.table("tableS4_eu_cd4T_cis_fdr05.tsv",header=T)
boston_eu$SNP <- as.character(boston_eu$SNP)
sum(mine$V8 %in% boston_eu$SNP)  # 650
  #  0.4129606

boston_ea <- read.table("tableS5_ea_cd4T_cis_fdr05.tsv",header=T)
boston_ea$SNP <- as.character(boston_ea$SNP)
sum(mine$V8 %in% boston_ea$SNP)  # 304
  #  0.1931385

sum(boston_eu$SNP %in% boston_aa$SNP) # 30877 (38237)
sum(boston_eu$SNP %in% boston_ea$SNP) # 34230/242467 14% (42756)
```


### creating manhattan plots
in `/home/greally-lab/Deepa_Andrew/eQTL/Validate2/TSS_density`
make Knit_manhat_eQTL.R file
```
library(knitr)
library(rmarkdown)
rmarkdown::render("eQTL_downstream_analysis2.Rmd")
```

```bash
qsub -S /bin/bash -N knit_manhat -j y -cwd -l h_vmem=70G << EOF
module load R/3.4.0/gcc.4.7.4
module load pandoc/1.19.2.1
module load texlive/2016
cd /home/greally-lab/Deepa_Andrew/eQTL/Validate2/TSS_density
R CMD BATCH --no-restore Knit_manhat_eQTL.R
EOF
```
