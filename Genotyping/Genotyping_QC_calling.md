# Genotype Workflow: Final Version
##  Andrew Johnston
## 02/27/18



### Table of Contents

1. [Sexchecking Samples](#sexcheck)
2. [QC Genotyping array pt. 1](#plateqc1)
3. [Remove multi-mapping probes](#mmprobes)
4. [QC Genotyping array pt. 2](#plateqc2)
5. [Convert Alleles to Ref Strand](#convertstrand)
6. [Plink Filtering pt. 1](#plinkfilt1)
7. [Matching genotype with sequencing data](#genomatch)
8. [Plink Filtering pt. 1](#plinkfilt2)
9. [PCA of all samples](#eigen)
10. [eQTL Filtering](#eQTLfilt)
11. [mQTL Filtering](#mQTLfilt)


<a name="prep"/>

### Preparing Genotype data

Received raw and analyzed (Genome Studio analyzed and genotype call) files from the Mt. Sinai core. The samples were split and ran on two different plates. The first plate had 90 samples and the second plate had 22. The first plate data (full data table and final report) was placed in `/home/greally-lab/Deepa_Andrew/MEGAarray2/`. The data from the second plate was placed in `/home/greally-lab/Deepa_Andrew/MEGAChips_repeats_Jan2018/` Specifically, the full data table and final report can be found in `/home/greally-lab/Deepa_Andrew/MEGAChips_repeats_Jan2018/DeepaRastogi_22RepeatSamples_Jan2018`.

The sexcheck via LRR signal and initial QC/curation of SNPs from the array plates was performed in parallel. Described here is the QC described for the second plate but the first plate was QC'd in the same manner leading up to the merging of the two plates.


<a name="sexcheck"/>

1. examine and sexcheck the second plate

Examine the full data table
```bash
wc -l Full_Data_Table_18.txt
 # 1,779,820 lines in Full_Data_Table_18.txt; 1,779,819 probes plus the header line

awk 'NR >1 && $6<.7 {print $6}' Full_Data_Table_18.txt | wc -l
 # 80,654 probes were below the "good" gentrain score

```

Check to see if sex (as determined by the mean probe signal from the X and Y chromosomes) matches the metadata for the samples. First we determine the sex using the array data. Please skip to [plate QC](#plateqc) to generate a file necessary for this analysis.

getting the LRR (normalized signal intensity for the X and Y chromosomes):
```bash
 # get the correct chr names for X and Y (who knows if naming convention differs)
awk '{print $4}' Full_Data_Table_18.txt | sort -u | less  
 # X and Y used in this dataset
awk '$4=="Y" || $4=="X"' Full_Data_Table_18.txt > Full_Data_Table_XY.txt
```

Then I need to get the columns for LRR for each sample. I did this using R
```R
 # Get columns names to filter
Full_Data_Table_header <- read.table("Full_Data_Table_18_header.txt", sep="\t")
Full_Data_Table_header <- t(Full_Data_Table_header)
 # check file
head(Full_Data_Table_header, 20)
 # change row names to just be ordered numbers.
rownames(Full_Data_Table_header)<- 1:nrow(Full_Data_Table_header)

 # Now I want to filter the columns for LRR
fil_LRR_cols <- c("Name","Chr", "Position","Log R Ratio")
idx_LRR_name <- grep(pattern = paste(fil_LRR_cols, collapse="|"), Full_Data_Table_header)
Full_Data_Table_header_fil_LRR <- Full_Data_Table_header[idx_LRR_name,]
length(Full_Data_Table_header_fil_LRR) # 25

 # Paste a $ sign in front for easy awk filtering
Col_LRR <- paste("$", names(Full_Data_Table_header_fil_LRR), sep = "")
Col_LRR <- as.data.frame(Col_LRR)
head(Col_LRR)

 # write out table
write.table(Col_LRR, "XY_LRR_columns.txt", quote=F,
            append = F, sep = "\t", row.names = F, col.names = F)
```
Using `XY_LRR_columns.txt` I filtered the columns
```bash
 # Now to filter using the columns found in R
awk '{ for (i=1;i<=NF;i++ ) printf $i "," }' XY_LRR_columns.txt
 # now to filter the columns using awk
awk -F "\t" 'BEGIN {OFS="\t";} {print $2,$4,$5,$20,$38,$56,$74,$92,$110,$128,$146,$164,$182,$200,$218,$236,$254,$272,$290,$308,$326,$344,$362,$380,$398}' Full_Data_Table_XY.txt > Full_data_XY_LRR.txt

 # get the Snps that were used in downstream analysis
grep -v "NC" Full_data_zcall_noMM_noNA.txt > Full_data_zcall_noMM_noNA_noNC.txt
wc -l Full_data_zcall_noMM_noNA_noNC.txt
 # 1,590,579 Full_data_zcall_noMM_noNA_noNC.txt
awk '{print $1}' Full_data_zcall_noMM_noNA_noNC.txt > Full_data_zcall_noMM_noNA_noNC.snps
```

Now I will make a plot of the signal intensity LRR (Log R Ratio) of the X and Y chromomes using R:
```R
library(data.table)
library(scales)
XY_table <- fread("Full_data_XY_LRR.txt")
XY_table <- as.data.frame(XY_table)
nrow(XY_table) # 55575
 # filter out the Multimapping and NA snps
snp_pass_qc <- fread("Full_data_zcall_noMM_noNA.snps", header=T, sep="\t")
snp_pass_qc <- as.data.frame(snp_pass_qc)

idx_pass <- which(XY_table$V1 %in% snp_pass_qc$Name)
length(idx_pass) # 49036

XY_table_pass <- XY_table[idx_pass,]

X_LRR <- XY_table_pass[XY_table_pass$V2=="X",]
dim(X_LRR)
Y_LRR <- XY_table_pass[XY_table_pass$V2=="Y",]
 # Y_LRR <- replace(Y_LRR,is.na(Y_LRR), 0) # using 0 will skew
Y_LRR<- na.omit(Y_LRR)
dim(Y_LRR) # 2506
options(scipen=999)
X_LRR_mean <- apply(X_LRR[,-c(1:3)], 2, function(x) {mean(x)} )
Y_LRR_mean <- apply(Y_LRR[,-c(1:3)], 2, function(x) {mean(x)} )

plot(x=X_LRR_mean, y=Y_LRR_mean)
 # appears to be a good mix of male and female.

 # get sample information
samp_info <- read.csv("plate2_sex_table.csv", header=F)
head(samp_info)
 #combine LRR means and sample info (sex)
samp_XY_LRR <- cbind(samp_info, X_LRR_mean, Y_LRR_mean)
head(samp_XY_LRR)

 # color males blue and females red
samp_XY_LRR$colors <- "blue"
samp_XY_LRR$colors[samp_XY_LRR$V2 == 1] <- "red"

 # make plot
pdf(file = "LRR_plot_plate2.pdf", width = 7, height = 5, family="ArialMT")
plot(x=samp_XY_LRR$X_LRR_mean, y=samp_XY_LRR$Y_LRR_mean,
     col = alpha(samp_XY_LRR$colors, .5), pch=16, cex=.75, xlab="Mean X intensity",
     ylab = "Mean Y intensity", main= "Samples match designated sex: Plate #2")
dev.off()
```
You can see how the recorded gender matches the X and Y intensities seen on the plate here: [plate2](Figures/LRR_plot_plate2.pdf). A plot generated for the first plate can be found here: [plate1](Figures/LRR_plot_plate1.pdf). A7 was a mismatch because the sample was a duplication of A9 who was female.

<a name="plateqc1"/>

2. QC Genotyping array

I used R to extract the correct column numbers for transforming the data into a format ready for zcall to convert into plink format

```bash
 # I need to extract the header to filter for the correct columns needed
head -n 1 Full_Data_Table_18.txt > Full_Data_Table_18_header.txt
```

```R
 # Get columns to filter for zcall table
Full_Data_Table_header <- read.table("Full_Data_Table_18_header.txt", sep="\t")
Full_Data_Table_header <- t(Full_Data_Table_header)
 # check file
head(Full_Data_Table_header, 20)
 # change row names to just be ordered numbers.
rownames(Full_Data_Table_header)<- 1:nrow(Full_Data_Table_header)

 # grab the columns wanted.
fil_cols <- c("Name","Chr", ".GType","Position",".X",".Y")
idx_name <- grep(pattern = paste(fil_cols, collapse="|"), Full_Data_Table_header)
Full_Data_Table_header_fil_1 <- Full_Data_Table_header[idx_name,]

 # raw scores and the custom GType need to be removed.
xraw_names <- c(".X Raw",".Y Raw", "Custom GType")
idx_raw <- grep(pattern= paste(xraw_names, collapse = "|"), Full_Data_Table_header_fil_1)
Full_Data_Table_header_fil_2 <- Full_Data_Table_header_fil_1[-idx_raw]

 # Visual check
head(Full_Data_Table_header_fil_2,12)
tail(Full_Data_Table_header_fil_2)

 # Sample Check
 # what are all of the included samples? Used theta because each sample only
 # has one column matching it
header_Theta <- grep(".Theta", Full_Data_Table_header, value = T)
head(header_Theta)
length(header_Theta) # 22
str(header_Theta)
length(header_Theta)
header_Theta_substr <- do.call(rbind,strsplit(header_Theta, "[.]"))[,1]
header_Theta_substr
 # "A18" "A19" "A22" "A55" "A56" "A57" "A58" "A59" "A60" "A61" "A62" "B43" "B44" "B45" "B46"
 # "B47" "B48" "B49" "B50" "B51" "B52" "B53"
 # all of which are what we want

 # Paste a $ sign in front for easy awk filtering
Col_zcall <- paste("$", names(Full_Data_Table_header_fil_2), sep = "")
Col_zcall <- as.data.frame(Col_zcall)
head(Col_zcall)

write.table(Col_zcall, "Gtype_XY_columns.txt", quote=F,
            append = F, sep = "\t", row.names = F, col.names = F)
```

I then used `awk '{ for (i=1;i<=NF;i++ ) printf $i "," }' Gtype_XY_columns.txt` to generate the text necessary to build the following `awk` command, used to filter the large table:

```bash
awk -F "\t" 'BEGIN {OFS="\t";} {print $2,$4,$5,$11,$17,$18,$29,$35,$36,$47,$53,$54,$65,$71,$72,$83,$89,$90,$101,$107,$108,$119,$125,$126,$137,$143,$144,$155,$161,$162,$173,$179,$180,$191,$197,$198,$209,$215,$216,$227,$233,$234,$245,$251,$252,$263,$269,$270,$281,$287,$288,$299,$305,$306,$317,$323,$324,$335,$341,$342,$353,$359,$360,$371,$377,$378,$389,$395,$396}' Full_Data_Table_18.txt > Full_data_zcall.txt
```

<a name="mmprobes"/>

3. Remove Multimapping probes

Remove multimapping probes (combination of my own alignment analyses + missing and multiple probe files from [wellcome trust institute](http://www.well.ox.ac.uk/~wrayner/strand/) (both build37 and build 38) + the multiple mappers from [illumina](https://support.illumina.com/content/dam/illumina-support/documents/downloads/productfiles/multiethnic-global/multi-ethnic-global-8-b1-mapping-comment.zip) )

a. Determine multimapping probes (my method)
First I need to extract the probe sequences and create a fasta file
```bash
awk -F "," 'NR>8 && NR<1779828 {print ">"$2"\n"$6}' Multi-EthnicGlobal_A1.csv > MEGA_probe.fasta
```

a\ is the append command for sed, http://www.funtoo.org/Sed_by_Example,_Part_2 I used this because I was worried about the line breaks (especially, the last line)

```bash
sed -i -e '$a\' MEGA_probe.fasta
```
Then I aligned the probes to an hg19 build, sorted the reads and indexed them.

```bash
qsub -S /bin/bash -N Probe_align -j y -q highmem.q -cwd -pe smp-highmem 10 -l h_vmem=5.6G << EOF
module load samblaster/v.0.1.22/gcc.4.4.7
module load bwa/0.7.15/gcc.4.4.7
module load samtools/1.2/gcc.4.4.7
array="MEGA_probe"
bwa mem -t 10 -V -M -a /home/greally-lab/indexes/hg19/bwa/hg19 ${array}.fasta > ${array}.sam
samtools view -Sb ${array}.sam | samtools sort -@ 3 -T temp_ -O sam -o ${array}_sort.sam
samtools view -Sb ${array}_sort.sam | samtools index -
rm ${array}.sam
EOF
```

I defined multi-mappers as having perfect 50 bp matching at more than 1 region of the genome

```bash
awk 'NR>96 && $6=="50M" {print $1}' MEGA_probe_sort.sam | sort | uniq -d | wc -l
 #81,024 probes are multi-mappers.
 #Saving the multi-mapping probe names
awk 'NR>96 && $6=="50M" {print $1}' MEGA_probe_sort.sam | sort | uniq -d > BWA_Duplicated.txt
```

b. Extract Illumina Defined Multimappers downloaded the [Infinium Multi-Ethnic Global-8 v1.0 Mapping Comments File](https://support.illumina.com/content/dam/illumina-support/documents/downloads/productfiles/multiethnic-global/multi-ethnic-global-8-b1-mapping-comment.zip) from the illumina support files for the MEGA array. I filtered for snps with a comment.

```bash
awk -F "\t" '$2!=""' Multi-EthnicGlobal_B1_MappingComment.txt > Illumina_Multimappers.txt
wc -l Illumina_Multimappers.txt
 # 15,783 (includes header)
 # 15,772 of my multi-mappers are also in illumina's multi-mapper file, I'm only missing 10 that they annotated
```

c. Retrieving the Wellcome Trust Inst. Multimappers

The Wellcome Trust also found multimapping and probes without a full match alignement (deemed missing) for the array on 3 different references (hg18, hg19, and hg38). Found [here](http://www.well.ox.ac.uk/~wrayner/strand/). 57% of my multimappers are in their hg19 multi-mapping file. You can see a count of the probes recorded as multimapping or missing.

```bash
wc -l Multi-EthnicGlobal_A1*
    # 384 Multi-EthnicGlobal_A1-b36.miss
    # 70121 Multi-EthnicGlobal_A1-b36.multiple
    # 103 Multi-EthnicGlobal_A1-b37.miss
    # 70083 Multi-EthnicGlobal_A1-b37.multiple
    # 461 Multi-EthnicGlobal_A1-b38.miss
    # 76635 Multi-EthnicGlobal_A1-b38.multiple
```

d. combining the multimapping and missing probes (minus hg36/hg18)

```R
library(data.table)

 # Read in my multi-mapper SNP IDs
MM_snp_ID <- read.table("BWA_Duplicated.txt")
head(MM_snp_ID)
nrow(MM_snp_ID) # 81024

 # Read in illumina multi-mapper SNP IDs
il_MM_snp_ID <- read.table("Illumina_Multimappers.txt", sep="\t", header=T)
head(il_MM_snp_ID)
nrow(il_MM_snp_ID) # 15782

 # Read in Wellcome multi-mapper SNP IDs
WThg37_MM_snp_ID <- fread("Multi-EthnicGlobal_A1-b37.multiple", sep="\t", header=F, fill=T)
head(WThg37_MM_snp_ID)
nrow(WThg37_MM_snp_ID) # 70,083
WThg38_MM_snp_ID <- fread("Multi-EthnicGlobal_A1-b38.multiple", sep="\t", header=F, fill=T)
head(WThg38_MM_snp_ID)
nrow(WThg38_MM_snp_ID) # 76,635

 # Read in Wellcome missing SNP IDs
WThg37_miss_snp_ID <- fread("Multi-EthnicGlobal_A1-b37.miss", sep=" ", header=T, fill=T)
head(WThg37_miss_snp_ID$best)
nrow(WThg37_miss_snp_ID) # 102
WThg38_miss_snp_ID <- fread("Multi-EthnicGlobal_A1-b38.miss", sep=" ", header=T, fill=T)
head(WThg38_miss_snp_ID$best)
nrow(WThg38_miss_snp_ID) # 460

 # combine all of the SNP lists
unfil_SNP <- c(MM_snp_ID$V1, il_MM_snp_ID$Name, WThg37_MM_snp_ID$V1,
                    WThg38_MM_snp_ID$V1, WThg37_miss_snp_ID$best,
                    WThg38_miss_snp_ID$best)
sum(duplicated(unfil_SNP)) # 84433
length(unfil_SNP) # 244,096
fil_SNP <- unfil_SNP[!(duplicated(unfil_SNP))]
length(fil_SNP) # 108,412

 # write out the probe names to be filtered
write.table(fil_SNP,"Combined_multi_miss_snp.txt", append = F,
            quote=F, row.names = F, col.names = F, sep = "\t")
```

The multimapping probe IDs are stored in `Combined_multi_miss_snp.txt`. Using this file, the following R script filtered the full data table:

```R
 # filter the multi-mapping SNPs from the large Zcall table.
options(stringsAsFactors = F)
library(data.table)

 # Read in the full table
Zcall_data_table <- fread("Full_data_zcall.txt", header=T)
Zcall_data_table <- as.data.frame(Zcall_data_table)

 # Read in the multi-mapper SNP IDs
fil_snp_ID <- read.table("Combined_multi_miss_snp.txt")
head(fil_snp_ID)
nrow(fil_snp_ID) # 108,412

idx_MM_snp <- which(Zcall_data_table[,1] %in% fil_snp_ID$V1)

Zcall_data_table_filtered <- Zcall_data_table[-idx_MM_snp,]
write.table(Zcall_data_table_filtered,"Full_data_zcall_noMM.txt", append = F, quote=F, row.names = F,
            col.names = T, sep = "\t")
```

<a name="plateqc2"/>

4. Genotyping QC pt. 2

a. Removing SNPs with values of NA for either X or Y for any sample because zcall cannot determine the genotype call for these SNPs.

```bash
awk '$0 ~ "NA" {sum+=1} END {print sum}' Full_data_zcall_noMM.txt
 # 1012 rows have NA values for X or Y

 # removing the SNPs in which samples have NA values for X or Y
grep -v "NA" Full_data_zcall_noMM.txt > Full_data_zcall_noMM_noNA.txt

wc -l Full_data_zcall_noMM_noNA.txt
 # 1,670,396 Full_data_zcall_noMM_noNA.txt; Snps (header included) left in the file
```

b. Removing SNPs with no calls (NC) in greater than 5% of the samples. Therefore, if more than 1 NC per snp for the second plate, but more than 4 NCs for the second plate.


How many of these snps are mostly NCs?
Count the NCs in each of these rows and look at the distribution.
```bash
 # grab the NC rows.
grep "NC" Full_data_zcall_noMM_noNA.txt > Full_data_zcall_noMM_noNA_NC.txt
```
Then use R to look at the distribution and then remove "poor quality" SNPs

```R
 # Look at the number of non-calls (NCs) in the genotype data

 # set options
options(scipen = 999, stringsAsFactors = F)
library(data.table)

 # read in the NC rows
Nc_rows <- fread("Full_data_zcall_noMM_noNA_NC.txt", header=FALSE)
Nc_rows <- as.data.frame(Nc_rows)
nrow(Nc_rows) # 79817

 # calculate the number of NCs per
NC_per_row <- apply(X = Nc_rows, MARGIN = 1, FUN = function(x){sum(x=="NC")})
length(NC_per_row)

pdf(file = "NC_per_SNP_dist.pdf", width = 7, height = 5, family="arialMT")
hist(NC_per_row, main="Distribution of NCs per SNPs with an NC", xlab = "# of Samples")
 # there are a good 30k SNPs for which all samples have NC.
dev.off()

 # get the SNPs to remove (all SNPs with >1 NC)
Nc_rows$NC <- NC_per_row
rm_NC_Snps <- Nc_rows[Nc_rows$NC>1, 1]
head(rm_NC_Snps)
nrow(rm_NC_Snps) # 47,727
 # length(rm_NC_Snps)

 # Read in the full table
Zcall_data_table <- fread("Full_data_zcall_noMM_noNA.txt", header=T)
Zcall_data_table <- as.data.frame(Zcall_data_table)

 # filter out the >1 NC SNPs
idx_NC_snp <- which(Zcall_data_table[,1] %in% rm_NC_Snps$V1)
 # idx_NC_snp <- which(Zcall_data_table[,1] %in% rm_NC_Snps)


Zcall_data_table_filtered <- Zcall_data_table[-idx_NC_snp,]
write.table(Zcall_data_table_filtered,"Full_data_zcall_noMM_noNA_NCfil.txt", append = F,
            quote=F, row.names = F,
            col.names = T, sep = "\t")
```

How many SNPs are left?
```bash
wc -l Full_data_zcall_noMM_noNA_NCfil.txt
 # 1,622,669 Full_data_zcall_noMM_noNA_NCfil.txt
```

c. Run zcall's QCReport
QCReport filters out samples with high NC rate.

```bash
module load zCall/050817
PY_PATH="/public/apps/zCall/zCall-master/Version3_GenomeStudio/GenomeStudio"
qsub -S /bin/bash -N CallRate_qc -cwd -l h_vmem=20G -j y << EOF
module load zCall/050817
 # hard code the path
PY_PATH="/public/apps/zCall/zCall-master/Version3_GenomeStudio/GenomeStudio"
echo $PY_PATH
python ${PY_PATH}/qcReport.py -R Full_data_zcall_noMM_noNA_NCfil.txt -C 0.99 > Full_data_zcall_noMM_noNA_NCfil_qcreport.txt
echo "Full_data_zcall_noMM_noNA.txt"
awk '{print NF}' Full_data_zcall_noMM_noNA.txt | sort -nu | tail -n 1
 # 69
echo "Full_data_zcall_noMM_noNA_NCfil_qcreport.txt"
awk '{print NF}' Full_data_zcall_noMM_noNA_NCfil_qcreport.txt | sort -nu | tail -n 1
 #  69 # Did I lose any samples? # No I didn't
EOF
```

d. Filter out SNPs not following HWE (each plate separately)
It was suggested to filter out and p_hwe < 1e-5)

```bash
 # using the filtered and call rate >95% files to generate plink compatible files to calculate HWE pvalues.
PY_PATH="/public/apps/zCall/zCall-master/Version3_GenomeStudio/GenomeStudio"

 # convert plate2 table to plink compatible file
qsub -S /bin/bash -N convertTPED -cwd -l h_vmem=15G -j y << EOF
module load zCall/050817
PY_PATH="/public/apps/zCall/zCall-master/Version3_GenomeStudio/GenomeStudio"
python ${PY_PATH}/convertReportToTPED.py -O Full_data_zcall_noMM_noNA_NCfil -R Full_data_zcall_noMM_noNA_NCfil.txt
EOF

 # use plink to filter out the hwe <0.000001 variants
qsub -S /bin/bash -N hwe_check -cwd -l h_vmem=5G -j y << EOF
module load plink/1.90b
plink --tfile Full_data_zcall_noMM_noNA_NCfil --split-x 'hg19' --hwe .000001 'midp' --make-bed --out Full_data_hwe_check_2
EOF

 #perform the same analysis on the first plate
cd /home/greally-lab/Deepa_Andrew/MEGAarray2

 # convert plate1 table to plink compatible file
qsub -S /bin/bash -N convertTPED -cwd -l h_vmem=15G -j y << EOF
module load zCall/050817
PY_PATH="/public/apps/zCall/zCall-master/Version3_GenomeStudio/GenomeStudio"
python ${PY_PATH}/convertReportToTPED.py -O Full_data_zcall_noMM_noNA_NCfil -R Full_data_zcall_noMM_noNA_NCfil.txt
EOF

 # use plink to filter out the hwe <0.000001 variants
qsub -S /bin/bash -N hwe_check -cwd -l h_vmem=5G -j y << EOF
module load plink/1.90b
plink --tfile Full_data_zcall_noMM_noNA_NCfil --split-x 'hg19' --hwe .000001 'midp' --make-bed --out Full_data_hwe_check_2
EOF
```

After HWE filtering, the plates had 1,647,596 variants (plate 1) and 1,622,657 variants (plate 2)

```R
 # code to filter out those SNPs breaking HWE p < 0.000001

 #libraries
library(data.table)

 #read-in the original SNP file
Data_table <- fread("Full_data_zcall_noMM_noNA_NCfil.txt")
Data_table <- as.data.frame(Data_table)

 # read-in filtered SNPs
hwe_passing_tab <- fread("Full_data_hwe_check_2.bim")
hwe_passing_tab <- as.data.frame(hwe_passing_tab)

 #extract only those variants which pass filtering
idx_pass <- which(Data_table$Name %in% hwe_passing_tab$V2)
Data_table_hwe <- Data_table[idx_pass,]

 #now do this for the first array plate
 #read-in the original SNP file
Data_table_1 <- fread("../../MEGAarray2/Full_data_zcall_noMM_noNA_NCfil.txt")
Data_table_1 <- as.data.frame(Data_table_1)

 # read-in filtered SNPs
hwe_passing_tab_1 <- fread("../../MEGAarray2/Full_data_hwe_check_2.bim")
hwe_passing_tab_1 <- as.data.frame(hwe_passing_tab_1)

 #extract only those variants which pass filtering
idx_pass_1 <- which(Data_table_1$Name %in% hwe_passing_tab_1$V2)
Data_table_1_hwe <- Data_table_1[idx_pass_1,]

 # now merge the two tables
Data_table_com <- merge(Data_table_1_hwe, Data_table_hwe, by="Name", sort=FALSE, all=FALSE)
dim(Data_table_com) # 1602747     341

 # clean up columns
idx_rmCol <- which(colnames(Data_table_com) %in% c("Chr.y","Position.y"))
Data_table_com <- Data_table_com[,-idx_rmCol]

 # write out the file
write.table(Data_table_com,"Full_data_zcall_noMM_noNA_NCfil_noHWE_mergePlates.txt", append = F, quote=F, row.names = F, col.names = T, sep = "\t")
```

Now I'll convert the merged file to plink and begin plink filtering

```bash
module load zCall/050817
PY_PATH="/public/apps/zCall/zCall-master/Version3_GenomeStudio/GenomeStudio"
qsub -S /bin/bash -N convertTPED -cwd -l h_vmem=15G -j y << EOF
module load zCall/050817
PY_PATH="/public/apps/zCall/zCall-master/Version3_GenomeStudio/GenomeStudio"
python ${PY_PATH}/convertReportToTPED.py -O Merged_plates_final -R Full_data_zcall_noMM_noNA_NCfil_noHWE_mergePlates.txt
EOF
```
This generates a generic .tfam file which needs to be edited with actual metadata.

```R
 #fixing the tfam with appropriate sample data

 # read in the zcall generated tfam file
tfam <- read.table("Merged_plates_final.tfam")
head(tfam)
nrow(tfam) #112
 # read in sample information from deepa
fam_deepa <- read.csv("Sample_table_forPLINK_all.csv")
head(fam_deepa)

 # get correct order of samples for fam file
idx_new <- match(tfam$V1, fam_deepa$Sample.ID)
fam_deepa_2 <- fam_deepa[idx_new,]
head(fam_deepa_2)
dim(fam_deepa_2)
fam_deepa_2
 # completely separated by sex except for A7 (which is a mismatched sample)

 # switch the sample.ID column with individual ID and write out
write.table(fam_deepa_2[,c(2,1,4:7)], "Merged_plates_final.tfam", quote=F,
            append = F, sep = "\t", row.names = F, col.names = F)
```

<a name="convertstrand"/>

5. Convert alleles to reference (+) Strand

Converting Alleles to the + strand, this way the variants will be in terms of the reference strand and will match sequencing datasets.

First thing is to convert all of the alleles to the + strand (reference strand), so that the allele calls match with
```bash
 # The plink filtering generates a ton of files so I'm moving the files to plink_filtering.
mkdir Plink_filtering
cp Merged_plates_final.* Plink_filtering/
cd Plink_filtering/

 #make the plink file

qsub -S /bin/bash -N Plink_splitX -cwd -l h_vmem=5G -j y << EOF
module load plink/1.90b
plink --tfile Merged_plates_final --split-x 'hg19' --make-bed --out SplX
EOF
 # 1602747 variants loaded from .bim file.
 # 112 people (64 males, 48 females) loaded from .fam.
 # Warning: 4265 het. haploid genotypes present (see SplX.hh ); many commands
 # Among remaining phenotypes, 53 are cases and 59 are controls.
```

Once in plink format, I will now convert the alleles. First we need to generate the filtered snp lists.
```
 # getting the list of snps that passed initial filters and merge
 # should be 1,602,747 snps.
awk '{print $2}' Merged_plates_final.tped > Merged_plates_final_snps.txt
 # getting the snplist (Name	SNP	ILMN Strand	Customer Strand)from the final report
cd ../
awk -v OFS='\t' 'NR>10 && $2=="A18" {print $1,$21,$22,$23}'  DeepaRastogi_22RepeatSamples_01_30_2018_FinalReport.txt > FinalReport.snplist
 #sanity check
wc -l FinalReport.snplist
 # 1779819 FinalReport.snplist
```

Now I need to generate a file which GACT will use to convert the alleles.

```R
 # Generating a filtered Annotation file for gact
library(data.table)
 # read in snp names that passed initial filters
noMM_noNA <- read.table("Plink_filtering/Merged_plates_final_snps.txt")
head(noMM_noNA)
dim(noMM_noNA) # 1602747       1
 # read in the MEGA metadata csv file
MEGA_csv <- read.csv("Multi-EthnicGlobal_A1.csv", skip=7)
head(MEGA_csv)
dim(MEGA_csv) # 1779843      21
 # clean the columns
annotation <- MEGA_csv[, c(10,1,2,3,4,21)]
head(annotation)
 # filter the annotation file
annotation_fil <- annotation[(annotation$Name %in% noMM_noNA$V1),]
head(annotation_fil)
dim(annotation_fil) # 1602747       6
 # write out the new filtered table
write.table(annotation_fil, "../../../../anjohnst/Programs/gact/references/Multi-EthnicGlobal_A1_filt/Multi-EthnicGlobal_A1_filt.annotation.txt", quote=F,
            append = F, sep = "\t", row.names = F, col.names = T)

 # Now filter the strand file too
 # read in the strand file
b37 <- fread("../../../../anjohnst/Programs/gact/references/Multi-EthnicGlobal_A1/Multi-EthnicGlobal_A1-b37.strand")
b37 <- as.data.frame(b37)
head(b37)
 # filter the strand data
b37_filt <- b37[(b37$V1 %in% noMM_noNA$V1),]
write.table(b37_filt, "../../../../anjohnst/Programs/gact/references/Multi-EthnicGlobal_A1_filt/Multi-EthnicGlobal_A1_filt-b37.strand", quote=F,
            append = F, sep = "\t", row.names = F, col.names = F)

 # Now to create the GenGen Conversion file
 # read in final report snplist
Final_report_snps <- fread("FinalReport.snplist", header = F)
Final_report_snps <- as.data.frame(Final_report_snps)
head(Final_report_snps)
dim(Final_report_snps) # 1779819       4
 # filter the snplist by snps passing filter
Final_report_snps_fil <-  Final_report_snps[(Final_report_snps$V1 %in% noMM_noNA$V1),]
dim(Final_report_snps_fil) #  1602747       4
 # recode the plus to P and minus to M
Final_report_snps_fil$ILMN <- gsub("PLUS","P",Final_report_snps_fil$V3, fixed = T)
Final_report_snps_fil$ILMN2 <- gsub("MINUS","M",Final_report_snps_fil$ILMN, fixed = T)
Final_report_snps_fil$CUST <- gsub("PLUS","P",Final_report_snps_fil$V4, fixed = T)
Final_report_snps_fil$CUST2 <- gsub("MINUS","M",Final_report_snps_fil$CUST, fixed = T)
 # check
head(Final_report_snps_fil, 40)
 # clean columns and rename them for output
Final_report_snps_fil_out <- Final_report_snps_fil[,c(1,2,6,8)]
colnames(Final_report_snps_fil_out) <- c("Name", "SNP", "ILMN Strand", "Customer Strand")
head(Final_report_snps_fil_out)
 # write the file
write.table(Final_report_snps_fil_out, "../../../../anjohnst/Programs/gact/references/Multi-EthnicGlobal_A1_filt/Multi-EthnicGlobal_A1_filt.snptable", quote=F,
            append = F, sep = "\t", row.names = F, col.names = T)
```

Since all of the necessary files are generated, we finally run GACT in it's folder
```bash
cd ~/Programs/gact
cp /home/greally-lab/Deepa_Andrew/MEGAChips_repeats_Jan2018/DeepaRastogi_22RepeatSamples_Jan2018/Plink_filtering/SplX.bim .
module load GACT/070617
./gact.sh SplX.bim  SplX_convert.bim  AB PLUS b37 b37 Multi-EthnicGlobal_A1_filt
cp SplX_convert.bim /home/greally-lab/Deepa_Andrew/MEGAChips_repeats_Jan2018/DeepaRastogi_22RepeatSamples_Jan2018/Plink_filtering/.
```

Use GACT generated bim to create update allele file for plink Need to now make a file to update the alleles with plink --update-alleles updates variant allele codes. Its input should have the following five fields:

1. Variant ID
2. first old allele code
3. second old allele code
4. New code for the first named allele
5. New code for the second named allele

```R
orig_bim <- fread("Plink_filtering/SplX.bim")
head(orig_bim)
orig_bim <- as.data.frame(orig_bim)
nrow(orig_bim) # 1602747
 # read in gact generated bim file
convert_bim <- fread("Plink_filtering/SplX_convert.bim")
head(convert_bim)
convert_bim <- as.data.frame(convert_bim)
nrow(convert_bim) # 1602747

 # The convert bim and original bim don't exactly match up. I have no idea why not.
sum(orig_bim$V2 == convert_bim$V2) # 1484150

 # make sure the original and converted bims match up in regards to snp
idx_match <- match(orig_bim$V2, convert_bim$V2)
length(idx_match)
sum(orig_bim$V2 == convert_bim$V2[idx_match]) # 1664244

 # merge the files
merge_bim <- cbind(orig_bim[,c(2,5,6)], convert_bim[idx_match,c(5,6)])
head(merge_bim)
dim(merge_bim) # 1602747  5
write.table(merge_bim, "Plink_filtering/convert_update_allele.txt", quote=F,
            append = F, sep = "\t", row.names = F, col.names = F)
```

There's an issue where monomorphic SNPs are switched using gact so that the ped file will show 0's instead of the snp. These SNPs are uninformative since the MAF is <.05; but we'll get around this error nonetheless.

```bash
 # switches the monomorphic snps back to 0 then the snp
awk -v OFS='\t' '{if ($2=="0") print $1,$2,$3,$5,$4; else print $0}' convert_update_allele.txt > convert_update_allele2.txt
 # sanity check that non-monomorphic snps are ok
awk  '$5=="0"' convert_update_allele.txt | awk '$3=="A" || $3 =="B"' | wc -l # 488603
awk  '$5=="0"' convert_update_allele2.txt | awk '$3=="A" || $3 =="B"' | wc -l # 0
```

Finally, we recode the snps
```bash
qsub -S /bin/bash -N convert_recode -cwd -l h_vmem=10G -j y << EOF
module load plink/1.90b
plink --bfile SplX --update-alleles convert_update_allele2.txt --make-bed --recode tabx --out SplX_converted
EOF
```

<a name="plinkfilt1"/>

6. Plink filtering pt. 1
First let's perform a sex check and then remove the heterozygous Y chromo calls that persist after accounting for the PARs on the sex chromosomes. There are 'awk '{print $3}' SplX_converted.hh | uniq | wc -l # 3010' het snps that persist.

```bash
qsub -S /bin/bash -N checkSex -cwd -l h_vmem=5G -j y << EOF
module load plink/1.90b
module load R/3.4.0/gcc.4.7.4
 # use --maf 0.1 for sex checks as indicated in "Illumina human exome genotyping array clustering and quality control" - Nature Protocol - Guo et al 2014
plink --bfile SplX_converted --maf 0.1 --check-sex --out SplX_converted.sex_check
Rscript Gender.R SplX_converted.sex_check.sexcheck SplX_converted.sex_check.jpeg
awk '{print $3}' SplX_converted.hh | uniq | wc -l # 6,797 SNPs to remove
plink --bfile SplX_converted --set-hh-missing --make-bed --out SplX_converted_noHH
plink --bfile SplX_converted_noHH --maf 0.1 --check-sex --out SplX_converted_noHH.sex_check
 # Gender.R is from https://github.com/slzhao/ExonChipProcessing/
Rscript Gender.R SplX_converted_noHH.sex_check.sexcheck SplX_converted_noHH.sex_check.jpeg
EOF

 #how many SNPs remain on X and Y?
awk '$1==25' SplX_converted_noHH.bim | wc -l # 8 on chrXY
awk '$1==24' SplX_converted_noHH.bim | wc -l # 80 on chrY
awk '$1==23' SplX_converted_noHH.bim | wc -l # 33990 on chrX
```

A7 is a sample mismatch/duplication that will be removed. A30 and A1 are borderline female, so we'll proceed. Next we simply want to match sequencing libraries with individuals. Therefore, before filtering out anyone, I perform a quick and dirty snp filter and then create a vcf file for matching.


<a name="genomatch"/>

7. Matching genotype with sequencing data

Let's match genotypes with sequenced data.

```bash
 # quick snp filter
qsub -S /bin/bash -N Filter_snp_match -cwd -l h_vmem=10G -j y << EOF
module load plink/1.90b
plink --bfile SplX_converted_noHH --maf 0.1 --mind 0.1 --geno 0.1 --hwe 0.000001 'midp' --make-bed --out SplX_converted_noHH_snpfilt_match
EOF

 # make the vcf file
qsub -S /bin/bash -N plink2vcf_match -cwd -q highmem.q -l h_vmem=10G -j y << EOF
module load plink/1.90b
plink --bfile SplX_converted_noHH_snpfilt_match --autosome --recode vcf --out SplX_converted_noHH_snpfilt_match_hg19
EOF

 # Warning: At least one VCF allele code violates the official specification;
 # other tools may not accept the file.  (Valid codes must either start with a
 # '<', only contain characters in {A,C,G,T,N,a,c,g,t,n}, be an isolated '*', or
 # represent a breakend.)
```

The genotype data is in hg19 but the gene expression and DNA methylation data is in hg38. Therefore, the vcf file will have to be lifted over. I got the vcf lift over script from https://github.com/knmkr/lift-over-vcf. It takes forever and is single threaded. You can split up the vcf, run and merge again.

```bash
 # bgzip and tabix the vcf
mkdir match_vcf # for file clutter
module load htslib/1.2.1/gcc.4.4.7
bgzip -c SplX_converted_noHH_snpfilt_match_hg19.vcf > match_vcf/SplX_converted_noHH_snpfilt_match_hg19.vcf.gz
cd match_vcf
tabix -p vcf SplX_converted_noHH_snpfilt_match_hg19.vcf.gz

 # split vcf and liftover each chromosome seperately
for chr in {1..22};
do
qsub -S /bin/bash -N vcf_lift_${chr} -cwd -l h_vmem=2G -j y << EOF
module load htslib/1.2.1/gcc.4.4.7
tabix SplX_converted_noHH_snpfilt_match_hg19.vcf.gz ${chr} > chr${chr}.vcf
module load python/2.7.8/gcc.4.4.7
module load vcftools/11172016/gcc.4.4.7
cat chr${chr}.vcf | python ~/Programs/vcf-liftOver/lift_over.py --chain hg19ToHg38 > chr${chr}_hg38.vcf
EOF
done

 # merge the lifted vcf files
for chr in {1..22};
do
echo chr${chr}_hg38.vcf >> lifted_vcfs.txt
done

 # vcf files need headers
awk 'NR<29' SplX_converted_noHH_snpfilt_match_hg19.vcf > SplX_converted_noHH_snpfilt_match_hg19.header.txt

for chr in {1..22};
do
cat SplX_converted_noHH_snpfilt_match_hg19.header.txt chr${chr}_hg38.vcf > chr${chr}_hg38_head.vcf
mv chr${chr}_hg38_head.vcf chr${chr}_hg38.vcf
done

 # concat the files
module load vcftools/11172016/gcc.4.4.7
vcf-concat -f lifted_vcfs.txt > SplX_converted_noHH_snpfilt_match_hg38.vcf

 # Note: the HELP-tagging data is aligned to a non "chr" reference,
 # therefore the chromosome vcf will only be used for eQTL and
 # a non manipulated vcf will be used for mQTL anaylsis and matching

 # add a "chr" in front of each chromosome number.
awk '{if($0 !~ /^#/) print "chr"$0; else print $0}' SplX_converted_noHH_snpfilt_match_hg38.vcf > SplX_converted_noHH_snpfilt_match_hg38_chr.vcf

  # sort the vcf
module load vcftools/11172016/gcc.4.4.7
vcf-sort SplX_converted_noHH_snpfilt_match_hg38_chr.vcf > SplX_converted_noHH_snpfilt_match_hg38_chr_sort.vcf

 # sort the non-chr vcf for mQTLs
vcf-sort SplX_converted_noHH_snpfilt_match_hg38.vcf > SplX_converted_noHH_snpfilt_match_hg38_sorted.vcf
 #also need to remove alt chromosomes
 #19_KI270938v1_alt
 #7_KI270803v1_alt
 #8_KI270821v1_alt
grep -v "_alt" SplX_converted_noHH_snpfilt_match_hg38_sorted.vcf > SplX_converted_noHH_snpfilt_match_hg38_sorted_noAlt.vcf
```

Additionally, I manually changed all one digit sample names from 1 to 01, etc. I also added chr to the contigs. This will help the scripting implemented in the qtltools match anaylsis scripts. 'vi SplX_converted_noHH_snpfilt_match_hg38_chr_sort.vcf'

Then I finished compressing and indexing the vcf for use in QTLtools.

```bash
module load htslib/1.2.1/gcc.4.4.7
bgzip SplX_converted_noHH_snpfilt_match_hg38_chr_sort.vcf
tabix -p vcf SplX_converted_noHH_snpfilt_match_hg38_chr_sort.vcf.gz
 # check to see if the tbi is readable by QTLtools
module load bcftools/1.3.1/gcc.4.4.7
bcftools view SplX_converted_noHH_snpfilt_match_hg38_chr_sort.vcf.gz | less -S

bgzip SplX_converted_noHH_snpfilt_match_hg38_sorted_noAlt.vcf
tabix -p vcf SplX_converted_noHH_snpfilt_match_hg38_sorted_noAlt.vcf.gz
 # check to see if the tbi is readable by QTLtools
module load bcftools/1.3.1/gcc.4.4.7
bcftools view SplX_converted_noHH_snpfilt_match_hg38_sorted_noAlt.vcf.gz | less -S
```

With this VCF, I moved onto matching the vcf with the RNAseq libraries. See the [eQTL_pipeline](../eQTL/QTL_pipeline.md) markdown file to see how the bam files were collected.

Here is the call to the QTLtools mbv command on all bams. Additionally,  `--filter-mapping-quality 150` is suggested but not used because I have already filtered for uniquely mapped reads during the bam file preparation

```bash
cd /home/greally-lab/Deepa_Andrew/eQTL/Bams/filt_bams
for f1 in *.bam;
do
SAMPLE="$(echo ${f1} | cut -c1-3)"
qsub -S /bin/bash -N match_${SAMPLE} -cwd -l h_vmem=10G -j y << EOF
module load QTLtools/1.1/gcc.4.7.4
QTLtools mbv --bam $f1 --vcf /home/greally-lab/Deepa_Andrew/MEGAChips_repeats_Jan2018/DeepaRastogi_22RepeatSamples_Jan2018/Plink_filtering/match_vcf/SplX_converted_noHH_snpfilt_match_hg38_chr_sort.vcf.gz --filter-mapping-quality 0 --out ../match_files2/${SAMPLE}.match.txt
EOF
done

 # the results were analyzed using the plot_sample_match.R R script
cd ../match_files2/
for f1 in *.match.txt;
do
SAMPLE="$(echo ${f1} | cut -c1-3)"
qsub -S /bin/bash -N plotMatch_${SAMPLE} -cwd -l h_vmem=1G -j y << EOF
module load R/3.4.0/gcc.4.7.4
Rscript plot_sample_match.R $f1
EOF
done
```

After looking at the results, some samples were swapped and others had poor quality RNAseq libraries, which made the matching of the correct genotype impossible.


Unfortunately, A47, B04, B23, and B17 (and B03 which appeared to match everything) all have a really low amount of reads that overlie gene annotation (poor quality RNAseq libraries). This was apparent from the QTLtools quant where all 4 samples had a low (~500k or less reads within exonic regions). Also, looking back at the original STAR ensembl mapped files, the percent of total reads that were protein coding was only 2-4% for each of the samples. These samples will be discarded from the analysis.

Additionally, B52, B48, A52, and A11 all seemed of poor quality and did not match well to their respective genotyped DNA.

Two other samples A19 and B25 were removed due to the quality of the RNAseq libraries. See the [eQTL pipeline](../eQTL/eQTL_pipeline2.md) for the full eQTL pipeline. To see the finishing of variant filtering for eQTL analysis, go to the [eQTL_prep](#eQTL).


Now with the non-chr VCF, I moved onto matching the vcf with the HELP-tagging libraries. See the [mQTL_pipeline](../mQTL/QTL_pipeline.md) markdown file to see how the bam files were collected.

```bash
 # quick snp filter
qsub -S /bin/bash -N Filter_snp_match2 -cwd -l h_vmem=10G -j y << EOF
module load plink/1.90b
plink --bfile SplX_converted_noHH --mind 0.1 --geno 0.1 --hwe 0.000001 'midp' --make-bed --out SplX_converted_noHH_snpfilt_match2
EOF

 # make the vcf file
qsub -S /bin/bash -N plink2vcf_match2 -cwd -l h_vmem=10G -j y << EOF
module load plink/1.90b
plink --bfile SplX_converted_noHH_snpfilt_match2 --autosome --recode vcf --out SplX_converted_noHH_snpfilt_match2_hg19
EOF

 # Warning: At least one VCF allele code violates the official specification;
 # other tools may not accept the file.  (Valid codes must either start with a
 # '<', only contain characters in {A,C,G,T,N,a,c,g,t,n}, be an isolated '*', or
 # represent a breakend.)
```

The methylation data (bam files are aligned to hg19; and without chr prefixs) so I went ahead and changed the single digit sample names to two digits and implemented the following
```bash
mkdir ../../../../mQTL/Genotype2/
cp SplX_converted_noHH_snpfilt_match2_hg19.vcf ../../../../mQTL/Genotype2/

 # gzip and tabix
module load htslib/1.2.1/gcc.4.4.7
bgzip SplX_converted_noHH_snpfilt_match2_hg19.vcf
tabix -p vcf SplX_converted_noHH_snpfilt_match2_hg19.vcf.gz
 # check to see if the tbi is readable by QTLtools
module load bcftools/1.3.1/gcc.4.4.7
bcftools view SplX_converted_noHH_snpfilt_match2_hg19.vcf.gz | less -S

```

Now with this vcf let's match the files.

```bash
cd /home/greally-lab/Deepa_Andrew/mQTL/Bams/sorted
for f1 in *.bam;
do
SAMPLE="$(echo ${f1} | cut -c1-3)"
qsub -S /bin/bash -N match_${SAMPLE} -cwd -l h_vmem=10G -j y << EOF
module load QTLtools/1.1/gcc.4.7.4
QTLtools mbv --bam $f1 --vcf  /home/greally-lab/Deepa_Andrew/mQTL/Genotype2/SplX_converted_noHH_snpfilt_match2_hg19.vcf.gz --filter-mapping-quality 0 --filter-keep-duplicates --filter-base-quality 5 --filter-minimal-coverage 10 --out ../match_files2/${SAMPLE}.match.txt
EOF
done

 # cd and grab plot script
cd ../match_files2
cp ../../../eQTL/Bams/match_files2/plot_sample_match.R .

 #Now to analyze the results in R
for f1 in *.match.txt;
do
SAMPLE="$(echo ${f1} | cut -c1-3)"
qsub -S /bin/bash -N plotMatch_${SAMPLE} -cwd -l h_vmem=1G -j y << EOF
module load R/3.4.0/gcc.4.7.4
Rscript plot_sample_match.R $f1
EOF
done

 # redid match with lower coverage threshold for poor-er quality BAMS
for f1 in A07.bam A22.bam A23.bam A24.bam A52.bam B18.bam B21.bam B23.bam B24.bam;
do
SAMPLE="$(echo ${f1} | cut -c1-3)"
qsub -S /bin/bash -N match_${SAMPLE} -cwd -l h_vmem=10G -j y << EOF
module load QTLtools/1.1/gcc.4.7.4
QTLtools mbv --bam $f1 --vcf  /home/greally-lab/Deepa_Andrew/mQTL/Genotype2/SplX_converted_noHH_snpfilt_match2_hg19.vcf.gz --filter-mapping-quality 0 --filter-keep-duplicates --filter-base-quality 5 --filter-minimal-coverage 5 --out ../match_files2/${SAMPLE}.match.txt
EOF
done

for f1 in A07.match.txt A22.match.txt A23.match.txt A24.match.txt A52.match.txt B18.match.txt B21.match.txt B23.match.txt B24..match.txt;
do
SAMPLE="$(echo ${f1} | cut -c1-3)"
qsub -S /bin/bash -N plotMatch_${SAMPLE} -cwd -l h_vmem=1G -j y << EOF
module load R/3.4.0/gcc.4.7.4
Rscript plot_sample_match.R $f1
EOF
done
```

Finally, the samples were matched. Here's the final check

```bash
bgzip SplX_converted_noHH_snpfilt_match_hg38_sorted_noAlt.vcf
tabix -p vcf SplX_converted_noHH_snpfilt_match_hg38_sorted_noAlt.vcf.gz
 # check to see if the tbi is readable by QTLtools
module load bcftools/1.3.1/gcc.4.4.7
bcftools view SplX_converted_noHH_snpfilt_match_hg38_sorted_noAlt.vcf.gz | less -S
```


<a name="plinkfilt2"/>

8. Plink filtering pt. 2

Next we want to look at the relationships among the samples in the dataset in order to exclude relatives.

```bash
qsub -S /bin/bash -N checkRelate -cwd -l h_vmem=5G -j y << EOF
module load plink/1.90b
module load R/3.4.0/gcc.4.7.4
plink --bfile SplX_converted_noHH --maf 0.1 --not-chr 23,24,25,26 --indep-pairwise 50 5 0.1 --out SplX_converted_noHH_indep
plink --bfile SplX_converted_noHH --extract SplX_converted_noHH_indep.prune.in --genome -het --out SplX_converted_noHH_indep_genome
Rscript PlotHeterozygosity.R SplX_converted_noHH_indep_genome.het
EOF

 # how many samples show too much homozygosity? 7 samples
awk '$6 >.1 || $6 < -.1' SplX_converted_noHH_indep_genome.het | more
 # FID  IID       O(HOM)       E(HOM)        N(NM)            F
 # 19  A23        40279    3.543e+04        61724       0.1843
 #  2   A2        39044    3.543e+04        61726       0.1373
 # 61  B17        38777    3.543e+04        61719       0.1273
 # 60  B16        38479    3.537e+04        61609       0.1185
 # 26  A30        40462    3.543e+04        61727       0.1912
 # 41  A48        38679    3.544e+04        61728       0.1234
 #  1   A1        39467    3.543e+04        61711       0.1537

Rscript plot_relatedness.R SplX_converted_noHH_indep_genome.genome
```
Check out the relatedness plot [here](Figures/Relatedness_both_plates.pdf) It's a bit hard to discern the three identical twin pairs (B11/B10, A58/A59, and sample duplication A7/A9) in the bottom left corner.


Relatedness can also be observed through MDS plots
```bash
qsub -S /bin/bash -N MDS_plots -cwd -q highmem.q -l h_vmem=5G -j y << EOF
module load plink/1.90b
module load R/3.4.0/gcc.4.7.4
plink --bfile SplX_converted_noHH --read-genome SplX_converted_noHH_indep_genome.genome --cluster --mds-plot 2 --out SplX_converted_noHH_indep_mds_2
EOF

 # create a plot with this script
Rscript MDS_plink_plot.R
```
I colored the MDS plots by [related individuals](Figures/MDS_plot_col_by_related.pdf) and by [high homozygosity](Figures/MDS_plot_col_by_homozygosity.pdf).

Additionally, one can look at identical by state (IBS) test results to observe relatedness and possible stratification within the cases/controls.
```bash
qsub -S /bin/bash -N ibs_test -cwd -q highmem.q -l h_vmem=5G -j y << EOF
module load plink/1.90b
plink --bfile SplX_converted_noHH --maf 0.1 --not-chr 23,24,25,26 --indep-pairwise 50 5 0.1 --ibs-test
EOF
```

This is the output. No significant differences between case and control.
```
  Between-group IBS (mean, SD)   = 0.695855, 0.0162041
  In-group (case) IBS (mean, SD) = 0.696752, 0.01796
  In-group (ctrl) IBS (mean, SD) = 0.69544, 0.0196716
  Approximate proportion of variance between group = 2.34964e-05
  IBS group-difference empirical p-values:
     T1: Case/control less similar                p = 0.246668
     T2: Case/control more similar                p = 0.753332

   T3: Case/case less similar than ctrl/ctrl    p = 0.854841
   T4: Case/case more similar than ctrl/ctrl    p = 0.145159

   T5: Case/case less similar                   p = 0.881521
   T6: Case/case more similar                   p = 0.118479

   T7: Control/control less similar             p = 0.208618
   T8: Control/control more similar             p = 0.791382

   T9: Case/case less similar than case/ctrl    p = 0.874031
   T10: Case/case more similar than case/ctrl    p = 0.125969

   T11: Ctrl/ctrl less similar than case/ctrl    p = 0.372256
   T12: Ctrl/ctrl more similar than case/ctrl    p = 0.627744
```

Here's a table of the related individuals as found via genotyping.


Will choose one of the following pairs to remove:

| Pair	|		Expected 	|	Actual			|
|:-----:|:-----------------:|:-----------------:|
|B7/A5	|	full siblings	|	not related 	|
|A2/A3	|	unrelated		|	half siblings   |
|B11/B10|	unrelated		|	identical twins |
|A40/A41| 	full siblings	|	half siblings	|
|B32/B33|	full siblings	|	full siblings	|
|B14/B13|	full siblings	|	full siblings	|
|A14/B23| 	unrelated 		|	full siblings	|
|A7/A9	|	unrelated		|	sample duplication|
|B20/A18|	unrelated		|	full siblings	|
|A58/A59|	full siblings	|	identical twins	|


Look at the missingness of each snp and individual to see which sample in a related pair is "better" to remove.
```bash
qsub -S /bin/bash -N missingness -cwd -q highmem.q -l h_vmem=5G -j y << EOF
module load plink/1.90b
plink --bfile SplX_converted_noHH --missing --out SplX_converted_noHH_miss
EOF
```

Figured out which of the sibling paired samples to remove by using `grep "$Sample" SplX_converted_noHH_miss.imiss` to look at the missing genotype calls for each individual. I prioritized samples that were B as well as those with all datasets (DNA methylation and expression) in addition to lower missingness.

3	A3
56  B11
34  A41          
76  B33  
58  B14
14  A14
86	A18
92	A58
7	A7


I saved these FID and IIDs to a file called `Het_related_sample_remove.txt`

Removing these samples with the following command as well as only selecting autosomal snps.

```bash
qsub -S /bin/bash -N filter_ind -cwd -q highmem.q -l h_vmem=5G -j y << EOF
module load plink/1.90b
plink --bfile SplX_converted_noHH --mind 0.1 --remove Het_related_sample_remove.txt --autosome --make-bed --out SplX_converted_noHH_indfilt
EOF
 # 1568630 variants and 103 people pass filters and QC.
 # Among remaining phenotypes, 50 are cases and 53 are controls.
```

This also removed any individuals with more than 10% missing allele calls (no individuals)
Now let's filter the SNPs based on:

1. missing individual rates < .1 (done)
2. missing genotypes for a SNP < .1  
3. HWE equilibrium pvalue < 0.000001
4. minor allele frequency of > .1

```bash
qsub -S /bin/bash -N Filter_snp -cwd -l h_vmem=10G -j y << EOF
module load plink/1.90b
plink --bfile SplX_converted_noHH_indfilt --maf 0.1 --mind 0.1 --geno 0.1 --hwe 0.000001 'midp' --make-bed --out SplX_converted_noHH_indfilt_snpfilt
EOF

 # 103 people (59 males, 44 females) loaded from .fam.
 # --hwe: 6 variants removed due to Hardy-Weinberg exact test. (no error messages)
 # 1065845 variants removed due to minor allele threshold(s)
 # 502779 variants and 103 people pass filters and QC.
 # Among remaining phenotypes, 50 are cases and 53 are controls.

```

Create a histogram of minor allele frequencies

```bash
qsub -S /bin/bash -N get_maf_freq -cwd -l h_vmem=10G -j y << EOF
module load plink/1.90b
plink --bfile SplX_converted_noHH_indfilt_snpfilt --freq --out SplX_converted_noHH_indfilt_snpfilt.freq
EOF

module load R/3.4.0/gcc.4.7.4
Rscript plot_hist_maf.R SplX_converted_noHH_indfilt_snpfilt.freq.frq 20 SplX_converted_noHH_indfilt_snpfilt.frq.pdf
```
Here is the [histogram](Figures/Allele_freq_all.pdf) output.

I'll also attempt with a MAF of .05 since the sample size has increased.
```bash
qsub -S /bin/bash -N Filter_snp_2 -cwd -l h_vmem=10G -j y << EOF
module load plink/1.90b
plink --bfile SplX_converted_noHH_indfilt --maf 0.05 --mind 0.1 --geno 0.1 --hwe 0.000001 'midp' --make-bed --out SplX_converted_noHH_indfilt_snpfilt_geno05
EOF

 # --hwe: 6 variants removed due to Hardy-Weinberg exact test. (no error messages)
 # 915812 variants removed due to minor allele threshold(s)
 # 652812 variants and 103 people pass filters and QC.
```

9. Observe the samples in regards to 1000G data

Now that the variants have been filtered. Now is a good time to observe which populations our study draws from in relation to 1000G genotype data from multiple global populations.

```bash
 # will create a tab version of the ped
qsub -S /bin/bash -N Plink_recode_eigen -cwd -l h_vmem=10G -j y << EOF
module load plink/1.90b
plink --bfile SplX_converted_noHH_indfilt_snpfilt --recode tabx --out SplX_converted_noHH_indfilt_snpfilt_eigen
EOF

 # eigenstrat cannot handle SNP names larger than 39 characters
 # This script cuts all SNP names to 39 or less characters.
qsub -S /bin/bash -N Plink_recode_eigen_map -cwd -l h_vmem=10G -j y << EOF
module load R/3.4.0/gcc.4.7.4
Rscript Eigen_map.R SplX_converted_noHH_indfilt_snpfilt_eigen.map
EOF

 # This is the textfile used to call the eigen convert command.
vi EigenConvertPedtoStrat.txt
 #~~~~~~~ file for eigen command ~~~~~~~#
genotypename:    SplX_converted_noHH_indfilt_snpfilt_eigen.ped
snpname:         SplX_converted_noHH_indfilt_snpfilt_eigen.map
indivname:       SplX_converted_noHH_indfilt_snpfilt_eigen.ped
outputformat:    EIGENSTRAT
genotypeoutname: SplX_converted_noHH_indfilt_snpfilt_eigen.eigenstratgeno
snpoutname:      SplX_converted_noHH_indfilt_snpfilt_eigen.snp
indivoutname:    SplX_converted_noHH_indfilt_snpfilt_eigen.ind
familynames:     YES
 #~~~~~~~ End file ~~~~~~~#
qsub -S /bin/bash -N Eigen_convert -q highmem.q -cwd -l h_vmem=10G -j y << EOF
module load eigensoft/6.0.1
convertf -p EigenConvertPedtoStrat.txt
EOF

 # merge my data with 1000 genome eigenstrat to generate pca map

vi eigein_mergit.txt
 # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
geno1: ../../../../Fabien_Andrew/plink_1000_genome.eigenstratgeno
snp1:  ../../../../Fabien_Andrew/plink_1000_genome.snp
ind1:  ../../../../Fabien_Andrew/plink_1000_genome.ind
geno2: SplX_converted_noHH_indfilt_snpfilt_eigen.eigenstratgeno
snp2:  SplX_converted_noHH_indfilt_snpfilt_eigen.snp
ind2:  SplX_converted_noHH_indfilt_snpfilt_eigen.ind
genooutfilename:   merge_1000_deepa.eigenstratgeno
snpoutfilename:    merge_1000_deepa.snp
indoutfilename:    merge_1000_deepa.ind
 # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

qsub -S /bin/bash -N Eigen_merge -cwd -l h_vmem=20G -j y << EOF
module load eigensoft/6.0.1
mergeit -p eigein_mergit.txt
EOF

 # run the eigenstrat PCA
qsub -S /bin/bash -N Eigen_pca -cwd -l h_vmem=20G -j y << EOF
module load eigensoft/6.0.1
smartpca.perl -i merge_1000_deepa.eigenstratgeno  -a merge_1000_deepa.snp -b merge_1000_deepa.ind -o merge_1000_deepa.pca -p merge_1000_deepa.plot -l merge_1000_deepa.log -e merge_1000_deepa.eval -m 0
EOF
```




<a name="eQTLfilt"/>

10. Create VCF file for eQTL analysis

Based on the genotype matching to bam files PLUS the bam file QC check. Let's filter the samples, then SNPs, and finally create a VCF for eQTLs

Now let's finally write out the variants to generate the vcf files:
There's no RNAseq data for sample A8 and A62, so I need to remove them. Additionally, samples A11, A48, A52, B3, B4, and B17 did not match up with any RNAseq dataset or the RNAseq library was of poor quality (high duplication, low % protein coding). Additionally, the following bam files were of poor quality despite genotype matching: A19 and B25. Related individuals also need to be removed

```bash
vi eQTL_remove.txt
 ##################
3	A3
7	A7
8	A8
11	A11
14  A14
86	A18
87	A19
34	A41
40	A47
45	A52
92	A58
96	A62
50	B3  
51	B4
56  B11
58	B14      
61	B17
69	B25
76  B33
102	B48
106	B52
 ##################


 # make a new directory for the file clutter
mkdir eQTL_filt

 #filter out these individuals as well as only take autosomes for analysis
qsub -S /bin/bash -N filter_ind_eQTL -cwd -l h_vmem=10G -j y << EOF
module load plink/1.90b
plink --bfile SplX_converted_noHH --remove eQTL_remove.txt --autosome --make-bed --out eQTL_filt/SplX_converted_noHH_eQTL
EOF
 # 1568630 variants and 91 people pass filters and QC.
 # Among remaining phenotypes, 44 are cases and 47 are controls.
cd eQTL_filt/
```

Now we filter the snps before performing some QC again.

1. missing individual rates < .1 (done)
2. missing genotypes for a SNP < .1  
3. HWE equilibrium pvalue < 0.000001
4. minor allele frequency of > .1

```bash
qsub -S /bin/bash -N Filter_snp -cwd -l h_vmem=10G -j y << EOF
module load plink/1.90b
plink --bfile SplX_converted_noHH_eQTL --maf 0.1 --mind 0.1 --geno 0.1 --hwe 0.000001 'midp' --make-bed --out SplX_converted_noHH_eQTL_snpfilt
EOF
 # 91 people (53 males, 38 females) loaded from .fam. # male heavy
 # 0 people removed due to missing genotype data (--mind).
 # --hwe: 4 variants removed due to Hardy-Weinberg exact test.
 # I did get the "hwe observation counts vary by more than 10%." warning
 # so I may have to split up the subpopulation and call hwe
 # the maf filter should overcome this
 # 1071278 variants removed due to minor allele threshold(s)
 # 497348 variants and 91 people pass filters and QC.
```

Next we want to perform a little QC to double check that the samples look good for using in QTL analysis. We look at homozygosity, relatedness, and population stratification.

First, we look at the homozygosity among the samples in the eQTL samples

```bash
qsub -S /bin/bash -N checkHH -cwd -l h_vmem=5G -j y << EOF
module load plink/1.90b
module load R/3.4.0/gcc.4.7.4
plink --bfile SplX_converted_noHH_eQTL_snpfilt --maf 0.1 --indep-pairwise 50 5 0.1 --out SplX_converted_noHH_eQTL_indep
plink --bfile SplX_converted_noHH_eQTL_snpfilt --extract SplX_converted_noHH_eQTL_indep.prune.in --genome -het --out SplX_converted_noHH_eQTL_indep_genome
Rscript ../PlotHeterozygosity.R SplX_converted_noHH_eQTL_indep_genome.het
EOF

 # how many samples show too much homozygosity? 7 samples
awk '$6 >.1 || $6 < -.1' SplX_converted_noHH_eQTL_indep_genome.het | more
 # FID  IID       O(HOM)       E(HOM)        N(NM)            F
 #  19  A23        37964     3.34e+04        58616        0.181
 #   2   A2        36746     3.34e+04        58621       0.1326
 #  60  B16        36140    3.333e+04        58501       0.1115
 #  26  A30        38137     3.34e+04        58618       0.1878
 #  41  A48        36317     3.34e+04        58622       0.1155
 #   1   A1        36990    3.339e+04        58598       0.1428

Rscript ../plot_relatedness.R SplX_converted_noHH_eQTL_indep_genome.genome
```
Check out the relatedness plot [here](Figures/Relatedness_eQTL.pdf). There are no longer related individuals.


Relatedness can also be observed through MDS plots
```bash
qsub -S /bin/bash -N MDS_plot -cwd -l h_vmem=5G -j y << EOF
module load plink/1.90b
module load R/3.4.0/gcc.4.7.4
plink --bfile SplX_converted_noHH_eQTL_snpfilt --read-genome SplX_converted_noHH_eQTL_indep_genome.genome --cluster --mds-plot 2 --out SplX_converted_noHH_eQTL_indep_genome_mds
EOF

 # create a plot with this script
Rscript MDS_eQTL_plot.R
```
I colored the MDS plots by [high homozygosity](Figures/MDS_plot_col_by_homozygosity_eQTL.pdf).

Additionally, one can look at identical by state (IBS) test results to observe relatedness and possible stratification within the cases/controls.
```bash
qsub -S /bin/bash -N ibs_test -cwd -l h_vmem=5G -j y << EOF
module load plink/1.90b
plink --bfile SplX_converted_noHH_eQTL_snpfilt --indep-pairwise 50 5 0.1 --ibs-test
EOF
```

This is the output. No significant differences between case and control when compared to control/control differences and case/case differences.

```
Between-group IBS (mean, SD)   = 0.694755, 0.0158478
In-group (case) IBS (mean, SD) = 0.695376, 0.014735
In-group (ctrl) IBS (mean, SD) = 0.693721, 0.0172639
Approximate proportion of variance between group = 6.71662e-05
IBS group-difference empirical p-values:
     T1: Case/control less similar                p = 0.648964
     T2: Case/control more similar                p = 0.351036

   T3: Case/case less similar than ctrl/ctrl    p = 0.886411
   T4: Case/case more similar than ctrl/ctrl    p = 0.113589

   T5: Case/case less similar                   p = 0.838852
   T6: Case/case more similar                   p = 0.161148

   T7: Control/control less similar             p = 0.0936391
   T8: Control/control more similar             p = 0.906361

   T9: Case/case less similar than case/ctrl    p = 0.735463
  T10: Case/case more similar than case/ctrl    p = 0.264537

  T11: Ctrl/ctrl less similar than case/ctrl    p = 0.0870891
  T12: Ctrl/ctrl more similar than case/ctrl    p = 0.912911
```

Create the VCF
```bash
qsub -S /bin/bash -N plink2vcf_eQTL -cwd -q highmem.q -l h_vmem=10G -j y << EOF
module load plink/1.90b
plink --bfile SplX_converted_noHH_eQTL_snpfilt --recode vcf --out SplX_converted_noHH_eQTL_snpfilt_hg19
EOF

 # Warning: At least one VCF allele code violates the official specification;
 # other tools may not accept the file.  (Valid codes must either start with a
 # '<', only contain characters in {A,C,G,T,N,a,c,g,t,n}, be an isolated '*', or
 # represent a breakend.)
```

The genotype data is in hg19 but the gene expression and DNA methylation data is in hg38. Therefore, the vcf file will have to be lifted over. I got the vcf lift over script from https://github.com/knmkr/lift-over-vcf. It takes forever and is single threaded. You can split up the vcf, run and merge again.
```bash
 # bgzip and tabix the vcf
mkdir lift_vcf # for file clutter
module load htslib/1.2.1/gcc.4.4.7
bgzip -c SplX_converted_noHH_eQTL_snpfilt_hg19.vcf > lift_vcf/SplX_converted_noHH_eQTL_snpfilt_hg19.vcf.gz

 # vcf files will need headers
awk 'NR<29' SplX_converted_noHH_eQTL_snpfilt_hg19.vcf > lift_vcf/SplX_converted_noHH_eQTL_snpfilt_hg19.header.txt

 # manually changed 1 digit samples to 2
 # example: A1 to A01

cd lift_vcf
tabix -p vcf SplX_converted_noHH_eQTL_snpfilt_hg19.vcf.gz

 # split vcf and liftover each chromosome seperately
for chr in {1..22};
do
qsub -S /bin/bash -N lift_${chr} -cwd -l h_vmem=4G -j y << EOF
module load htslib/1.2.1/gcc.4.4.7
tabix SplX_converted_noHH_eQTL_snpfilt_hg19.vcf.gz ${chr} > chr${chr}.vcf
module load python/2.7.8/gcc.4.4.7
module load vcftools/11172016/gcc.4.4.7
cat chr${chr}.vcf | python ~/Programs/vcf-liftOver/lift_over.py --chain hg19ToHg38 > chr${chr}_hg38.vcf
EOF
done

 # create file for concat'ing the lifted vcf files
for chr in {1..22};
do
echo chr${chr}_hg38.vcf >> lifted_vcfs.txt
done

 # place the header on each
for chr in {1..22};
do
cat SplX_converted_noHH_eQTL_snpfilt_hg19.header.txt chr${chr}_hg38.vcf > chr${chr}_hg38_head.vcf
mv chr${chr}_hg38_head.vcf chr${chr}_hg38.vcf
done

 # concat the files
module load vcftools/11172016/gcc.4.4.7
vcf-concat -f lifted_vcfs.txt > SplX_converted_noHH_eQTL_snpfilt_hg38.vcf

wc -l SplX_converted_noHH_eQTL_snpfilt_hg38.vcf
 # 497,064 SplX_converted_noHH_eQTL_snpfilt_hg38.vcf

 # sort the vcf
module load vcftools/11172016/gcc.4.4.7
vcf-sort SplX_converted_noHH_eQTL_snpfilt_hg38.vcf > SplX_converted_noHH_eQTL_snpfilt_hg38_sort.vcf
```

Then I finished compressing and indexing the vcf for use in QTLtools.

```bash
module load htslib/1.2.1/gcc.4.4.7
bgzip SplX_converted_noHH_eQTL_snpfilt_hg38_sort.vcf
tabix -p vcf SplX_converted_noHH_eQTL_snpfilt_hg38_sort.vcf.gz
 # check to see if the tbi is readable by QTLtools
module load bcftools/1.3.1/gcc.4.4.7
bcftools view SplX_converted_noHH_eQTL_snpfilt_hg38_sort.vcf.gz | less -S

cp SplX_converted_noHH_eQTL_snpfilt_hg38_sort.vcf* ../../../../../eQTL/Genotype2/.
```
With the genotype file for eQTL analysis completed, the rest of the eQTL analysis can be found in the [eQTL pipeline markdown](../eQTL/eQTL_pipeline2.md)


<a name="mQTLfilt">

11. Create VCF file for mQTL analysis

Based on the genotype matching to the sorted HELP-tagging bam files. Let's filter the samples, then SNPs, and finally create a VCF for mQTLs. See the `sample_sequencing_info.xlsx` file.

```bash
vi mQTL_remove.txt
 ##################
3	A3
7	A7
14	A14
86	A18
88	A22
19	A23
20	A24
32	A38
34	A41
40	A47
44	A51
45	A52
92	A58
93	A59
56  B11
58	B14     
64	B20
65	B21
66	B22
68	B24
69	B25
76  B33
77	B34
82	B39
83	B40
85	B42
97	B43
99	B45
100	B46
101	B47
102	B48
103	B49
105	B51
106	B52
107	B53
 ##################

 # make a new directory for the file clutter
mkdir mQTL_filt

 #filter out these individuals as well as only take autosomes for analysis
qsub -S /bin/bash -N filter_ind_mQTL -cwd -l h_vmem=10G -j y << EOF
module load plink/1.90b
plink --bfile SplX_converted_noHH --remove mQTL_remove.txt --autosome --make-bed --out mQTL_filt/SplX_converted_noHH_mQTL
EOF
 # 1568630 variants and 77 people pass filters and QC.
 # Among remaining phenotypes, 32 are cases and 45 are controls.
cd mQTL_filt/
```

Now we filter the snps before performing some QC again.

1. missing individual rates < .1 (done)
2. missing genotypes for a SNP < .1  
3. HWE equilibrium pvalue < 0.000001
4. minor allele frequency of > .1

```bash
qsub -S /bin/bash -N Filter_snp -cwd -l h_vmem=10G -j y << EOF
module load plink/1.90b
plink --bfile SplX_converted_noHH_mQTL --maf 0.1 --mind 0.1 --geno 0.1 --hwe 0.000001 'midp' --make-bed --out SplX_converted_noHH_mQTL_snpfilt
EOF
 # 77 people (43 males, 34 females) loaded from .fam. # male heavy
 # 0 people removed due to missing genotype data (--mind).
 # --hwe: 1 variants removed due to Hardy-Weinberg exact test.
 # 1069034 variants removed due to minor allele threshold(s)
 # 499595 variants and 77 people pass filters and QC.
```

Next we want to perform a little QC to double check that the samples look good for using in QTL analysis. We look at homozygosity, relatedness, and population stratification.

First, we look at the homozygosity among the samples in the mQTL samples

```bash
qsub -S /bin/bash -N checkHH -cwd -l h_vmem=5G -j y << EOF
module load plink/1.90b
module load R/3.4.0/gcc.4.7.4
plink --bfile SplX_converted_noHH_mQTL_snpfilt --maf 0.1 --indep-pairwise 50 5 0.1 --out SplX_converted_noHH_mQTL_indep
plink --bfile SplX_converted_noHH_mQTL_snpfilt --extract SplX_converted_noHH_mQTL_indep.prune.in --genome -het --out SplX_converted_noHH_mQTL_indep_genome
Rscript ../PlotHeterozygosity.R SplX_converted_noHH_mQTL_indep_genome.het
EOF

 # how many samples show too much homozygosity? 7 samples
awk '$6 >.1 || $6 < -.1' SplX_converted_noHH_mQTL_indep_genome.het | more
 # FID  IID       O(HOM)       E(HOM)        N(NM)            F
 #   2   A2        34378    3.125e+04        55272       0.1304
 #  61  B17        34100    3.124e+04        55265        0.119
 #  60  B16        33820    3.118e+04        55155       0.1101
 #  26  A30        35553    3.125e+04        55274       0.1792
 #  41  A48        33894    3.124e+04        55267       0.1103
 #   1   A1        34520    3.124e+04        55254       0.1367
Rscript ../plot_relatedness.R SplX_converted_noHH_mQTL_indep_genome.genome
```
Check out the replatedness plot [here](Figures/Relatedness_mQTL.pdf). There are no longer related individuals.


Relatedness can also be observed through MDS plots
```bash
qsub -S /bin/bash -N MDS_plot -cwd -l h_vmem=5G -j y << EOF
module load plink/1.90b
module load R/3.4.0/gcc.4.7.4
plink --bfile SplX_converted_noHH_mQTL_snpfilt --read-genome SplX_converted_noHH_mQTL_indep_genome.genome --cluster --mds-plot 2 --out SplX_converted_noHH_mQTL_indep_genome_mds
EOF

 # create a plot with this script
Rscript MDS_mQTL_plot.R
```
I colored the MDS plots by [high homozygosity](Figures/MDS_plot_col_by_homozygosity_mQTL.pdf).
The plot clearly shows that the high homozygosity individuals are on the "asian-european" axis


Additionally, one can look at identical by state (IBS) test results to observe relatedness and possible stratification within the cases/controls.
```bash
qsub -S /bin/bash -N ibs_test -cwd -l h_vmem=5G -j y << EOF
module load plink/1.90b
plink --bfile SplX_converted_noHH_mQTL_snpfilt --indep-pairwise 50 5 0.1 --ibs-test
EOF
```

This is the output. No significant differences between case and control when compared to control/control differences and case/case differences.

```
 Between-group IBS (mean, SD)   = 0.695049, 0.0160264
 In-group (case) IBS (mean, SD) = 0.695231, 0.015462
 In-group (ctrl) IBS (mean, SD) = 0.694243, 0.0164481
 Approximate proportion of variance between group = 0.000219488
 IBS group-difference empirical p-values:
   T1: Case/control less similar                p = 0.824072
   T2: Case/control more similar                p = 0.175928

   T3: Case/case less similar than ctrl/ctrl    p = 0.763292
   T4: Case/case more similar than ctrl/ctrl    p = 0.236708

   T5: Case/case less similar                   p = 0.683483
   T6: Case/case more similar                   p = 0.316517

   T7: Control/control less similar             p = 0.200518
   T8: Control/control more similar             p = 0.799482

   T9: Case/case less similar than case/ctrl    p = 0.525795
   T10: Case/case more similar than case/ctrl    p = 0.474205

   T11: Ctrl/ctrl less similar than case/ctrl    p = 0.178288
   T12: Ctrl/ctrl more similar than case/ctrl    p = 0.821712

```

Create the VCF
```bash
qsub -S /bin/bash -N plink2vcf_mQTL -cwd -q highmem.q -l h_vmem=10G -j y << EOF
module load plink/1.90b
plink --bfile SplX_converted_noHH_mQTL_snpfilt --recode vcf --out SplX_converted_noHH_mQTL_snpfilt_hg19
EOF

 # Warning: At least one VCF allele code violates the official specification;
 # other tools may not accept the file.  (Valid codes must either start with a
 # '<', only contain characters in {A,C,G,T,N,a,c,g,t,n}, be an isolated '*', or
 # represent a breakend.)
```

The genotype data is in hg19 but the gene expression and DNA methylation data is in hg38. Therefore, the vcf file will have to be lifted over. I got the vcf lift over script from https://github.com/knmkr/lift-over-vcf. It takes forever and is single threaded. You can split up the vcf, run and merge again.

```bash
 # bgzip and tabix the vcf
mkdir lift_vcf # for file clutter
module load htslib/1.2.1/gcc.4.4.7
bgzip -c SplX_converted_noHH_mQTL_snpfilt_hg19.vcf > lift_vcf/SplX_converted_noHH_mQTL_snpfilt_hg19.vcf.gz

 # vcf files will need headers
awk 'NR<29' SplX_converted_noHH_mQTL_snpfilt_hg19.vcf > lift_vcf/SplX_converted_noHH_mQTL_snpfilt_hg19.header.txt

 # manually changed 1 digit samples to 2
 # example: A1 to A01

cd lift_vcf
tabix -p vcf SplX_converted_noHH_mQTL_snpfilt_hg19.vcf.gz

 # split vcf and liftover each chromosome seperately
for chr in {1..22};
do
qsub -S /bin/bash -N lift_${chr} -cwd -l h_vmem=4G -j y << EOF
module load htslib/1.2.1/gcc.4.4.7
tabix SplX_converted_noHH_mQTL_snpfilt_hg19.vcf.gz ${chr} > chr${chr}.vcf
module load python/2.7.8/gcc.4.4.7
module load vcftools/11172016/gcc.4.4.7
cat chr${chr}.vcf | python ~/Programs/vcf-liftOver/lift_over.py --chain hg19ToHg38 > chr${chr}_hg38.vcf
EOF
done

 # create file for concat'ing the lifted vcf files
for chr in {1..22};
do
echo chr${chr}_hg38.vcf >> lifted_vcfs.txt
done

 # place the header on each
for chr in {1..22};
do
cat SplX_converted_noHH_mQTL_snpfilt_hg19.header.txt chr${chr}_hg38.vcf > chr${chr}_hg38_head.vcf
mv chr${chr}_hg38_head.vcf chr${chr}_hg38.vcf
done

 # concat the files
module load vcftools/11172016/gcc.4.4.7
vcf-concat -f lifted_vcfs.txt > SplX_converted_noHH_mQTL_snpfilt_hg38.vcf

wc -l SplX_converted_noHH_mQTL_snpfilt_hg38.vcf
 # 499,312 SplX_converted_noHH_mQTL_snpfilt_hg38.vcf

 # sort the vcf
module load vcftools/11172016/gcc.4.4.7
vcf-sort SplX_converted_noHH_mQTL_snpfilt_hg38.vcf > SplX_converted_noHH_mQTL_snpfilt_hg38_sort.vcf
```

Then I finished compressing and indexing the vcf for use in QTLtools.

```bash
module load htslib/1.2.1/gcc.4.4.7
bgzip SplX_converted_noHH_mQTL_snpfilt_hg38_sort.vcf
tabix -p vcf SplX_converted_noHH_mQTL_snpfilt_hg38_sort.vcf.gz
 # check to see if the tbi is readable by QTLtools
module load bcftools/1.3.1/gcc.4.4.7
bcftools view SplX_converted_noHH_mQTL_snpfilt_hg38_sort.vcf.gz | less -S

cp SplX_converted_noHH_mQTL_snpfilt_hg38_sort.vcf* ../../../../../mQTL/Genotype2/.
```
With the genotype file for mQTL analysis completed, the rest of the mQTL analysis can be found in the [mQTL pipeline markdown](../mQTL/mQTL_pipeline2.md)
