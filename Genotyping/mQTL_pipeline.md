# mQTL pipeline for Asthma Project
## Andrew Johnston
## 08/22/17					
## Workflow to perform mQTL analysis

Note: I will be using the vcf file before lifting over coordinates from hg19 to hg38 (so will use hg19 coords)

1. [Download files](#download)
2. [Genotype data](#geno)
3. [Matching methylation and genotype data](#match)
4. [Quantify methylation data](#quant)
5. [PCA](#pca)
6. [Nominal mQTL calling](#nommqtl)
7. [Permutation mQTL calling](#permmqtl)
8. [Annotating mQTLs to genes](#mqtlgene)
9. [Plotting mQTLs](#eigen)
10. [Overlap with Blueprint study](#otherstudy)
11. [mQTL Filtering](#mQTLfilt)


<a name="download"/>

1. Download bam files
I need to download the methylation bam files so that I can perform QTLtools match mode and make sure that all of the data is properly paired with the genotyping (and by transitive property, RNAseq data). They're found in these jobs:
1. J575
2. J341
3. J72
4. J71
5. J25


```
cd /home/greally-lab/Deepa_Andrew/mQTL/Bams/raw_bams/
 #J575
wget http://waspsystem.einstein.yu.edu/wsfile/get/file/cdcc55cd-0fdf-4b13-942a-189a2a138c55,7f66c734-62ce-4279-915c-dc31b13c9be9,1d0b6850-7dc7-4d72-8ae0-d10fa9ca15a0,f1b118bb-ea86-43d0-a2fa-cc3568970c14,253c63d3-8fb3-4a68-87ca-29c738ad670b,f4604bb9-10ad-43ec-9e32-c5857eba6a40,702929d7-41aa-4f3b-aea4-622de931466c,c092b745-4a34-4383-a7be-319fd61e0744,c1d4aa48-4f55-4e57-8492-95daa2352327,b3926083-cdc3-4b8b-9d1b-eca36cbce8f9,da28cee4-9698-449b-ad8d-59b66b56a8ab,e02995b6-94c5-4566-aad0-163b6289ee14,65566d2f-0e83-4242-b652-c4eb9a406a9c,9b4c38c9-628b-4a0f-96d8-75117cb51527,321a6527-2a5a-4088-8d75-f8d74b2caca9,c31e8f24-473a-4e30-8f67-e327dd018861,837c70a2-43fb-4577-9e22-c62919abe65e,8232e54c-d984-4a0b-8b45-5cdd57ff683b,6cb1b5e8-007d-43a9-ab06-104f6a162178,60a0e623-7f34-40e1-be81-a9afb380446e,fb627a90-e452-41a0-94c1-f9f9b5dd5aca,ec88291b-3047-475c-a80b-a2ddd05c25b4,77a0a215-3257-4298-b1f8-8b2eaec4c1b4,eafb0ce7-d737-4d74-9c64-2d35ee7f9343,be411b88-ac71-49fe-9d85-7ad5a38c448e,37f387ae-c451-4b7c-a2a0-d23169be8177,40164d61-cd1e-4db4-aa9f-38be3ed8c0e5,b409f907-1df2-4041-b6b9-1d431c12f7d5,1a46da77-b22d-41f4-89a9-fe36e65bcf84,22f8416b-6a32-4bc7-8549-56dfefb61f73,a2551bbc-ad43-456d-a187-b12f015365c2,80b63255-cd45-4977-b51f-e5929ee02b0e,f785f5f0-2086-4b8e-905d-02a16cc3b23b,ff762b9a-b0d0-418e-8125-9c601dd188ed,a981a05d-e231-46d2-b75b-d698063a8887,aa951f88-a4bf-42a9-9c98-d57dfbe125d6,8971497f-1354-4f48-8c4b-4f04be9d4ae4,169feeb3-fb1b-4ac0-ac71-2425aa97a0e8,cec2deca-3bb9-4f0b-8f5e-6ea4ae345f4a,2bb99bca-02a9-42f3-a54b-ba6d4c3adddd,ae635ff4-e0b1-4d72-978b-315047656c05,655d6e9b-5a8e-4705-aa2f-c16bd19d3c3f,5c4ac0d5-1a66-4d3f-8fea-d50844a8e380,174d2ca7-4f7f-435d-9151-2cd8d3cd2093,6a1207b9-5b82-447b-b7cf-355a7d6ba8f1,fe69f87a-cf64-4f48-a2a9-64eb53200ba5/multi_files.tar
tar xf multi_files.tar
 #J341
wget http://waspsystem.einstein.yu.edu/wsfile/get/file/2453aec1-2417-4404-9ff6-2723d1ae3e6c,828d9149-a789-4729-9256-d11ac6a76c64,16cb913f-1255-4e04-aa8c-77aa1aedd414,ecc9e887-7b67-4e0f-bb24-961bc447e4b0,151c8172-f3c0-4f33-b094-b72dab530be7,83a53f8a-4da3-46af-8b54-f9fc028a0f3a,50999a5f-c2f7-4793-a032-87f0624ba617,30f739a5-5d3c-45e9-bc68-2bf6550fb29a,bd3452f2-5b5b-46cc-b514-13b999634536,b5a08fb7-da92-4a46-87eb-00772b7a160c,85826ab6-b107-4484-ab6c-41a07a6da898,a836ab3d-ceaf-416e-a6b9-1ed8155ac732,dc2a4f5d-9903-47fb-ab72-5f49d9c5fdd5,eacdf100-c40a-4ed4-bcd6-1d81cd7000aa,c522274f-fd4b-41a2-91f0-7cda7eb74902,5173b95f-4b48-45d3-ad89-e1d392e22567,cb8f27c4-7769-4e16-a62f-8a82c2d81110,a7215e6a-fc5b-4d3f-874c-9bb1a58cca8c,ea7a3885-3c92-4f1a-bbac-af4b6e4bec33,3a576eb2-e11e-4578-8a1e-5117ff88a73e,6aed41a3-2fe6-442c-983c-62639ed756b0,530bacd0-ed6d-452e-ae47-06809ee72531,a8f05fc4-fae1-42d9-8021-2b44ced6ec52,a7c46a7e-7ccf-48a3-a35c-c0a8e4704fde,ad0c1c50-1eff-45ee-b4af-db65fe1b4efc,263f3731-f488-4757-a55a-e33f42123153,25a528f3-bee2-4d9c-a10d-e564ead2e282,4aa7bae6-35b5-467a-a1aa-5c8fe85c116a,190a9fd9-1a22-4383-a221-0f91cd46905c,0b13e539-250b-4128-abec-04ff0325230d,5a0a0ea8-2462-4316-ac09-87d1fd09fe30,bb0f11f9-559f-4764-a361-8d2e86a12687/multi_files_341.tar
tar xf multi_files_341.tar
 # J72
wget http://waspsystem.einstein.yu.edu/wsfile/get/file/7f7b80f0-8228-480e-a195-0f045b4bd88a/a1.C5TJVACXX.C6_SA_ICGATGT.bwa.bam
wget http://waspsystem.einstein.yu.edu/wsfile/get/file/4eb1a17d-f709-45f9-a2e1-6a4d2095d43a/a10.C5TJVACXX.C3_SA_IGCCAAT.bwa.bam
wget http://waspsystem.einstein.yu.edu/wsfile/get/file/aad3682c-f3f4-4116-960c-54b6d7e9b7b8/a11.C5TJVACXX.C2_SA_IACAGTG.bwa.bam
wget http://waspsystem.einstein.yu.edu/wsfile/get/file/793dc3a0-2670-4b96-b05e-60ee6d2e165f/a12.C5TJVACXX.C2_SA_IGCCAAT.bwa.bam
wget http://waspsystem.einstein.yu.edu/wsfile/get/file/b78aa72b-952b-44f6-b15a-87bfc92a63b5/a13.C5TJVACXX.C5_SA_ICGATGT.bwa.bam
wget http://waspsystem.einstein.yu.edu/wsfile/get/file/944a14c6-b144-4083-9f55-f8858079f53b/a15.C5TJVACXX.C5_SA_ITGACCA.bwa.bam
wget http://waspsystem.einstein.yu.edu/wsfile/get/file/fbe1221a-10c0-412a-963b-46ea9e3ccc73/a16.C5TJVACXX.C5_SA_IACAGTG.bwa.bam
wget http://waspsystem.einstein.yu.edu/wsfile/get/file/3af2b862-a3ed-402d-80df-e16509f396c8/a17.C5TJVACXX.C5_SA_IGCCAAT.bwa.bam
wget http://waspsystem.einstein.yu.edu/wsfile/get/file/14fbae40-3c41-47f5-bbc8-852540127f3a/a25.C5TJVACXX.C6_SA_IGCCAAT.bwa.bam
wget http://waspsystem.einstein.yu.edu/wsfile/get/file/754ac109-01bb-47d7-979b-059f6a8d8108/a19.C5TJVACXX.C6_SA_IACAGTG.bwa.bam
wget http://waspsystem.einstein.yu.edu/wsfile/get/file/68721134-4f22-40a6-bd70-a2be548e0b59/a4.C5TJVACXX.C3_SA_ICGATGT.bwa.bam
wget http://waspsystem.einstein.yu.edu/wsfile/get/file/91d6efa3-b954-4283-bcac-009546927d31/a5.C5TJVACXX.C3_SA_ITGACCA.bwa.bam
wget http://waspsystem.einstein.yu.edu/wsfile/get/file/d03d4b3a-46e0-48d9-811d-8731f067dd20/a6.C5TJVACXX.C3_SA_IACAGTG.bwa.bam
wget http://waspsystem.einstein.yu.edu/wsfile/get/file/70fbfa2e-4228-477c-91a3-6f1493fefee5/a7.C5TJVACXX.C2_SA_ICGATGT.bwa.bam
wget http://waspsystem.einstein.yu.edu/wsfile/get/file/3867fb17-4a18-4be5-873c-8e39286d58b7/a8.C5TJVACXX.C6_SA_ITGACCA.bwa.bam
wget http://waspsystem.einstein.yu.edu/wsfile/get/file/98126e1c-1991-4d29-9876-260c3051a8c4/a9.C5TJVACXX.C2_SA_ITGACCA.bwa.bam
wget http://waspsystem.einstein.yu.edu/wsfile/get/file/2e2307ff-fe14-4148-a5a3-7ed41e10c4d1/b10.C5TJVACXX.C2_SA_IAGTTCC.bwa.bam
wget http://waspsystem.einstein.yu.edu/wsfile/get/file/5481cc25-b858-4bc6-9708-be36ef187d68/b11.C5TJVACXX.C2_SA_IGATCAG.bwa.bam
wget http://waspsystem.einstein.yu.edu/wsfile/get/file/b57d1b1b-ec59-439f-9cb8-9b23474d24e8/b12.C5TJVACXX.C5_SA_ICAGATC.bwa.bam
wget http://waspsystem.einstein.yu.edu/wsfile/get/file/f1cc8687-6bf9-49d2-88ac-c321354fb02d/b13.C5TJVACXX.C5_SA_IAGTTCC.bwa.bam
wget http://waspsystem.einstein.yu.edu/wsfile/get/file/f8650fab-5750-4785-9a95-5bab2a4a052b/b14.C5TJVACXX.C5_SA_IGATCAG.bwa.bam
wget http://waspsystem.einstein.yu.edu/wsfile/get/file/e969049e-8241-401b-98e3-2e5e93923111/b15.C5TJVACXX.C2_SA_ICCGTCC.bwa.bam
wget http://waspsystem.einstein.yu.edu/wsfile/get/file/1b8f5cf7-45aa-4cf3-bb0e-a68c7fcb8ba4/b16.C5TJVACXX.C6_SA_IGATCAG.bwa.bam
wget http://waspsystem.einstein.yu.edu/wsfile/get/file/5aa9f2fa-78c7-4cb3-8ebd-6a8336996d14/b17.C5TJVACXX.C5_SA_ICCGTCC.bwa.bam
wget http://waspsystem.einstein.yu.edu/wsfile/get/file/cd1fa99c-adb6-4079-a226-74c6f7e377c3/b26.C5TJVACXX.C6_SA_ICCGTCC.bwa.bam
wget http://waspsystem.einstein.yu.edu/wsfile/get/file/53edbf33-fa86-4afe-a2be-a9baebaba5ba/b3.C5TJVACXX.C3_SA_ICAGATC.bwa.bam
wget http://waspsystem.einstein.yu.edu/wsfile/get/file/45537f26-6da8-4fff-8ef1-7d43cd86ab0e/b4.C5TJVACXX.C3_SA_IAGTTCC.bwa.bam
wget http://waspsystem.einstein.yu.edu/wsfile/get/file/c5a43437-6723-46fd-8f49-8b3b474b3953/b5.C5TJVACXX.C3_SA_IGATCAG.bwa.bam
wget http://waspsystem.einstein.yu.edu/wsfile/get/file/02b6fcb6-279d-471a-aa8b-f14ef6232b45/b6.C5TJVACXX.C6_SA_ICAGATC.bwa.bam
wget http://waspsystem.einstein.yu.edu/wsfile/get/file/4b2edad3-b89a-452b-acc9-df9f9e585039/b7.C5TJVACXX.C3_SA_ICCGTCC.bwa.bam
wget http://waspsystem.einstein.yu.edu/wsfile/get/file/971ed51e-8674-46a7-bccf-105f3beb85ed/b8.C5TJVACXX.C6_SA_IAGTTCC.bwa.bam
wget http://waspsystem.einstein.yu.edu/wsfile/get/file/ac02abdc-4ef4-4af5-a8db-b67366be9140/b9.C5TJVACXX.C2_SA_ICAGATC.bwa.bam

 # J71
wget http://waspsystem.einstein.yu.edu/wsfile/get/file/166468ba-259c-4bee-8d03-6003bedef64d,ca73df34-d4d4-4335-a23b-7db116ad5ea1,ce197143-bb6a-4a21-af85-f2cb793897ac,ecee1a22-309f-4c58-ae36-d35776f4bbb2,db940605-7871-4293-a533-879f35ffaa47,bb669bc5-78d1-49c5-8816-da50d8db5948,16f91eb7-9914-4a9e-b062-52716bb4382d,98e3af54-8ac3-4bca-a1f7-aea52fd8856e/multi_files_J71.tar
tar xf multi_files_J71.tar
rm multi_files_J71.tar
 # J25
 # wget for the tar folder failed so here they are separately
wget http://waspsystem.einstein.yu.edu/wsfile/get/file/62802592-f3f1-4c4e-a4e6-119402dd193b/a21.C56CUACXX.C2_SA_ICGATGT.bwa.bam
wget http://waspsystem.einstein.yu.edu/wsfile/get/file/253b0d23-568b-4292-91ed-d6cdb8c9810f/a22.C56CUACXX.C2_SA_ITGACCA.bwa.bam
wget http://waspsystem.einstein.yu.edu/wsfile/get/file/22f2d8d8-ae6d-4b5c-a7f2-3a3af774e588/b24.C56CUACXX.C2_SA_ICCGTCC.bwa.bam
wget http://waspsystem.einstein.yu.edu/wsfile/get/file/2441e461-7d38-450e-b17b-c6b5c2141ea6/b23.C56CUACXX.C2_SA_IGATCAG.bwa.bam
wget http://waspsystem.einstein.yu.edu/wsfile/get/file/ac587ed5-01f8-45e4-80f9-1b68d243bbe2/b21.C56CUACXX.C2_SA_IAGTTCC.bwa.bam
wget http://waspsystem.einstein.yu.edu/wsfile/get/file/3b6724db-69e9-49b4-935a-973d57f511a5/a24.C56CUACXX.C2_SA_IGCCAAT.bwa.bam
wget http://waspsystem.einstein.yu.edu/wsfile/get/file/0831bf29-d5d6-45f6-8ce5-9c5a8f214f41/a23.C56CUACXX.C2_SA_IACAGTG.bwa.bam
wget http://waspsystem.einstein.yu.edu/wsfile/get/file/f3e0cedb-8708-4196-9b9a-a4473bca1717/b18.C56CUACXX.C2_SA_ICAGATC.bwa.bam
```

To make downstream analysis easier, I renamed some of the files manually.

```
mv a1.C5TJVACXX.C6_SA_ICGATGT.bwa.bam A01.bam
mv a2.C5TJVACXX.C1_SA_ICAGATC.bwa.bam A02.bam
mv a3.C5TJVACXX.C1_SA_IAGTTCC.bwa.bam A03.bam
mv a4.C5TJVACXX.C3_SA_ICGATGT.bwa.bam A04.bam
mv a5.C5TJVACXX.C3_SA_ITGACCA.bwa.bam A05.bam
mv a6.C5TJVACXX.C3_SA_IACAGTG.bwa.bam A06.bam
mv a7.C5TJVACXX.C2_SA_ICGATGT.bwa.bam A07.bam
mv a8.C5TJVACXX.C6_SA_ITGACCA.bwa.bam A08.bam
mv a9.C5TJVACXX.C2_SA_ITGACCA.bwa.bam A09.bam

for f1 in a*.C5TJVACXX*
do
SAMPLE="$(echo ${f1} | cut -c 2-3)"
mv $f1 A${SAMPLE}.bam
echo A${SAMPLE}.bam
done

for f1 in a*.C56CUACXX*
do
SAMPLE="$(echo ${f1} | cut -c 2-3)"
mv $f1 A${SAMPLE}.bam
echo A${SAMPLE}.bam
done

mv b1.C5TJVACXX.C1_SA_IGATCAG.bwa.bam B01.bam
mv b2.C5TJVACXX.C1_SA_ICCGTCC.bwa.bam B02.bam
mv b3.C5TJVACXX.C3_SA_ICAGATC.bwa.bam B03.bam
mv b4.C5TJVACXX.C3_SA_IAGTTCC.bwa.bam B04.bam
mv b5.C5TJVACXX.C3_SA_IGATCAG.bwa.bam B05.bam
mv b6.C5TJVACXX.C6_SA_ICAGATC.bwa.bam B06.bam
mv b7.C5TJVACXX.C3_SA_ICCGTCC.bwa.bam B07.bam
mv b8.C5TJVACXX.C6_SA_IAGTTCC.bwa.bam B08.bam
mv b9.C5TJVACXX.C2_SA_ICAGATC.bwa.bam B09.bam

for f1 in b*.C5TJVACXX*
do
SAMPLE="$(echo ${f1} | cut -c 2-3)"
mv $f1 B${SAMPLE}.bam
echo B${SAMPLE}.bam
done

for f1 in b*.C56CUACXX*
do
SAMPLE="$(echo ${f1} | cut -c 2-3)"
mv $f1 B${SAMPLE}.bam
echo B${SAMPLE}.bam
done
```

I need to move all of the sorted bam files to the same folder. The older Help tagging jobs were sorted bam files. The newer ones had both sorted and unsorted files. I was able to easily move the older ones because I had tried indexing the files before moving them.

```
mv *_sorted.bam ../sorted/

for f1 in *.bai
do
SAMPLE="$(echo ${f1} | cut -d '.' -f1)"
echo $SAMPLE
mv ${SAMPLE}.bam ../sorted/
mv $f1 ../sorted/
done
```

Now I need to index the bam files.
```
for f1 in *.bam;
do
SAMPLE="$(echo ${f1} | cut -c1-3)"
qsub -S /bin/bash -N index_${SAMPLE} -cwd -q highmem.q -l h_vmem=10G -j y << EOF
module load samtools/1.4.1/gcc.4.4.7
samtools index $f1
EOF
done
```

<a name="geno"/>

2. Ready the genotype data.
This was performed in the [MEGA_pipeline3.md](../Genotype\ Filtering/MEGA_pipeline3.md), specifically in the mQTL filtering section.

The finalized VCF file that I'm using is `/home/greally-lab/Deepa_Andrew/mQTL/Genotype2/SplX_converted_noHH_mQTL_snpfilt_hg38_sort.vcf.gz`.

With the finalized file I ran PCA on it.

```bash
qsub -S /bin/bash -N PCA_geno -cwd -l h_vmem=10G -j y << EOF
module load QTLtools/1.1/gcc.4.7.4
QTLtools pca --vcf SplX_converted_noHH_mQTL_snpfilt_hg38_sort.vcf.gz --scale --center --maf 0.1 --distance 50000 --out genotypes.pca
EOF
 # I call maf .1 but the variants were previously reduced to this.
 # --distance 50000 to only consider variant sites separated by at least 50kb
 	# how does QTLtools choose which SNPs within 50kb to use?
```

Next I analyzed the PCA data in R.

<a name="match"/>

3. Match methylation data to genotype data

The major step of matching the HELP-tagging bam files to the genotyped individuals (hg19) was performed in the [MEGA_pipeline3.md](../Genotype\ Filtering/MEGA_pipeline3.md). All of the mismatched bam-sample need to be properly switched in the CpG file and the unused samples (poor quality and related) need to be removed.

NOTE: the HELP-tagging bam files are hg19 while the CpG files are hg38.


<a name="quant">

4. Generate phenotype (DNA methylation) file

The bam files were run through the HELP-tagging pipeline to generate angle values for each of the CpGs. The coordinates for the CpGs were lited over to hg38. A table of all of the raw CpG values can be found here: `/home/greally-lab/Deepa_Andrew/mQTL/Quants/all_meth_nona_rma33b27b53.txt`

First we switched sample IDs appropriately and removed unused samples
```bash
cd /home/greally-lab/Deepa_Andrew/mQTL/Quants
module load R/3.4.0/gcc.4.7.4
```
We prepared the methylation quantification as follows:

first I needed to figure out which CpGs overlapped with repeat masker regions:
```bash
awk 'NR>1 {OFS="\t"; print $2,$3,$3+1,$1}' CpGID_Prom_GB_Enh2_annot.txt > CpG_xtid.bed

module load bedtools2
bedtools intersect -a CpG_xtid.bed -b RepeatMasker_hg38.bed | awk '{print $4'} | uniq > CpGs_over_repeats.txt
wc -l CpG_xtid.bed
 # 2,455,516 CpG_xtid.bed
wc -l CpGs_over_repeats.txt
 # 1,272,032
```


```R
 # set options
options(stringsAsFactors = F)
options(scipen = 999)

 # load in file
library(data.table)
meth <- fread("07-2017-all_meth_nona.txt", header=TRUE)
meth <- as.data.frame(meth)

 # switch samples
 # swap A33 and B27
idx_B27 <- which(colnames(meth) %in% "b27")
idx_A33 <- which(colnames(meth) %in% "a33")
meth_a33 <- meth[,idx_A33]
meth[,idx_A33] <- meth[,idx_B27]
meth[,idx_B27] <- meth_a33

 # B44 to A56
idx_A56 <- which(colnames(meth) %in% "a56")
idx_B44 <- which(colnames(meth) %in% "b44")
meth[,idx_B44] <- meth[,idx_A56]

 # swap A55,A56,A57
 # idx_A55 <- which(colnames(meth) %in% "a55") # A55 doesn't exist because it never had a methylation library produced.. or so we thought?!?!!
idx_A56 <- which(colnames(meth) %in% "a56")
idx_A57 <- which(colnames(meth) %in% "a57")
idx_A58 <- which(colnames(meth) %in% "a58")
idx_A59 <- which(colnames(meth) %in% "a59")
meth$a55 <- meth[,idx_A57]
meth[,idx_A56] <- meth[,idx_A58]
meth[,idx_A57] <- meth[,idx_A59]

 # swap B05 and B04
idx_B04 <- which(colnames(meth) %in% "b04")
idx_B05 <- which(colnames(meth) %in% "b05")
meth_B04<- meth[,idx_B04]
meth[,idx_B04] <- meth[,idx_B05]
meth[,idx_B05] <- meth_B04

 # B13 from B14
idx_B13 <- which(colnames(meth) %in% "b13")
idx_B14 <- which(colnames(meth) %in% "b14")
meth[,idx_B13] <- meth[,idx_B14]

 # remove unused samples
 ## make a's and b's capital
colnames(meth) <- gsub("a","A",colnames(meth))
colnames(meth) <- gsub("b","B",colnames(meth))
metadata <- read.csv("../Bams/PCA/sample_metadata_mQTL_3.csv")
samps_keep <- as.character(metadata$Sample)
idx_keep <- which(colnames(meth) %in% samps_keep)

meth_filt <- meth[,c(1,idx_keep)]
dim(meth_filt) # [1] 1433614      78

write.table(meth_filt, "07-2017-all_meth_nona_mQTL.txt", sep="\t", col.names=TRUE, row.names=F, append=F, quote=F)
saveRDS(meth_filt, "meth_filt.rds")

meth_filt <- readRDS("meth_filt.rds")

 # plot the densities for the samples
library(RColorBrewer)

pdf(file = "Density_plot_b4_filt.pdf", width=6, height=4.5, family="ArialMT")
d2 <- list()
for(i in 1:(ncol(meth_filt)-1)){
  d2[[i]] <- density(meth_filt[,i+1])
}
plot(d2[[1]], ylim=c(0,.8))
colBrew <- colorRampPalette(brewer.pal(11,"Spectral"))(length(d2))
for(i in 2:length(d2)){
  lines(d2[[i]], col=colBrew[i])
}
dev.off()

mspI <- read.table(file = "mspi_hg19.hcount", skip = 1)
dim(mspI) # 1954997 6 - checks out
 # wc -l mspi_hg19.hcount
 # 1954998 mspi_hg19.hcount
head(mspI, 10)

MspI_ok <- mspI[which(mspI$V6 > 4), 1]
head(MspI_ok)
length(MspI_ok) #1,713,004

meth_filt_mspI <-meth_filt[which(meth_filt$X.tid %in% MspI_ok),]
dim(meth_filt_mspI)
 # [1] 1303656      78
head(meth_filt_mspI)

 #remove all HpaII sites which overpalled with RepeatMask
CpGs_over_repeats <- read.table("CpGs_over_repeats.txt")
dim(CpGs_over_repeats)
 # [1] 1272032      1
meth_filt_mspI_repeats <- meth_filt_mspI[-which(meth_filt_mspI$X.tid %in% CpGs_over_repeats[,1]),]
dim(meth_filt_mspI_repeats)
 #763,495 CpGs after all of the filters

 # Filter out sex chromosomes
 ## need to add chr information to CpGs
anno <- fread("CpG_xtid.bed")
head(anno)
merge_meth_meth_filt_mspI_repeats <- merge(x=meth_filt_mspI_repeats, y=anno, by.x="X.tid", by.y="V4")
dim(merge_meth_meth_filt_mspI_repeats)
head(merge_meth_meth_filt_mspI_repeats)
meth_filt_mspI_repeats_chrX <-subset(merge_meth_meth_filt_mspI_repeats,merge_meth_meth_filt_mspI_repeats$V1 != "chrX")
dim(meth_filt_mspI_repeats_chrX)
 #[1] 746,584     81
meth_filt_mspI_repeats_chrXY <-subset(meth_filt_mspI_repeats_chrX,meth_filt_mspI_repeats_chrX$V1 != "chrY") # no chrY to begin with I guess.
dim(meth_filt_mspI_repeats_chrXY)
 # [1] 746584     81
table(meth_filt_mspI_repeats_chrXY$V1)

 # make density plots again
meth_final <- meth_filt_mspI_repeats_chrXY[,-c(79:81)]
saveRDS(meth_final, "meth_final.rds")

d3 <- list()
for(i in 1:(ncol(meth_final)-1)){
  d3[[i]] <- density(meth_final[,i+1])
}

pdf(file = "Density_plot_after_filt.pdf", width=6, height=4.5, family="ArialMT")
plot(d3[[1]], ylim=c(0,.8))
colBrew2 <- colorRampPalette(brewer.pal(11,"Spectral"))(length(d3))
for(i in 2:length(d3)){
  lines(d3[[i]], col=colBrew2[i])
}
dev.off()

idx_zero <- NULL
zero_val <- NULL
for(i in 1:length(d3)){
	dt = data.table(d3[[i]][[1]], val = d3[[i]][[1]])
	setattr(dt, "sorted", "V1")
	x=0
	idx_zero <- c(idx_zero,dt[J(x), .I, roll = "nearest", by = .EACHI][[2]])
	zero_val <- c(zero_val,d3[[i]][[2]][idx_zero[i]])
}

idx_bad <- which(zero_val > .15)
length(idx_bad) # 11

colnames(meth_final)[idx_bad+1]
 # [1] "A53" "A54" "A56" "A57" "A60" "A61" "A62" "B41" "B44" "B50" "A55"
 # I'm going to ignore this for now
 # they are all from batches 10 and higher. It could be due to different hands?


 #cluster the samples
meth_final_t=t(meth_final[,-1])
hc_clust <- hclust(dist(meth_final_t), "ward.D2")
par(mfrow=c(1,1))

pdf(file = "hc_clust_after_filt.pdf", width=12, height=4, family="ArialMT")
plot(hc_clust, main="hc clustering of filtered meth data")
dev.off()

 # create the final bed file:
meth_filt_mspI_repeats_chrXY_final <- cbind(meth_filt_mspI_repeats_chrXY[,c(79:81)],
paste("CpG", meth_filt_mspI_repeats_chrXY[,1],sep="_"),
paste("CpG", meth_filt_mspI_repeats_chrXY[,1],sep="_"),
"+",
meth_filt_mspI_repeats_chrXY[,-c(1,79:81)])

colnames(meth_filt_mspI_repeats_chrXY_final)[1:6] <- c("Chr", "start", "end", "pid", "gid", "strand")

 # Now I need to make sure that I am calling my samples the same name and get
 # them in the same order
geno_pca <- read.table("../Genotype2/genotypes.pca.pca", header=T)
head(geno_pca)
geno_IDs <- substring(colnames(geno_pca)[-1], 2, 8)

idx_match <- match(colnames(meth_filt_mspI_repeats_chrXY_final[,-c(1:6)]), do.call(rbind, strsplit(geno_IDs, "_"))[,2])
same_ids <- which(colnames(meth_filt_mspI_repeats_chrXY_final) %in% do.call(rbind, strsplit(geno_IDs, "_"))[,2])
length(same_ids) #77
colnames(meth_filt_mspI_repeats_chrXY_final)[-same_ids]

idx_match2 <- match(do.call(rbind, strsplit(geno_IDs, "_"))[,2],colnames(meth_filt_mspI_repeats_chrXY_final[,-c(1:6)]))

 # add the correct genotype names to expression matrix
mat4 <- meth_filt_mspI_repeats_chrXY_final
colnames(mat4)[-c(1:6)] <- geno_IDs[idx_match]
order(colnames(mat4)[-c(1:6)])
 # reorder the expression columns to match the order in the genotype pca file.
mat4_lite <- mat4[,-c(1:6)]
mat4_lite <- mat4_lite[,idx_match2]
head(mat4_lite)

 # combine the annotation of genes with the reordered sample columns.
mat5 <- cbind(mat4[,c(1:6)], mat4_lite)
dim(mat5) # [1] 746,584     83
head(mat5,1)

 # write out the expression table
write.table(mat5, "all_meth_nona_mQTL_filt.bed", col.names = T, row.names = F,
            append = F, quote = F, sep = "\t")
mat5 <- fread("all_meth_nona_mQTL_filt.bed", header=T)
mat5 <- as.data.frame(mat5)
sum(rowVars(mat5[,-c(1:6)])==0) # 14239 (this would break it)
idx_zeroVar <- which(rowVars(mat5[,-c(1:6)])==0)
table(rowSums(mat5[idx_zeroVar,-c(1:6)])) # 14238 are all hypermethylated and 1 is not methylated in all samples, these need to be removed.
mat6 <- mat5[-idx_zeroVar,]
write.table(mat6, "all_meth_nona_mQTL_filt_2.bed", col.names = T, row.names = F,
            append = F, quote = F, sep = "\t")
```

Make into proper phenotype file
```bash
 # if needed:
 # vi all_meth_nona_mQTL_filt.bed #  add a "#" in front of the first line aka comment it out
 # need to remove chr prefix
sed 's/^chr\|%$//g' all_meth_nona_mQTL_filt_2.bed > all_meth_nona_mQTL_filt_2_nochr.bed

 # bgzip and tabix
module load htslib/1.2.1/gcc.4.4.7
bgzip all_meth_nona_mQTL_filt_2_nochr.bed && tabix -p bed all_meth_nona_mQTL_filt_2_nochr.bed.gz
```

Finally, perform PCA analysis on the methylation data
```bash
 # methylation PCA
qsub -S /bin/bash -N PCA_meth -cwd -l h_vmem=100G -q highmem.q -j y << EOF
module load QTLtools/1.1/gcc.4.7.4
QTLtools pca --bed all_meth_nona_mQTL_filt_nochr.bed.gz --scale --center --out all_meth_nona_mQTL_filt
EOF
```

<a name="pca">

5. Explore Expression PCA

Using the sample metadata, observe which covariates are contributing the most to the various expression PCs. Visualizing using PCA heatmaps.

meth PC: 1-10
meth PC:1-4,7-10
meth PC:1-4,7-10,12,14
Geno PC: 1,2 OR none


<a name="nommqtl"/>

6. Nominal eQTL passes

I needed to try out different combinations of covariates based on what I saw in the PCA heatmaps to make educated guesses on which PCs to include to optimize the eQTL calling.

First I selected:
Meth PC: 1-10,13
Geno PC: 1,2

```bash
cd /home/greally-lab/Deepa_Andrew/mQTL/Genotype2
awk 'NR<4' genotypes.pca.pca > Cov_geno_2PCs.txt

mkdir ../Nominal
cd ../Nominal
cp ../Quants/all_meth_nona_mQTL_filt.pca .

 # creating cov file for first 10 methylation PCs plus the first 2 geno PCs
awk 'NR>1 && NR<12' all_meth_nona_mQTL_filt.pca | cat /home/greally-lab/Deepa_Andrew/mQTL/Genotype2/Cov_geno_2PCs.txt - > Cov_2PC_geno_10PC_meth.txt

 # grabbing meth PCs 1-4
awk 'NR>1 && NR<6' all_meth_nona_mQTL_filt.pca > Cov_1-4_meth.txt
 # grabbing meth PCs 7-10
awk 'NR>7 && NR<12' all_meth_nona_mQTL_filt.pca > Cov_7-10_meth.txt
cat /home/greally-lab/Deepa_Andrew/mQTL/Genotype2/Cov_geno_2PCs.txt Cov_1-4_meth.txt Cov_7-10_meth.txt > Cov_2PC_geno_1-4_7-10_PC_meth.txt

 # grabbing meth PCs 12 and 14
awk 'NR==13' all_meth_nona_mQTL_filt.pca | cat Cov_2PC_geno_1-4_7-10_PC_meth.txt - > Cov_2PC_geno_1-4_7-10_12_PC_meth.txt

awk 'NR==15' all_meth_nona_mQTL_filt.pca | cat Cov_2PC_geno_1-4_7-10_12_PC_meth.txt - > Cov_2PC_geno_1-4_7-10_12_14_PC_meth.txt

awk 'NR!=2 && NR!=3' Cov_2PC_geno_10PC_meth.txt > Cov_10PC_meth.txt

awk 'NR==7' all_meth_nona_mQTL_filt.pca | cat Cov_2PC_geno_1-4_7-10_12_PC_meth.txt - > Cov_2PC_geno_1-4_6-10_12_PC_meth.txt

awk 'NR==6' all_meth_nona_mQTL_filt.pca | cat Cov_2PC_geno_1-4_7-10_12_PC_meth.txt - > Cov_2PC_geno_1-5_7-10_12_PC_meth.txt

 # the files need to be gzip'd
module load htslib/1.2.1/gcc.4.4.7
bgzip Cov_2PC_geno_10PC_meth.txt
bgzip Cov_10PC_meth.txt
bgzip Cov_2PC_geno_1-4_7-10_PC_meth.txt
bgzip Cov_2PC_geno_1-4_7-10_12_PC_meth.txt
bgzip Cov_2PC_geno_1-4_7-10_12_14_PC_meth.txt
bgzip Cov_2PC_geno_1-4_6-10_12_PC_meth.txt
bgzip Cov_2PC_geno_1-5_7-10_12_PC_meth.txt

```

Now let's run the nominal passes for the various covariate files to find the optimal QTL calling covariate file.

Let's first start with the first two genotype PCs (ancestry) and the first 10 methylation PCs
```bash
 # running the permutation
mkdir geno_2_10_meth
cd geno_2_10_meth
for j in $(seq 1 200);
do
qsub -S /bin/bash -N mQTL_${j}_chunk_nom -cwd -l h_vmem=5.6G -j y << EOF
module load QTLtools/1.1/gcc.4.7.4
QTLtools cis --vcf ../../Genotype2/SplX_converted_noHH_mQTL_snpfilt_hg38_sort.vcf.gz --bed ../../Quants/all_meth_nona_mQTL_filt_nochr.bed.gz --cov ../Cov_2PC_geno_10PC_meth.txt.gz --nominal 0.01  --chunk $j 200 --out chunk_$j\_200.txt --normal --window 1000000
EOF
done

cat chunk*_200.txt > nominals_full.txt
wc -l nominals_full.txt
 # 3,730,904 mQTLs (pval <.01) found in nominals_full.txt

 # top eQTL for gene
awk '$14==1'  nominals_full.txt | wc -l
 # 615,187 top eQTLs for gene; so only 10828 genes had a variant associated with it.
```

Now let's remove the first two geno and keep the first 10 methylation PCs
```bash
 # running the permutation
mkdir no_geno_10_meth
cd no_geno_10_meth
for j in $(seq 1 200);
do
qsub -S /bin/bash -N mQTL_${j}_chunk_nom -cwd -l h_vmem=5.6G -j y << EOF
module load QTLtools/1.1/gcc.4.7.4
QTLtools cis --vcf ../../Genotype2/SplX_converted_noHH_mQTL_snpfilt_hg38_sort.vcf.gz --bed ../../Quants/all_meth_nona_mQTL_filt_nochr.bed.gz --cov ../Cov_10PC_meth.txt.gz --nominal 0.01  --chunk $j 200 --out chunk_$j\_200.txt --normal --window 1000000
EOF
done

cat chunk*_200.txt > nominals_full.txt
wc -l nominals_full.txt
 # 3,234,029 eQTLs (pval <.01) found in nominals_full.txt

 # top eQTL for gene
awk '$14==1'  nominals_full.txt | wc -l
 # 625,292 top eQTLs for gene; so only 10828 genes had a variant associated with it.
```

Adding the ancestral PCs helps calling? Yes for SNPs assocaited but reduces the number of CpGs with an associated variant. I'm going to keep the geno PCs because of the known ancestry difference within the population.


Now let's run the nominal pass by only adding the PCs that are significantly influenced by batch and covariates

```bash
cd ../
mkdir geno_2_14_710_meth
cd geno_2_14_710_meth

 # running the permutation
for j in $(seq 1 200);
do
qsub -S /bin/bash -N mQTL_${j}_chunk_nom -cwd -l h_vmem=5.6G -j y << EOF
module load QTLtools/1.1/gcc.4.7.4
QTLtools cis --vcf ../../Genotype2/SplX_converted_noHH_mQTL_snpfilt_hg38_sort.vcf.gz --bed ../../Quants/all_meth_nona_mQTL_filt_nochr.bed.gz --cov ../Cov_2PC_geno_1-4_7-10_PC_meth.txt.gz --nominal 0.01  --chunk $j 200 --out chunk_$j\_200.txt --normal --window 1000000
EOF
done

cat chunk*_200.txt > nominals_full.txt
wc -l nominals_full.txt
 # 3,753,460 mQTLs (pval <.01) found in nominals_full.txt

 # top eQTL for gene
awk '$14==1'  nominals_full.txt | wc -l
 # 624,677 top eQTLs for gene; so only 10764 genes had a variant associated with it.
```

Now let's run the nominal pass by only adding the PCs that are significantly influenced by batch and covariates, plus a couple more

```bash
cd ../
mkdir geno_2_14_710_12_meth
cd geno_2_14_710_12_meth

 # running the permutation
for j in $(seq 1 200);
do
qsub -S /bin/bash -N mQTL_${j}_chunk_nom -cwd -l h_vmem=5.6G -j y << EOF
module load QTLtools/1.1/gcc.4.7.4
QTLtools cis --vcf ../../Genotype2/SplX_converted_noHH_mQTL_snpfilt_hg38_sort.vcf.gz --bed ../../Quants/all_meth_nona_mQTL_filt_nochr.bed.gz --cov ../Cov_2PC_geno_1-4_7-10_12_PC_meth.txt.gz --nominal 0.01  --chunk $j 200 --out chunk_$j\_200.txt --normal --window 1000000
EOF
done

cat chunk*_200.txt > nominals_full.txt
wc -l nominals_full.txt
 # 3753006 mQTLs (pval <.01) found in nominals_full.txt

 # top mQTL for CpG
awk '$14==1'  nominals_full.txt | wc -l
 # 624,929 top mQTLs for CpG
```

Now Now let's run the nominal pass by only adding the PCs that are significantly influenced by batch and covariates, plus a couple more  

```bash
cd ../
mkdir geno_2_14_710_12_14_meth
cd geno_2_14_710_12_14_meth

 # running the permutation
for j in $(seq 1 200);
do
qsub -S /bin/bash -N mQTL_${j}_chunk_nom -cwd -l h_vmem=5.6G -j y << EOF
module load QTLtools/1.1/gcc.4.7.4
QTLtools cis --vcf ../../Genotype2/SplX_converted_noHH_mQTL_snpfilt_hg38_sort.vcf.gz --bed ../../Quants/all_meth_nona_mQTL_filt_nochr.bed.gz --cov ../Cov_2PC_geno_1-4_7-10_12_14_PC_meth.txt.gz --nominal 0.01  --chunk $j 200 --out chunk_$j\_200.txt --normal --window 1000000
EOF
done

cat chunk*_200.txt > nominals_full.txt
wc -l nominals_full.txt
 # 3,753,006 mQTLs (pval <.01) found in nominals_full.txt

 # top eQTL for gene
awk '$14==1'  nominals_full.txt | wc -l
 # 624,929 top mQTLs for CpGs
```
PC 14 doesn't add anything.


Let's try adding back in PC 6 / 5

```bash
cd ../
mkdir geno_2_14_610_12_14_meth
cd geno_2_14_610_12_14_meth

 # running the permutation
for j in $(seq 1 200);
do
qsub -S /bin/bash -N mQTL_${j}_chunk_nom -cwd -l h_vmem=5.6G -j y << EOF
module load QTLtools/1.1/gcc.4.7.4
QTLtools cis --vcf ../../Genotype2/SplX_converted_noHH_mQTL_snpfilt_hg38_sort.vcf.gz --bed ../../Quants/all_meth_nona_mQTL_filt_nochr.bed.gz --cov ../Cov_2PC_geno_1-4_6-10_12_PC_meth.txt.gz --nominal 0.01  --chunk $j 200 --out chunk_$j\_200.txt --normal --window 1000000
EOF
done

cat chunk*_200.txt > nominals_full.txt
wc -l nominals_full.txt
 # 3,570,091 mQTLs (pval <.01) found in nominals_full.txt

 # top eQTL for gene
awk '$14==1'  nominals_full.txt | wc -l
 # 591,101 top mQTLs for CpGs
```

```bash
cd ../
mkdir geno_2_15_710_12_14_meth
cd geno_2_15_710_12_14_meth

 # running the permutation
for j in $(seq 1 200);
do
qsub -S /bin/bash -N mQTL_${j}_chunk_nom -cwd -l h_vmem=5.6G -j y << EOF
module load QTLtools/1.1/gcc.4.7.4
QTLtools cis --vcf ../../Genotype2/SplX_converted_noHH_mQTL_snpfilt_hg38_sort.vcf.gz --bed ../../Quants/all_meth_nona_mQTL_filt_nochr.bed.gz --cov ../Cov_2PC_geno_1-5_7-10_12_PC_meth.txt.gz --nominal 0.01  --chunk $j 200 --out chunk_$j\_200.txt --normal --window 1000000
EOF
done

cat chunk*_200.txt > nominals_full.txt
wc -l nominals_full.txt
 # 3,753,006 mQTLs (pval <.01) found in nominals_full.txt

 # top eQTL for gene
awk '$14==1'  nominals_full.txt | wc -l
 # 624,929 top mQTLs for CpGs
cd ../../
```

<a name="permmqtl"/>

7. Permutation mQTL calling

It seems that using the first 10 expression PCs and the first 2 genotype PCs gives the most mQTLs in the result. So I ran QTL analysis through 10,000 permutations to generate a "final" eQTL list.

```bash
mkdir Permutation
cd /home/greally-lab/Deepa_Andrew/eQTL/Permutation2
mkdir perm_pc10_
cd perm_g12_123471012

 # call the permutations
for j in $(seq 1 200);
do
qsub -S /bin/bash -N mQTL_${j}_chunk -cwd -R y -l h_vmem=5.6G -j y << EOF
module load QTLtools/1.1/gcc.4.7.4
QTLtools cis --vcf ../../Genotype2/SplX_converted_noHH_mQTL_snpfilt_hg38_sort.vcf.gz --bed ../../Quants/all_meth_nona_mQTL_filt_2_nochr.bed.gz --cov ../../Nominal/Cov_2PC_geno_1-4_6-10_12_PC_meth.txt.gz --permute 10000  --chunk $j 200 --out perms_$j\_200.txt --normal --window 1000000 --seed 123456
EOF
done

 # combine the permutation files together
cat perms_*200.txt > mQTL_final1.txt
wc -l mQTL_final1.txt
 # 732,350 mQTL_final1.txt
 # 732,345

sed "146976q;d" mQTL_final1.txt
grep "CpG_127393" perms*
perms_128_200.txt # took out CpG_127393
154494
sed "154494q;d" mQTL_final1.txt
CpG_701621
grep "CpG_701621" perms*
vi perms_129_200.txt
took out CpG_701640

module load R/3.4.0/gcc.4.7.4
cp ../../../eQTL/Permutation2/perm_pc10_manyperm/plot_beta.R .
Rscript plot_beta.R mQTL_final1.txt G_12_E_1-4_7-10_12_manyperm2
 # Found 2166 significant QTls.

gzip mQTL_final1.txt
cp ../../../eQTL/Permutation2/perm_pc10_manyperm/RunFDR_cis.R .
Rscript RunFDR_cis.R mQTL_final1.txt.gz 0.05 G_12_E_1-4_7-10_12_final

wc -l G_12_E_1-4_7-10_12_final.significant.txt
 #443 G_12_E_10_final.significant.txt
```

Now we have a list of 443 variants (top variants of CpGs), which pass the significance threshold. Now we'll make a conditional pass to find other significant varaints, which are independent of the top variant

```bash
cd ../
mkdir perm_g12_123471012_cond
cd perm_g12_123471012_cond

for j in $(seq 1 200);
do
qsub -S /bin/bash -N mQTL_${j}_chunk -cwd -R y -l h_vmem=5.6G -j y << EOF
module load QTLtools/1.1/gcc.4.7.4
QTLtools cis --vcf ../../Genotype2/SplX_converted_noHH_mQTL_snpfilt_hg38_sort.vcf.gz --bed ../../Quants/all_meth_nona_mQTL_filt_2_nochr.bed.gz --cov ../../Nominal/Cov_2PC_geno_1-4_6-10_12_PC_meth.txt.gz --mapping ../perm_g12_123471012/G_12_E_1-4_7-10_12_final.thresholds.txt --normal --window 1000000 --chunk $j 200 --out cond_$j\_200.txt
EOF
done

cat cond_*200.txt > mQTL_final2.txt
wc -l mQTL_final2.txt
 # 4124 mQTL_final2.txt
 # take only those that pass backward significance
awk '$20==1' mQTL_final2.txt > mQTL_final3.txt

 #how many SNPs
wc -l mQTL_final3.txt
 # 4124 mQTL_final3.txt

 # how many CpGs?
awk '$19==1' mQTL_final3.txt | wc -l #2168
awk '{print $1}' mQTL_final3.txt | uniq | wc -l #2166, which means there's a tie 2 places?
```

examine these 4124 mQTLs and 2166 mGenes.

<a name="mqtlgene"/>

8. Connecting CpGs to Genes

First I need to annotate the CpGs to a gene and then select out CpGs "linked" to the genes of interest. Also, how far from the CpG is the variant?

```bash
cd /home/greally-lab/Deepa_Andrew/mQTL/Validate/Anno_CpG
cp ../../Permutation/perm_g12_123471012_cond/mQTL_final3.txt .
cp ../../../SMITE_2018/CpG_Anno_wEnh50kb_info.txt .
```

Then see how I annotated the genes in `Exploring_mQTL.rmd`



First see if any of the pathway genes came up.
```
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
PAK1	ENSG00000149269
```

```bash
 # which genes were tested? (all above not commented on)
for f1 in ENSG00000160293 ENSG00000147459 ENSG00000077264 ENSG00000075651 ENSG00000179604 ENSG00000070831 ENSG00000173327 ENSG00000133895 ENSG00000121879 ENSG00000171608 ENSG00000051382 ENSG00000105647 ENSG00000145675 ENSG00000149269;
do
grep $f1 mQTL_final3_ENS.txt
done
 # PLD1 is hit; ENSG00000075651; all the same CpG but with three different variants
 # 464549	3	171676808	171676808	+	453	-77633	JHU_3.171316964	3	171599175	171599175	0	0.000000000121068	0.992401	0	1	0.000000000121068	0.992401	0	1	NA,ENSG00000075651,NA
 # 464549	3	171676808	171676808	+	453	-24822	rs4894707	3	171651986	171651986	0	0.0000000000846013	0.99158	1	1	0.0000000000846013	0.99158	1	1	NA,ENSG00000075651,NA
 # 464549	3	171676808	171676808	+	453	870	rs2124147	3	171677678	171677678	0	0.0000000000976631	0.946823	0	1	0.0000000000976631	0.946823	0	1	NA,ENSG00000075651,NA

```

grab all of the ENSIDs from the mGenes and pass through string.

```bash
 #grabbing IDs
tr , '\n' < mQTL_linked_genes_1.txt > mQTL_linked_genes_2.txt
wc -l mQTL_linked_genes_2.txt
 # 3194 mQTL_linked_genes_2.txt

 # making a unique only gene list
cat mQTL_linked_genes_2.txt | sort | uniq | wc -l # 1268
cat mQTL_linked_genes_2.txt | sort | uniq > mQTL_linked_genes_3.txt

```

how many of these genes are also within the eGenes (genes associated with an eQTL)?
```bash
grep -f mQTL_linked_genes_2.txt  ../../../eQTL/Validate2/hc_genes/eQTL_final3.txt | awk '{print $1}' | sort | uniq | wc -l # 49

grep -f mQTL_linked_genes_2.txt  ../../../eQTL/Validate2/hc_genes/eQTL_final3.txt | awk '{print $1}' | sort | uniq > e-mQTL_Genes.txt
```

We also want to check to see if genotype is contributing to the differential methylation of any of the genes found within the original 89 genes identifeid via RNAseq. These genes can be found in the following file: `/home/greally-lab/Deepa_Andrew/04-2018_Helptagging_cellprop_analysis/04_2018_hcgenes_rnaseqpaper_fix.txt`

```bash
grep -F -f /home/greally-lab/Deepa_Andrew/04-2018_Helptagging_cellprop_analysis/04_2018_hcgenes_rnaseqpaper_fix.txt mQTL_final3_ENS.txt > mQTL_final3_ENS_hcgenes.txt
 # ENSG00000198752 (CDC42BPB), ENSG00000175866 (BAIAP2), ENSG00000169184 (MN1), ENSG00000163803 (PLB1), ENSG00000135636 (DYSF), ENSG00000117115 (PADI2), ENSG00000075651 (PLD1), ENSG00000109265 (KIAA1211), ENSG00000173406 (DAB1), ENSG00000106069 (CHN2), ENSG00000221866 (PLXNA4), ENSG00000147488 (ST18)
```


<a name="plotmqtl"/>

9. Plotting mQTLs

make some graphs for the various hc genes found as mQTL-genes

```bash
 # creating the dir and relocating files for ease
mkdir ../SNP_plots
cd ../SNP_plots
mkdir plots
cp ../../Quants/all_meth_nona_mQTL_filt_2_nochr.bed.gz  .
gunzip -d all_meth_nona_mQTL_filt_2_nochr.bed.gz
cp ../../Bams/PCA/sample_metadata_mQTL_3.csv .
cp ../../Genotype2/SplX_converted_noHH_mQTL_snpfilt_hg38_sort.vcf.gz .
gunzip -d SplX_converted_noHH_mQTL_snpfilt_hg38_sort.vcf.gz
cp ../Anno_CpG/mQTL_final3_ENS.txt .

module load R/3.4.0/gcc.4.7.4
```
We wanted to look at the genes from the original paper that were found as mQTL-associated genes: `mQTL_associated_genes.txt`

```
hgnc_symbol ensembl_gene_id
CDC42BPB	ENSG00000198752
BAIAP2	ENSG00000175866
MN1	ENSG00000169184
PLB1	ENSG00000163803
DYSF	ENSG00000135636
PADI2	ENSG00000117115
PLD1	ENSG00000075651
KIAA1211	ENSG00000109265
DAB1	ENSG00000173406
CHN2	ENSG00000106069
PLXNA4	ENSG00000221866
ST18	ENSG00000147488
```

creating the plots for the original paper genes associated with mQTL-CpGs

```bash
while read -r b c; do
	qsub -S /bin/bash -N plot_${b}_eQTL -cwd -R y -l h_vmem=5.6G -j y << EOF
	module load R/3.4.0/gcc.4.7.4
  Rscript Plot_mQTL.R $b $c mQTL_final3_ENS.txt SplX_converted_noHH_mQTL_snpfilt_hg38_sort.vcf sample_metadata_mQTL_3.csv all_meth_nona_mQTL_filt_2_nochr.bed
EOF
done < mQTL_associated_genes.txt
```


We wanted to look at some other genes: `Genes_of_interest.txt`

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
	Rscript Plot_mQTL.R $b $c mQTL_final3_ENS.txt SplX_converted_noHH_mQTL_snpfilt_hg38_sort.vcf sample_metadata_mQTL_3.csv all_meth_nona_mQTL_filt_2_nochr.bed
EOF
done < Genes_of_interest.txt
```

I moved all of the plots into subdirectory `plots/`
```bash
mv ENSG00000* plots/
mkdir STDOUT
mv plot_*.o* STDOUT/
```


Wanted to look at me-genes  w/ diff methylayed CpGs: `DM-me-genes.txt`

4 me-genes contained differentially methylated CpGs on the gene promoter (ZNF7, ACER3, DYM, and C19orf53), 1 me-gene overlapped with differentially methylated CpGs on the gene body (PLA2G6) and 2 me-genes had differentially methylated CpGs on a cis-regulatory region (TPO, GDF10)

```
hgnc_symbol ensembl_gene_id
ZNF7  ENSG00000147789
ACER3 ENSG00000078124
DYM ENSG00000141627
C19orf53  ENSG00000104979
PLA2G6  ENSG00000184381
TPO ENSG00000115705
GDF10 ENSG00000266524
```

while read -r b c; do
	qsub -S /bin/bash -N plot_${b}_mQTL -cwd -R y -l h_vmem=5.6G -j y << EOF
	module load R/3.4.0/gcc.4.7.4
	Rscript Plot_mQTL.R $b $c mQTL_final3_ENS.txt SplX_converted_noHH_mQTL_snpfilt_hg38_sort.vcf sample_metadata_mQTL_3.csv all_meth_nona_mQTL_filt_2_nochr.bed
EOF
done < DM-me-genes.txt

scratch
```
gene_name <- "ACER3"
geneID <- "ENSG00000078124"
mqtlFile <- "mQTL_final3_ENS.txt"
vcfFile <- "SplX_converted_noHH_mQTL_snpfilt_hg38_sort.vcf"
covFile <- "sample_metadata_mQTL_3.csv"
methFile <- "all_meth_nona_mQTL_filt_2_nochr.bed"
```

me-genes that have differentially expressed genes
```
hgnc_symbol ensembl_gene_id
PPP2R2C ENSG00000074211
ADAMTS2 ENSG00000087116
PACSIN2 ENSG00000100266
PRKCG ENSG00000126583
MAP1B ENSG00000131711
ALDH1L1  ENSG00000144908
ITM2C ENSG00000135916
CA10  ENSG00000154975
KLHDC2 ENSG00000165516
LARGE2  ENSG00000165905
ANKRD33 ENSG00000167612
MYO1C ENSG00000197879
```

<a name="otherstudy"/>

10. Overlap with Blueprint

Validating our mQTLs against previous study (Blueprint data - 450k array/whole genome sequencing). Need to keep in mind that this dataset is Caucasian as well.

The columns of `tcel_meth_M_peer_10_all_summary.txt` are as follows:
```
column_number	column_label
1		chr:pos_ref_alt
2		rsid
3		phenotypeID
4		p.value
5		beta
6		Bonferroni.p.value
7		FDR
8		alt_allele_frequency
9		std.error_of_beta
```

```bash
cd /home/greally-lab/Deepa_Andrew/mQTL/Validate/Other_studies
wc -l tcel_meth_M_peer_10_all_summary.txt
 # 1,679,265,698
 # grab only FDR <0.05
awk '$7 <0.05' tcel_meth_M_peer_10_all_summary.txt > tcel_meth_M_peer_10_all_summary_fdr05.txt
awk '$7 <0.05' tcel_meth_M_peer_10_all_summary.txt |less

awk '($7+0) <0.05' tcel_meth_M_peer_10_all_summary.txt | wc -l
 # 12,656,798
```
That's a ton of mQTLs.

Seeing how many variants are in both



#### Making manhattan plots for mQTLs
in `/home/greally-lab/Deepa_Andrew/mQTL/Validate/SNP_plots`
make Knit_manhat_mQTL.R file
```
library(knitr)
library(rmarkdown)
rmarkdown::render("mQTL_downstream_analysis2.Rmd")
```

```bash
qsub -S /bin/bash -N knit_manhat -j y -cwd -q highmem.q -l h_vmem=300G << EOF
module load R/3.4.0/gcc.4.7.4
module load pandoc/1.19.2.1
module load texlive/2016
cd /home/greally-lab/Deepa_Andrew/mQTL/Validate/SNP_plots
R CMD BATCH --no-restore Knit_manhat_mQTL.R
EOF
```
