# Annotating CpGs
## 08/28/18
### Andrew D. Johnston

This file serves as a centralized copy of the code used to annotate the CpGs used in the study. The CpGs were marked as lying within either promoter, enhancer, or genebody - when overlaps of annotation arose, the respective order was applied.

In summary, ensembl protein coding gene annotations were used to mark promoters and genebodies while enhancers were defined by a combination of ATACseq results from CD4 T cells (Buenrostro et al. 2013) and the ensembl regulatory build for CD4 TCell Venous blood (Zerbino et al. 2015).


1. Annotating promoters and gene bodies

I downloaded the ensembl build 83 from the web portal.

```bash
I need to get the refseq annotations though via ensembl

/home/greally-lab/indexes/hg38/ensembl/Homo_sapiens.GRCh38.83.gtf

awk '$3~"gene" {print $0}' Homo_sapiens.GRCh38.83.gtf | awk '$18~"protein_coding" {print $0}' | wc -l # 19826

## went back and changed the promoters to +/- 1kb
awk '$3~"gene" {print $0}' Homo_sapiens.GRCh38.83.gtf | awk '$18~"protein_coding" {print $0}' | awk '$7=="+" {print $0}' | awk '{if ($5-$4 >= 1000) print "chr"$1"\t"$4-1000"\t"$4+1000"\t"$10"\t"".""\t"$7; else print "chr"$1"\t"$4-1000"\t"$5"\t"$10"\t"".""\t"$7}' | tr -d '";' > ensembl_prom_plus.bed

# are any starts < 0 ?
awk '$2 < 0 {print $0}' ensembl_prom_plus.bed | less #none


awk '$3~"gene" {print $0}' Homo_sapiens.GRCh38.83.gtf | awk '$18~"protein_coding" {print $0}' | awk '$7=="-" {print $0}' | awk '{if ($5-$4 >= 1000) print "chr"$1"\t"$5-1000"\t"$5+1000"\t"$10"\t"".""\t"$7; else print "chr"$1"\t"$4"\t"$5+1000"\t"$10"\t"".""\t"$7}' | tr -d '";' > ensembl_prom_neg.bed

awk '$2 < 0 {print $0}' ensembl_prom_neg.bed | less #none



#Get Gene bodies
# If genes are shorter than 1000 bp, only the promoter region of that gene was annotated.
awk '$3~"gene" {print $0}' Homo_sapiens.GRCh38.83.gtf | awk '$18~"protein_coding" {print $0}' | awk '$7=="+" {print $0}' | awk '{if ($5-$4 > 1000) print "chr"$1"\t"$4+1000"\t"$5"\t"$10"\t"".""\t"$7}' | tr -d '";' > ensembl_gene_plus.bed

#check
module load bedtools2
bedtools intersect -a ensembl_prom_plus.bed -b ensembl_gene_plus.bed -wo | wc -l
#1720 but these are other close-by genes overlapping when promoters are +/- 2kb
# 1507 when promoters are +/- 1kb
# I'll have to annotate the CpGs and prioritize them based on enhancer>prom>genebody

#example
awk '$3~"gene" {print $0}' Homo_sapiens.GRCh38.83.gtf | awk '$18~"protein_coding" {print $0}' | awk '$7=="+" {print $0}' | less

1       ensembl_havana  gene    966497  975865  .       +       .       gene_id "ENSG00000187583"

chr1    964497  968497  ENSG00000187583 .       +       chr1    962587  965715  ENSG00000187961 .       +       1218

1       ensembl_havana  gene    960587  965715  .       +       .       gene_id "ENSG00000187961"

#get neg gene bodies
awk '$3~"gene" {print $0}' Homo_sapiens.GRCh38.83.gtf | awk '$18~"protein_coding" {print $0}' | awk '$7=="-" {print $0}' | awk '{if ($5-$4 > 1000) print "chr"$1"\t"$4"\t"$5-1000"\t"$10"\t"".""\t"$7}' | tr -d '";' > ensembl_gene_neg.bed

cat ensembl_gene_plus.bed ensembl_gene_neg.bed > ensembl_gene.bed
cat ensembl_prom_plus.bed ensembl_prom_neg.bed > ensembl_prom.bed

sort -k1,1 -k2,2n ensembl_gene.bed > ensembl_gene_sorted.bed
wc -l ensembl_gene_sorted.bed #19,254
bedtools merge -i ensembl_gene_sorted.bed > ensembl_gene_merge.bed
wc -l ensembl_gene_merge.bed # 16042 ensembl_gene_merge.bed
# 3,212 gene annotations are lost in the merge and are therefore overlapping. (<6,424 are overlapping)


sort -k1,1 -k2,2n ensembl_prom.bed > ensembl_prom_sorted.bed
wc -l ensembl_prom_sorted.bed # 19826
bedtools merge -i ensembl_prom_sorted.bed > ensembl_prom_merge.bed
wc -l ensembl_prom_merge.bed # 17907

/home/greally-lab/indexes/hg38/ensembl

cd  /home/greally-lab/Deepa_Andrew/Deepa-helptagging
awk 'NR>1 {print $2"\t"$3"\t"$3+1"\t"$1}' help_annot > help_annot.bed

wc -l help_annot.bed
# 2,455,516 help_annot.bed

# bedtools intersect -u -wa -a help_annot.bed -b /home/greally-lab/indexes/hg38/ensembl/ensembl_prom_merge.bed | wc -l
bedtools intersect -wo -a help_annot.bed -b /home/greally-lab/indexes/hg38/ensembl/ensembl_prom_sorted.bed > help_annot_prom.bed
wc -l help_annot_prom.bed
#	previously.. but 268,104 help_annot_prom.bed, but I have no idea how this number changed.
#  219,281 help_annot_prom.bed

# bedtools intersect -u -wa -a help_annot.bed -b /home/greally-lab/indexes/hg38/ensembl/ensembl_gene_sorted.bed | wc -l
bedtools intersect -wo -a help_annot.bed -b /home/greally-lab/indexes/hg38/ensembl/ensembl_gene_sorted.bed > help_annot_gene.bed
wc -l help_annot_gene.bed
#	1,162,017 help_annot_gene.bed
```





(Roadmap Epigenomics Consortium et al. 2015)


references:
Buenrostro, J. D., Giresi, P. G., Zaba, L. C., Chang, H. Y., & Greenleaf, W. J. (2013). Transposition of native chromatin for fast and sensitive epigenomic profiling of open chromatin, DNA-binding proteins and nucleosome position. Nature Methods, 10(12), 1213–1218. http://doi.org/10.1038/nmeth.2688
Roadmap Epigenomics Consortium et al. (2015). Integrative analysis of 111 reference human epigenomes., 518(7539), 317–330. http://doi.org/10.1038/nature14248
The ensembl regulatory build. (2015). The ensembl regulatory build., 16, 56. http://doi.org/10.1186/s13059-015-0621-5

Integrative analysis of 111 reference human epigenomes. (2015). Integrative analysis of 111 reference human epigenomes., 518(7539), 317–330. http://doi.org/10.1038/nature14248

downloaded the the regulatory build from ftp://ftp.ensembl.org/pub/release-90/regulation/homo_sapiens/RegulatoryFeatureActivity/. I got the CD4 TCell Venous blood build.
