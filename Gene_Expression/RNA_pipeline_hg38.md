Command used to generate STAR aligned counts from trimmed fastq files using hg38 reference genome and gtf.

```bash
mkdir Mapped_Ensembl
for f1 in fastq_trimmed/*_trimmed.fq.gz   
do
FILE="${f1##*/}"
SAMPLE="$(echo ${FILE} | cut -d '_' -f5 | cut -d '.' -f1)"
echo ${SAMPLE}
echo ${FILE}
qsub -S /bin/bash -N ${SAMPLE}_alignEns -cwd -l h_vmem=5.6G -j y -pe smp 20 << EOF
module load samtools/1.2/gcc.4.4.7
module load STAR/2.5.2b/gcc.5.3.0
STAR --runThreadN 20 --genomeDir /home/greally-lab/indexes/hg38/Star --readFilesIn ${f1} --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --chimSegmentMin 20 --quantMode TranscriptomeSAM GeneCounts --quantTranscriptomeBan Singleend --sjdbGTFfile /home/greally-lab/indexes/hg38/ensembl/Homo_sapiens.GRCh38.83.gtf --sjdbGTFchrPrefix chr --sjdbOverhang 99 --outFileNamePrefix Mapped_Ensembl/${SAMPLE}
EOF
done
```
