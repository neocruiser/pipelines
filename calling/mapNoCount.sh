#! /user/bin/bash

:'
this script accomplish 5 things:
1. map all paired end samples to reference with bwa
2. sort the mapped contigs with samtools
3. remove duplicate contigs with picard
4. index contigs with samtools
5. count contigs with htseq
'

java -jar /home/neo/data/Trimmomatic-0.33/trimmomatic-0.33.jar PE \
/media/sf_docs/data/QPX-RNA-Seq/mmetsp0098.1.NY.fastq.gz \
/media/sf_docs/data/QPX-RNA-Seq/mmetsp0098.2.NY.fastq.gz \
/media/sf_docs/data/QPX-RNA-Seq/trimmed/mmetsp0098.1.trimmed.P.NY.fastq.gz \
/media/sf_docs/data/QPX-RNA-Seq/trimmed/mmetsp0098.1.trimmed.U.NY.fastq.gz \
/media/sf_docs/data/QPX-RNA-Seq/trimmed/mmetsp0098.2.trimmed.P.NY.fastq.gz \
/media/sf_docs/data/QPX-RNA-Seq/trimmed/mmetsp0098.2.trimmed.U.NY.fastq.gz \
ILLUMINACLIP:TrueSeq3-PE-3.fa:2:30:10 \
SLIDINGWINDOW:4:15 \
TRAILING:5 \
MINLEN:45

java -jar /home/neo/data/Trimmomatic-0.33/trimmomatic-0.33.jar PE \
/media/sf_docs/data/QPX-RNA-Seq/mmetsp00992.1.NY.fastq.gz \
/media/sf_docs/data/QPX-RNA-Seq/mmetsp00992.2.NY.fastq.gz \
/media/sf_docs/data/QPX-RNA-Seq/trimmed/mmetsp00992.1.trimmed.P.NY.fastq.gz \
/media/sf_docs/data/QPX-RNA-Seq/trimmed/mmetsp00992.1.trimmed.U.NY.fastq.gz \
/media/sf_docs/data/QPX-RNA-Seq/trimmed/mmetsp00992.2.trimmed.P.NY.fastq.gz \
/media/sf_docs/data/QPX-RNA-Seq/trimmed/mmetsp00992.2.trimmed.U.NY.fastq.gz \
ILLUMINACLIP:TrueSeq3-PE-3.fa:2:30:10 \
SLIDINGWINDOW:4:15 \
TRAILING:5 \
MINLEN:45

java -jar /home/neo/data/Trimmomatic-0.33/trimmomatic-0.33.jar PE \
/media/sf_docs/data/QPX-RNA-Seq/mmetsp001002.1.NY.fastq.gz \
/media/sf_docs/data/QPX-RNA-Seq/mmetsp001002.2.NY.fastq.gz \
/media/sf_docs/data/QPX-RNA-Seq/trimmed/mmetsp001002.1.trimmed.P.NY.fastq.gz \
/media/sf_docs/data/QPX-RNA-Seq/trimmed/mmetsp001002.1.trimmed.U.NY.fastq.gz \
/media/sf_docs/data/QPX-RNA-Seq/trimmed/mmetsp001002.2.trimmed.P.NY.fastq.gz \
/media/sf_docs/data/QPX-RNA-Seq/trimmed/mmetsp001002.2.trimmed.U.NY.fastq.gz \
ILLUMINACLIP:TrueSeq3-PE-3.fa:2:30:10 \
SLIDINGWINDOW:4:15 \
TRAILING:5 \
MINLEN:45

java -jar /home/neo/data/Trimmomatic-0.33/trimmomatic-0.33.jar PE \
/media/sf_docs/data/QPX-RNA-Seq/mmetsp001433.1.NY.fastq.gz \
/media/sf_docs/data/QPX-RNA-Seq/mmetsp001433.2.NY.fastq.gz \
/media/sf_docs/data/QPX-RNA-Seq/trimmed/mmetsp001433.1.trimmed.P.NY.fastq.gz \
/media/sf_docs/data/QPX-RNA-Seq/trimmed/mmetsp001433.1.trimmed.U.NY.fastq.gz \
/media/sf_docs/data/QPX-RNA-Seq/trimmed/mmetsp001433.2.trimmed.P.NY.fastq.gz \
/media/sf_docs/data/QPX-RNA-Seq/trimmed/mmetsp001433.2.trimmed.U.NY.fastq.gz \
ILLUMINACLIP:TrueSeq3-PE-3.fa:2:30:10 \
SLIDINGWINDOW:4:15 \
TRAILING:5 \
MINLEN:45

java -jar /home/neo/data/Trimmomatic-0.33/trimmomatic-0.33.jar PE \
/media/sf_docs/data/QPX-RNA-Seq/mmetsp0099.1.NY.fastq.gz \
/media/sf_docs/data/QPX-RNA-Seq/mmetsp0099.2.NY.fastq.gz \
/media/sf_docs/data/QPX-RNA-Seq/trimmed/mmetsp0099.1.trimmed.P.NY.fastq.gz \
/media/sf_docs/data/QPX-RNA-Seq/trimmed/mmetsp0099.1.trimmed.U.NY.fastq.gz \
/media/sf_docs/data/QPX-RNA-Seq/trimmed/mmetsp0099.2.trimmed.P.NY.fastq.gz \
/media/sf_docs/data/QPX-RNA-Seq/trimmed/mmetsp0099.2.trimmed.U.NY.fastq.gz \
ILLUMINACLIP:TrueSeq3-PE-3.fa:2:30:10 \
SLIDINGWINDOW:4:15 \
TRAILING:5 \
MINLEN:45

java -jar /home/neo/data/Trimmomatic-0.33/trimmomatic-0.33.jar PE \
/media/sf_docs/data/QPX-RNA-Seq/mmetsp00100.1.NY.fastq.gz \
/media/sf_docs/data/QPX-RNA-Seq/mmetsp00100.2.NY.fastq.gz \
/media/sf_docs/data/QPX-RNA-Seq/trimmed/mmetsp00100.1.trimmed.P.NY.fastq.gz \
/media/sf_docs/data/QPX-RNA-Seq/trimmed/mmetsp00100.1.trimmed.U.NY.fastq.gz \
/media/sf_docs/data/QPX-RNA-Seq/trimmed/mmetsp00100.2.trimmed.P.NY.fastq.gz \
/media/sf_docs/data/QPX-RNA-Seq/trimmed/mmetsp00100.2.trimmed.U.NY.fastq.gz \
ILLUMINACLIP:TrueSeq3-PE-3.fa:2:30:10 \
SLIDINGWINDOW:4:15 \
TRAILING:5 \
MINLEN:45


sample[1]=mmetsp0098
sample[2]=mmetsp001433
sample[3]=mmetsp00992
sample[4]=mmetsp001002
sample[5]=mmetsp0099
sample[6]=mmetsp00100

ir=/media/sf_docs/data/QPX-RNA-Seq/trimmed
dir=/media/sf_docs/data/mappingX
ddir=/media/sf_docs/data/rmdupX

extension=.trimmed.P.NY.fastq.gz
reference=/media/sf_docs/data/QPX-RNA-Seq/Steve_Roberts/QPXTranscriptome_v21/QPX_transcriptome_v2orf.fasta
count=/media/sf_docs/data/QPX-RNA-Seq/Steve_Roberts/QPXTranscriptome_v21/QPX_transcriptome_v2orf.gff3

for i in 1 2 3 4 5 6
do
    sample=${sample[${i}]}
    bwa mem ${reference} \
        ${ir}/${sample}.1${extension} \
        ${ir}/${sample}.2${extension} | \
        samtools view -Shu - | \
        samtools sort - ${dir}/${sample}.sorted

    samtools index ${ddir}/${sample}.sorted.bam

    java -Xmx2g -jar /home/neo/data/picard/picard.jar \
        MarkDuplicates \
        INPUT=${dir}/${sample}.sorted.bam \
        OUTPUT=${ddir}/${sample}.nodup.bam \
        METRICS_FILE=${ddir}/${sample}.dup.metrics \
        REMOVE_DUPLICATES=true \
        ASSUME_SORTED=true

done


sample[1]=mmetsp0098
sample[2]=mmetsp001433
sample[3]=mmetsp00992
sample[4]=mmetsp001002
sample[5]=mmetsp0099
sample[6]=mmetsp00100

ir=/media/sf_docs/data/QPX-RNA-Seq/trimmed
dir=/media/sf_docs/data/mappingX2
ddir=/media/sf_docs/data/rmdupX2
extension=.trimmed.P.NY.fastq.gz

reference=/media/sf_docs/data/mmetsp0098/contigs.fa
count=/media/sf_docs/data/mmetsp0098/MMETSP0098.gff3

for i in 1 2 3 4 5 6
do
    sample=${sample[${i}]}
    bwa mem ${reference} \
        ${ir}/${sample}.1${extension} \
        ${ir}/${sample}.2${extension} | \
        samtools view -Shu - | \
        samtools sort - ${dir}/${sample}.sorted

    samtools index ${ddir}/${sample}.sorted.bam

    java -Xmx2g -jar /home/neo/data/picard/picard.jar \
        MarkDuplicates \
        INPUT=${dir}/${sample}.sorted.bam \
        OUTPUT=${ddir}/${sample}.nodup.bam \
        METRICS_FILE=${ddir}/${sample}.dup.metrics \
        REMOVE_DUPLICATES=true \
        ASSUME_SORTED=true

done



sample[1]=mmetsp0098
sample[2]=mmetsp001433
sample[3]=mmetsp00992
sample[4]=mmetsp001002
sample[5]=mmetsp0099
sample[6]=mmetsp00100

ir=/media/sf_docs/data/QPX-RNA-Seq/trimmed
dir=/media/sf_docs/data/mappingX3
ddir=/media/sf_docs/data/rmdupX3
extension=.trimmed.P.NY.fastq.gz

reference=/media/sf_docs/data/genomeSRv015/QPX_v015.fasta
count=/media/sf_docs/data/genomeSRv015/QPX_v015.gff3

for i in 1 2 3 4 5 6
do
    sample=${sample[${i}]}
    bwa mem ${reference} \
        ${ir}/${sample}.1${extension} \
        ${ir}/${sample}.2${extension} | \
        samtools view -Shu - | \
        samtools sort - ${dir}/${sample}.sorted

    samtools index ${ddir}/${sample}.sorted.bam

    java -Xmx2g -jar /home/neo/data/picard/picard.jar \
        MarkDuplicates \
        INPUT=${dir}/${sample}.sorted.bam \
        OUTPUT=${ddir}/${sample}.nodup.bam \
        METRICS_FILE=${ddir}/${sample}.dup.metrics \
        REMOVE_DUPLICATES=true \
        ASSUME_SORTED=true

done
