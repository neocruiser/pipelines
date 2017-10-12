#! /user/bin/bash

:'
this script accomplish 5 things:
1. map all paired end samples to reference with bwa
2. sort the mapped contigs with samtools
3. remove duplicate contigs with picard
4. index contigs with samtools
5. count contigs with htseq
-M: bwa mark shorter hits as secondary, increase picard comaptibility
'

sample[1]=mmetsp0098
sample[2]=mmetsp001433
sample[3]=mmetsp00992
sample[4]=mmetsp001002
sample[5]=mmetsp0099
sample[6]=mmetsp00100

ir=/media/sf_docs/data/QPX-RNA-Seq/trimmed
dir=/media/sf_docs/data/mappingY
ddir=/media/sf_docs/data/rmdupY

extension=.trimmed.P.NY.fastq.gz
reference=/media/sf_docs/data/QPX-RNA-Seq/Steve_Roberts/QPXTranscriptome_v21/QPX_transcriptome_v2orf.fasta

RG[1]='@RG\tID:mmetsp0098\tSM:NY1\tPL:illumina\tLB:mmetsp0098\tPU:QPXtrxSRv21'
RG[2]='@RG\tID:mmetsp001433\tSM:NY1\tPL:illumina\tLB:mmetsp001433\tPU:QPXtrxSRv21'
RG[3]='@RG\tID:mmetsp00992\tSM:MA1\tPL:illumina\tLB:mmetsp00992\tPU:QPXtrxSRv21'
RG[4]='@RG\tID:mmetsp001002\tSM:VA1\tPL:illumina\tLB:mmetsp001002\tPU:QPXtrxSRv21'
RG[5]='@RG\tID:mmetsp0099\tSM:MA2\tPL:illumina\tLB:mmetsp0099\tPU:QPXtrxSRv21'
RG[6]='@RG\tID:mmetsp00100\tSM:VA2\tPL:illumina\tLB:mmetsp00100\tPU:QPXtrxSRv21'

    java -jar /home/neo/data/picard/picard.jar \
        CreateSequenceDictionary \
        R=${reference} \
        O=/media/sf_docs/data/QPX-RNA-Seq/Steve_Roberts/QPXTranscriptome_v21/QPX_transcriptome_v2orf.dict

    samtools faidx ${reference}

for i in 1 2 3 4 5 6
do
    sample=${sample[${i}]}
    RG=${RG[${i}]}
    bwa mem -M \
        -R ${RG} \
        -p ${reference} \
        ${ir}/${sample}.1${extension} \
        ${ir}/${sample}.2${extension} \
    > ${dir}/${sample}.sam

    java -Xmx2g -jar /home/neo/data/picard/picard.jar \
        SortSam \
        INPUT=${dir}/${sample}.sam \
        OUTPUT=${ddir}/${sample}.sorted.bam \
        SORT_ORDER=coordinate

    java -Xmx2g -jar /home/neo/data/picard/picard.jar \
        MarkDuplicates \
        INPUT=${ddir}/${sample}.sorted.bam \
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
dir=/media/sf_docs/data/mappingY2
ddir=/media/sf_docs/data/rmdupY2

extension=.trimmed.P.NY.fastq.gz
reference=/media/sf_docs/data/mmetsp0098/contigs.fa

RG[1]='@RG\tID:mmetsp0098\tSM:NY1\tPL:illumina\tLB:mmetsp0098\tPU:MMETSP0098'
RG[2]='@RG\tID:mmetsp001433\tSM:NY1\tPL:illumina\tLB:mmetsp001433\tPU:MMETSP0098'
RG[3]='@RG\tID:mmetsp00992\tSM:MA1\tPL:illumina\tLB:mmetsp00992\tPU:MMETSP0098'
RG[4]='@RG\tID:mmetsp001002\tSM:VA1\tPL:illumina\tLB:mmetsp001002\tPU:MMETSP0098'
RG[5]='@RG\tID:mmetsp0099\tSM:MA2\tPL:illumina\tLB:mmetsp0099\tPU:MMETSP0098'
RG[6]='@RG\tID:mmetsp00100\tSM:VA2\tPL:illumina\tLB:mmetsp00100\tPU:MMETSP0098'

    java -jar /home/neo/data/picard/picard.jar \
        CreateSequenceDictionary \
        R=${reference} \
        O=/media/sf_docs/data/mmetsp0098/contigs.dict

    samtools faidx ${reference}

for i in 1 2 3 4 5 6
do
    sample=${sample[${i}]}
    RG=${RG[${i}]}
    bwa mem -M \
        -R ${RG} \
        -p ${reference} \
        ${ir}/${sample}.1${extension} \
        ${ir}/${sample}.2${extension} \
    > ${dir}/${sample}.sam

    java -Xmx2g -jar /home/neo/data/picard/picard.jar \
        SortSam \
        INPUT=${dir}/${sample}.sam \
        OUTPUT=${ddir}/${sample}.sorted.bam \
        SORT_ORDER=coordinate

    java -Xmx2g -jar /home/neo/data/picard/picard.jar \
        MarkDuplicates \
        INPUT=${ddir}/${sample}.sorted.bam \
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
dir=/media/sf_docs/data/mappingY3
ddir=/media/sf_docs/data/rmdupY3

extension=.trimmed.P.NY.fastq.gz
reference=/media/sf_docs/data/genomeSRv015/QPX_v015.fasta

RG[1]='@RG\tID:mmetsp0098\tSM:NY1\tPL:illumina\tLB:mmetsp0098\tPU:QPXgenomSRv015'
RG[2]='@RG\tID:mmetsp001433\tSM:NY1\tPL:illumina\tLB:mmetsp001433\tPU:QPXgenomSRv015'
RG[3]='@RG\tID:mmetsp00992\tSM:MA1\tPL:illumina\tLB:mmetsp00992\tPU:QPXgenomSRv015'
RG[4]='@RG\tID:mmetsp001002\tSM:VA1\tPL:illumina\tLB:mmetsp001002\tPU:QPXgenomSRv015'
RG[5]='@RG\tID:mmetsp0099\tSM:MA2\tPL:illumina\tLB:mmetsp0099\tPU:QPXgenomSRv015'
RG[6]='@RG\tID:mmetsp00100\tSM:VA2\tPL:illumina\tLB:mmetsp00100\tPU:QPXgenomSRv015'

    java -jar /home/neo/data/picard/picard.jar \
        CreateSequenceDictionary \
        R=${reference} \
        O=/media/sf_docs/data/genomeSRv015/QPX_v015.dict

    samtools faidx ${reference}

for i in 1 2 3 4 5 6
do
    sample=${sample[${i}]}
    RG=${RG[${i}]}
    bwa mem -M \
        -R ${RG} \
        -p ${reference} \
        ${ir}/${sample}.1${extension} \
        ${ir}/${sample}.2${extension} \
    > ${dir}/${sample}.sam

    java -Xmx2g -jar /home/neo/data/picard/picard.jar \
        SortSam \
        INPUT=${dir}/${sample}.sam \
        OUTPUT=${ddir}/${sample}.sorted.bam \
        SORT_ORDER=coordinate

    java -Xmx2g -jar /home/neo/data/picard/picard.jar \
        MarkDuplicates \
        INPUT=${ddir}/${sample}.sorted.bam \
        OUTPUT=${ddir}/${sample}.nodup.bam \
        METRICS_FILE=${ddir}/${sample}.dup.metrics \
        REMOVE_DUPLICATES=true \
        ASSUME_SORTED=true


done

java -Xmx2g -jar /home/neo/data/picard/picard.jar MarkDuplicates INPUT=./rmdupY/mmetsp00100.sorted.bam OUTPUT=./rmdupY/mmetsp00100.nodup.bam METRICS_FILE=./rmdupY/mmetsp00100.dup.metrics REMOVE_DUPLICATES=true ASSUME_SORTED=true


java -Xmx2g -jar /home/neo/data/picard/picard.jar BuildBamIndex INPUT=./rmdupY/mmetsp0098.nodup.bam

java -jar /home/neo/data/GenomeAnalysisTK.jar -T CountReads -R /media/sf_docs/data/QPX-RNA-Seq/Steve_Roberts/QPXTranscriptome_v21/QPX_transcriptome_v2orf.fasta -fixMisencodedQuals -I ./rmdupY/mmetsp0098.nodup.bam

# output to a file than extract info of nb of reads

java -jar /home/neo/data/GenomeAnalysisTK.jar -T CountReads -R /media/sf_docs/data/QPX-RNA-Seq/Steve_Roberts/QPXTranscriptome_v21/QPX_transcriptome_v2orf.fasta -fixMisencodedQuals -I ./rmdupY/mmetsp0099.nodup.bam -rf DuplicateRead

# output to a file than extract info of nb of duplicated reads

java -jar /home/neo/data/GenomeAnalysisTK.jar -T RealignerTargetCreator -R /media/sf_docs/data/QPX-RNA-Seq/Steve_Roberts/QPXTranscriptome_v21/QPX_transcriptome_v2orf.fasta -fixMisencodedQuals -I ./rmdupY/mmetsp0099.nodup.bam -o ./rmdupY/mmetsp0099.target.intervals.list

java -jar /home/neo/data/GenomeAnalysisTK.jar -T IndelRealigner -R /media/sf_docs/data/QPX-RNA-Seq/Steve_Roberts/QPXTranscriptome_v21/QPX_transcriptome_v2orf.fasta -fixMisencodedQuals -I ./rmdupY/mmetsp0099.nodup.bam -targetIntervals ./rmdupY/mmetsp0099.target.intervals.list -o ./rmdupY/mmetsp0099.realign.bam


# first call = high filters
java -jar /home/neo/data/GenomeAnalysisTK.jar -T HaplotypeCaller -R /media/sf_docs/data/QPX-RNA-Seq/Steve_Roberts/QPXTranscriptome_v21/QPX_transcriptome_v2orf.fasta -I ./rmdupY/mmetsp0099.realign.bam --genotyping_mode DISCOVERY -stand_emit_conf 30 -stand_call_conf 45 -o ./rmdupY/mmetsp0099.hi.raw.snp.vcf

# recalibration
java -jar /home/neo/data/GenomeAnalysisTK.jar -T BaseRecalibrator -R /media/sf_docs/data/QPX-RNA-Seq/Steve_Roberts/QPXTranscriptome_v21/QPX_transcriptome_v2orf.fasta -I ./rmdupY/mmetsp0099.realign.bam -knownSites ./rmdupY/mmetsp0099.hi.raw.snp.vcf -o ./rmdupY/mmetsp0099.table

# recal (2)
java -jar /home/neo/data/GenomeAnalysisTK.jar -T BaseRecalibrator -R /media/sf_docs/data/QPX-RNA-Seq/Steve_Roberts/QPXTranscriptome_v21/QPX_transcriptome_v2orf.fasta -I ./rmdupY/mmetsp0099.realign.bam -knownSites ./rmdupY/mmetsp0099.hi.raw.snp.vcf -BQSR ./rmdupY/mmetsp0099.table -o ./rmdupY/mmetsp0099.postrecal.table

# plots
java -jar /home/neo/data/GenomeAnalysisTK.jar -T AnalyzeCovariates -R /media/sf_docs/data/QPX-RNA-Seq/Steve_Roberts/QPXTranscriptome_v21/QPX_transcriptome_v2orf.fasta -before ./rmdupY/mmetsp0099.table -after ./rmdupY/mmetsp0099.postrecal.table -plots ./rmdupY/mmetsp0099.recal.plots.pdf

# apply recal
java -jar /home/neo/data/GenomeAnalysisTK.jar -T PrintReads -R /media/sf_docs/data/QPX-RNA-Seq/Steve_Roberts/QPXTranscriptome_v21/QPX_transcriptome_v2orf.fasta -I ./rmdupY/mmetsp0099.realign.bam -BQSR ./rmdupY/mmetsp0099.table -o ./rmdupY/mmetsp0099.recal.bam


#second calling
java -jar /home/neo/data/GenomeAnalysisTK.jar -T HaplotypeCaller -R /media/sf_docs/data/QPX-RNA-Seq/Steve_Roberts/QPXTranscriptome_v21/QPX_transcriptome_v2orf.fasta -I ./rmdupY/mmetsp0099.recal.bam --genotyping_mode DISCOVERY -stand_emit_conf 30 -stand_call_conf 50 -o ./rmdupY/mmetsp0099.hi.raw.snp.vcf
