#!/bin/bash
#SBATCH --partition=LM
#SBATCH --nodes=1
#SBATCH --mem=3000GB
#SBATCH -t 05:00:00
#SBATCH --job-name="fastq"
#SBATCH --output="fastq.%j.%N.out"
#SBATCH --mail-type=ALL
#SBATCH --mail-user=sleiman.bassim@stonybrook.edu

module load java
module load picard
module load samtools

## DONT__CHANGE
scratch=/pylon2/oc4ifip/bassim
home=/home/bassim
version=$home/picard*/picard-2.8.2.jar
pbs=$SLURM_JOBID



for i in A; do
    for k in {1..3}; do
        target=$i$k.nodup
        _DIR=$scratch/nodule/clam_reads/$target
        _IR=$scratch/nodule/rmdup_QPX/$target

mkdir -p $_DIR

        # extract only the mapped reads in bam (-b)
        # remove (-F) the reads that do not map (4)
        samtools view -b -F2 $_IR.bam > $_DIR.clam.bam
#        samtools view -b -F 4 $_IR.bam > $_DIR.mapped.bam

java -Xmx1500g -jar $version SamToFastq I=$_DIR.clam.bam FASTQ=$_DIR.R1.fastq SECOND_END_FASTQ=$_DIR.R2.fastq

for j in {1..2}; do
    gzip $_DIR.R$j.fastq
done


    done
done
