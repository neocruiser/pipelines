#!/bin/bash
#SBATCH --partition=LM
#SBATCH --nodes=1
#SBATCH --mem=3000GB
#SBATCH -t 48:00:00
#SBATCH --job-name="gills"
#SBATCH --output="abundance.%j.%N.out"
#SBATCH --mail-type=ALL
#SBATCH --mail-user=sleiman.bassim@stonybrook.edu

## !! ##
# 20h for 26 eXpress datasets
module load java
module load bowtie
module load samtools

# CHANGE___FILE ID___METHOD eXpress kallisto salmon
transcriptome=gills
method=kallisto
lib=RF

## DONT___CHANGE
version=trinityrnaseq-2.2.0
nthreads=64
pbs=$SLURM_JOBID
scratch=/pylon2/oc4ifip/bassim
home=/home/bassim
target=trinity_out_dir_$transcriptome
project=$scratch/ganglia/trinity/$target/abundance_${method}

# __CHANGE__
assembly=$scratch/ganglia/trinity/$target/Trinity.fasta
#target=parasite.$transcriptome
#project=$scratch/ganglia/blat/$target/abundance_${method}
#assembly=$scratch/ganglia/blat/$target/parasite.$n.selected.500.fa
reads=$scratch/ganglia/raw.reads
abundance=$home/$version/util/align_and_estimate_abundance.pl

time=$home/time
jobid=$transcriptome.$method.abundance
start=$(date); echo "Job started at: $start" > $time/$jobid.time

## Express uses bowtie, so its slow
## Kallisto is fast delivering short summary
if [ "$method" == eXpress ]; then
    for f in br
    do
	for i in {1..24}
	do
	    mkdir -p $project/$f$i
	    perl $abundance --transcripts $assembly \
		--SS_lib_type $lib --seqType fq \
		--left $reads/$f.${i}.R1.fastq.gz \
		--right $reads/$f.${i}.R2.fastq.gz \
		--est_method $method \
		--aln_method bowtie \
		--trinity_mode \
		--thread_count $nthreads \
		--output_dir $project/$f$i \
		--output_prefix $f$i.$method \
		--prep_reference
    done
done
    else
    for f in br
    do
	for i in {1..24}
	do
	    mkdir -p $project/$f$i
	    perl $abundance --transcripts $assembly \
		--SS_lib_type $lib --seqType fq \
		--left $reads/$f/$f.${i}.R1.fastq.gz \
		--right $reads/$f/$f.${i}.R2.fastq.gz \
		--est_method $method \
      		--trinity_mode \
		--thread_count $nthreads \
		--output_dir $project/$f$i \
		--output_prefix $f$i.$method \
		--prep_reference
    done
done
fi
end=$(date); echo "Job ended at: $end" >> $time/$jobid.time
