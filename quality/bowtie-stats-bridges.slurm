#!/bin/bash
#SBATCH --partition=LM
#SBATCH --nodes=1
#SBATCH -t 48:00:00
#SBATCH --job-name="staPar1"
#SBATCH --output="bowtie.stats.parasite1.%j.%N.out"
#SBATCH --export=ALL
#SBATCH --mail-user=sleiman.bassim@stonybrook.edu

module load java
module load samtools

## CHANGE
n=1
transcriptome=$n.124816
lib=RF

## DONT__CHANGE
version=trinityrnaseq-2.2.0
scratch=/pylon2/oc4ifip/bassim
home=/home/bassim
#workdir=${scratch}/ganglia/trinity/trinity_out_dir_$transcriptome
#target=$workdir/Trinity.fasta

workdir=${scratch}/ganglia/blat/parasite.$transcriptome
target=$workdir/parasite.$n.selected.500.fa
join=$home/$version/util/bowtie_PE_separate_then_join.pl
counts=$home/$version/util/SAM_nameSorted_to_uniq_count_stats.pl

cd $workdir

# fastq raw files COMBINED (all R1 and all R2 files separately)
#sense=$(find ${scratch}/ganglia/merged.trimmed -name "r*all.R1*q")
#antisense=$(find ${scratch}/ganglia/merged.trimmed -name "r*all.R2*q")
sense=$(find $scratch/ganglia/raw.reads -name "*R1*gz" | paste -s -d,)
antisense=$(find $scratch/ganglia/raw.reads -name "*R2*gz" | paste -s -d,)

time=$home/time
jobid=$transcriptome.bowtie.stats_$SLURM_JOBID
start=$(date); echo "Job started at: $start" > $time/$jobid.time

# Run bowtie
$join --seqType fq --SS_lib_type $lib --left $sense --right $antisense --target $target --aligner bowtie -- -p 4 --all --best --strata -m 300 >& $workdir/${jobid}_joined
# run trinity integrated stat algorithm
$counts $workdir/bowtie_out/bowtie_out.nameSorted.bam >& $workdir/${jobid}_counts

end=$(date); echo "Job ended at: $end" >> $time/$jobid.time
