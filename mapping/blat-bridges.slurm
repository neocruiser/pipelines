#!/bin/bash
#SBATCH --partition=LM
#SBATCH --nodes=1
#SBATCH -t 48:00:00
#SBATCH --job-name="blat"
#SBATCH --output="blat.%j.%N.out" 
#SBATCH --export=ALL
#SBATCH --mail-user=sleiman.bassim@stonybrook.edu

## choose genome (or parasite)
genome=1
transcriptome=00000
ft=dnax
fq=rnax

## DONT___CHANGE
scratch=/pylon2/oc4ifip/bassim
home=/home/bassim
blat=$home/blat/blat
bit=$home/blat/faToTwoBit
genome_fasta=$(find $scratch/ganglia/genomes/parasite.$genome -name "*fna")
target=$genome_fasta.2bit
assembly=$scratch/ganglia/trinity/trinity_out_dir_$transcriptome/Trinity.fasta
dir=$scratch/ganglia/blat/parasite.$genome.$SLURM_JOBID
output=$dir/output.trx.$transcriptome.$SLURM_JOBID.psl

time=$home/time
jobid=$SLURM_JOBID.blat.trx.$transcriptome.genome.$genome
start=$(date); echo "Job started at: $start" > $time/$jobid.time

mkdir -p $dir
if [ ! -f $target ]; then
    $bit $genome_fasta $target
    $blat -t=$ft -q=$fq $target $assembly $output
else
    $blat -t=$ft -q=$fq $target $assembly $output
fi
end=$(date); echo "Job ended at: $end" >> $time/$jobid.time


