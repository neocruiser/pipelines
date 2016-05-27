#!/bin/bash
#SBATCH --partition=LM
#SBATCH --nodes=1
#SBATCH -t 48:00:00
#SBATCH --job-name="detonate"
#SBATCH --output="detonate.%j.%N.out"
#SBATCH --export=ALL
#SBATCH --mail-user=sleiman.bassim@stonybrook.edu


### SCRIPT NOT COMPLETE

module load R
module load samtools
module load detonate
module load bowtie

# CHANGE
transcriptome=111380
length=500
nthread=64

# DONT__CHANGE
scratch=/pylon2/oc4ifip/bassim
home=/home/bassim
workdir=$scratch/ganglia/trinity/trinity_out_dir_$transcriptome
output=$workdir/detonate
target=$workdir/Trinity.fasta

mkdir -p $workdir $detonate
cd $workdir

# fastq raw files COMBINED (all R1 and all R2 files separately)                                                                                                                                       
sense=$(find ${scratch}/ganglia/merged.trimmed -name "r*all.R1*q")
antisense=$(find ${scratch}/ganglia/merged.trimmed -name "r*all.R2*q")
#sense=$(find $scratch/ganglia/raw.reads -name "*R1*gz" | paste -s -d,)
#antisense=$(find $scratch/ganglia/raw.reads -name "*R2*gz" | paste -s -d,)

time=$home/time
jobid=$transcriptome.detonate_$SLURM_JOBID
start=$(date); echo "Job started at: $start" > $time/$jobid.time

$home/detonate-1.10/rsem-eval/rsem-eval-calculate-score \
--seed 3471609 \
--samtools-sort-mem $nthread \
--bowtie \
--strand-specific \
--num-threads $nthread \
--time \
--paired-end \
$sense $antisense $target \
gg60 \
$length
