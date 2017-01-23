#!/bin/bash
#SBATCH --partition=LM
#SBATCH --nodes=1
#SBATCH -t 00:30:00
#SBATCH --mem=3000GB
#SBATCH --job-name="filter"
#SBATCH --output="filter.%j.%N.out"
#SBATCH --mail-type=ALL
#SBATCH --mail-user=sleiman.bassim@stonybrook.edu

## !! ##
## filter_low_expr_transcripts.pl has been modified.
## there is a bug in the intial script
## new script was found in google group
## run chmod -R u+x after deploy
module load java
module load bowtie
module load samtools

# CHANGE___FILE ID___METHOD eXpress kallisto salmon
transcriptome=00000
method=eXpress
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
filter=$home/$version/util/filter_low_expr_transcripts.pl

#matrix=$project/trans_counts.TPM.not_cross_norm
matrix=$project/trans_counts.TMM.EXPR.matrix

time=$home/time
jobid=$transcriptome.$method.filter
start=$(date); echo "Job started at: $start" > $time/$jobid.time

cd $project
#$filter --matrix $matrix --transcripts $assembly --min_expr_any 10 --min_pct_dom_iso 50  --trinity_mode

$filter --matrix $matrix --transcripts $assembly --min_expr_any 20 --highest_iso_only --trinity_mode

end=$(date); echo "Job ended at: $end" >> $time/$jobid.time
