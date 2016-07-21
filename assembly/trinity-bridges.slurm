#!/bin/bash
#SBATCH --partition=LM
#SBATCH --nodes=1
#SBATCH --mem=3000GB
#SBATCH -t 24:00:00
#SBATCH --job-name="trin98799"
#SBATCH --output="trinity.%j.%N.out"
#SBATCH --mail-type=ALL
#SBATCH --mail-user=sleiman.bassim@stonybrook.edu

### Would finish in 40 hours for 400 million reads and 65h for 1 billion reads

module load java
module load bowtie
module load samtools

## direct temp files to scratch
#export TMPDIR=$LOCAL

# xsede directories
version=trinityrnaseq-2.2.0
scratch=/pylon2/oc4ifip/bassim
home=/home/bassim
backupdir=${scratch}/ganglia/trinity/trinity_out_dir_${SLURM_JOBID}
workdir=/dev/shm/trinity_out_dir_${SLURM_JOBID}
mkdir -p $workdir $backupdir
cd $workdir

# fastq raw files COMBINED (all R1 and all R2 files separately)
#sense=$(find ${scratch}/ganglia/merged.trimmed -name "r*all.R1*q")
#antisense=$(find ${scratch}/ganglia/merged.trimmed -name "r*all.R2*q")
sense=$(find $scratch/ganglia/raw.reads -name "*R1*gz" | paste -s -d,)
antisense=$(find $scratch/ganglia/raw.reads -name "*R2*gz" | paste -s -d,)

#############
# TRINITY
#############
lib=RF
JMb=3000G
bthreads=32
nthreads=32
heap=350G
gc=6
init=2G

#--normalize_by_read_set
#--normalize_max_read_cov 50
#--quality_trimming_params "LEADING:5 TRAILING:5 SLIDINGWINDOW:4:15 MINLEN:36"

$home/$version/Trinity --seqType fq --SS_lib_type $lib --left ${sense} --right ${antisense}  --normalize_max_read_cov 50 --quality_trimming_params "LEADING:5 TRAILING:5 SLIDINGWINDOW:4:15 MINLEN:36" --max_memory ${JMb} --CPU $nthreads --bflyCPU $bthreads --bflyHeapSpaceMax $heap --bflyHeapSpaceInit $init --bflyGCThreads $gc --min_contig_length 200 --output $workdir >& ${home}/trinity.${SLURM_JOBID}_output.log

mv $workdir/Trinity.fasta $backupdir
mv $workdir/Trinity.timing $backupdir
cd $workdir
perl -e 'for(<*>){((stat)[9]<(unlink))}'
rm -rf $workdir
