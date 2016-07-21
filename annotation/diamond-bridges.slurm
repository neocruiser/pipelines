#!/bin/bash
#SBATCH --partition=RM
#SBATCH --nodes=1
#SBATCH -t 01:00:00
#SBATCH --job-name="diamond"
#SBATCH --output="diamond.%j.%N.out"
#SBATCH --mail-type=ALL
#SBATCH --mail-user=sleiman.bassim@stonybrook.edu

module load blast
module load diamond
module load boost

# Files _CHANGE_
## OPTIONS: Pfam-A nr
db=nr
blast=blastx
maxseq=10
transcriptome=00000
index=graph.pow8.thd05.gen115

# DONT CHANGE #
###############
scratch=/pylon2/oc4ifip/bassim
home=/home/bassim
pbs=$SLURM_JOBID
project=$scratch/ganglia/trinity/trinity_out_dir_$transcriptome
query=$project/$index.fa
output=$project/diamond/diamond.$index.$db.$blast.$transcriptome.$pbs
tmp=/dev/shm/tmp_$pbs

mkdir -p $tmp $project/diamond

## extract Fasta sequences from selected networks
cat $project/Trinity.fasta | sed 's/.len.*$//g' | perl -ne 'if(/^>(\S+)/){$c=$i{$1}}$c?print:chomp;$i{$_}=1 if @ARGV' <(cat $home/$index.txt | cut -f13 | awk 'NR>1') - > $query

# blast libraries
#diamond=$home/diamond.0.7/diamond
#export PATH="$PATH:/gpfs/home/ballam/ncbi-blast-2.2.31+/bin"
export BLASTDB="/pylon2/oc4ifip/bassim/db/$db"

## Full blastx
time=$home/time
jobid=diamond.$db.$blast
start=$(date); echo "Job started at: $start" > $time/$jobid.$pbs.time

cd $scratch/db/$db

diamond $blast -d $db -q $query -a $output -t $tmp -k $maxseq
diamond view -a $output.daa -o $output.txt

rm -r $tmp

end=$(date); echo "Job ended at: $end" >> $time/$jobid.$pbs.time
