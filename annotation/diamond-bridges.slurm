#!/bin/bash
#SBATCH --partition=LM
#SBATCH --nodes=1
#SBATCH --mem=3000GB
#SBATCH -t 10:00:00
#SBATCH --job-name="diamond"
#SBATCH --output="diamond.%j.%N.out"
#SBATCH --mail-type=ALL
#SBATCH --mail-user=sleiman.bassim@stonybrook.edu

module load blast
module load diamond
module load boost

# Files _CHANGE_
## OPTIONS: Pfam-A nr
db=00000
blast=blastx
maxseq=10
transcriptome=00000
## file is either txt or csv
file=fa
index=lectins.marta

# DONT CHANGE #
###############
scratch=/pylon2/oc4ifip/bassim
backup=/pylon1/oc4ifip/bassim
home=/home/bassim
pbs=$SLURM_JOBID
project=$scratch/ganglia/trinity/trinity_out_dir_$transcriptome
query=$project/$index.fa
output=$project/diamond/diamond.$index.$db.$blast.$transcriptome.$pbs
tmp=/dev/shm/tmp_$pbs

mkdir -p $tmp $project/diamond

## extract Fasta sequences from selected networks

if [ $file == "txt" ]; then
    # if IDs are in txt file (1 column)
    cat $project/Trinity.fasta | sed 's/.len.*$//g' | perl -ne 'if(/^>(\S+)/){$c=$i{$1}}$c?print:chomp;$i{$_}=1 if @ARGV' <(cat $project/$index.txt) - > $query
elif [ $file == "csv" ]; then
    # If the IDs are stored in CSV file (1 column) imported from excel or cytoscape
    cat $project/Trinity.fasta | sed 's/.len.*$//g' | perl -ne 'if(/^>(\S+)/){$c=$i{$1}}$c?print:chomp;$i{$_}=1 if @ARGV' <(cat $home/$index.csv | sed -e '1 d' -e 's/","/\t/g' -e 's/"//g' | cut -f6 | sort - | uniq) - > $query
    else
    # If only fasta file is used -- nucleotide sequence
    query=$project/$index.fasta
fi


# blast libraries
#diamond=$home/diamond.0.7/diamond
#export PATH="$PATH:/gpfs/home/ballam/ncbi-blast-2.2.31+/bin"
export BLASTDB="$backup/db/$db"

## Full blastx
time=$home/time
jobid=diamond.$db.$blast
start=$(date); echo "Job started at: $start" > $time/$jobid.$pbs.time

cd $backup/db/$db

diamond $blast -d $db -q $query -a $output -t $tmp -k $maxseq
diamond view -a $output.daa -o $output.txt

rm -r $tmp

end=$(date); echo "Job ended at: $end" >> $time/$jobid.$pbs.time
