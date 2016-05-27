#!/bin/bash
#SBATCH -N 1
#SBATCH -p LM
#SBATCH -t 24:00:00
#SBATCH --job-name="standard"
#SBATCH --output="kraken.%j.out"
#SBATCH --mail-type=ALL
#SBATCH --mail-user=sleiman.bassim@stonybrook.edu

## FIRST____CHANGE
## fungi, virus, plasmid, archae, invertebrate
DBNAME=all

## DONT____CHANGE
# Load XSEDE modules
module load java
export CLASSPATH=$JAVA_HOME/bin:$CLASSPATH

set echo
scratch=/pylon2/oc4ifip/bassim
home=/home/bassim
location=$scratch/amelie/virome
project=${scratch}/amelie/kraken/$DBNAME
rm -r $project
mkdir -p $project
db=$scratch/db/kraken/$DBNAME
DEV=/dev/shm/$DBNAME
n=64

time=$home/time
start=$(date); echo "Job started at: $start" > $time/$SLURM_JOBID.$DBNAME.time
################
## RUN KRAKEN ##
################

cp -a $db /dev/shm

for i in {1..13}
do
    output=$project/sample.$i.$SLURM_JOBID
    result=$output/kraken
    mkdir -p $output

    kraken --threads $n --db $DEV --gzip-compressed --fastq-input --check-names --output $result.out --paired $location/*S${i}_R1* $location/*S${i}_R2*
    kraken-translate --db $DEV $result.out > $result.taxa

    kraken-report --db $DEV $result.out > $result.report
    kraken-mpa-report --db $DEV $result.out > $result.mpa.report
done

rm -r $DEV

# get the proportion of classified and unclassified reads
grep "class" $home/kraken.$SLURM_JOBID.out > $project/output.${DBNAME}_log
# concatenate classified reads/phyla accross all samples
cd $project
out=summary.$DBNAME.txt
final=summary.$DBNAME.txt
rm $final

for d in sample.{1..13}.*; do
    for f in $d/kraken.mpa.report; do
        abbrv=$(echo $d | cut -f2 -d '.')
        sed -e 's/$/\t'$abbrv'/g' <(cat $f) >> $out
    done
done


end=$(date); echo "Job ended at: $end" >> $time/$SLURM_JOBID.$DBNAME.time
