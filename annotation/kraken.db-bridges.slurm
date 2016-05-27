#!/bin/bash
#SBATCH -N 1
#SBATCH -p LM
#SBATCH -t 48:00:00
#SBATCH --job-name="build"
#SBATCH --output="kraken.%j.%N.out"
#SBATCH --mail-type=ALL
#SBATCH --mail-user=sleiman.bassim@stonybrook.edu

## database names: bacteria, archae, virus, plasmid, invertebrate; or standard (all)
DBNAME=all
scratch=/pylon2/oc4ifip/bassim
home=/home/bassim
db=$scratch/db/kraken
DEV=/dev/shm/$DBNAME
lib=$DEV/library

rm -r $DEV
cp -a $db/$DBNAME /dev/shm
#mkdir -p $lib && cd $lib
#kraken-build --download-taxonomy --db $DEV
#wget ftp://ftp.ncbi.nlm.nih.gov/refseq/release/$DBNAME/* >& $home/kraken.build.out_log.$SLURM_JOBID
for i in archae bacteria fungi invertebrate plasmid virus
do
cd $lib/$i
gunzip -v *.gz >& $home/kraken.build.out_log.$SLURM_JOBID
done
cd $lib
chmod -R u+x $DEV
kraken-build --threads 64 --build --jellyfish-hash-size 12800M --db $DEV >& $home/kraken.build.out_log.$SLURM_JOBID
kraken-build --clean --db $DEV
mv $DEV $db/$DBNAME.$SLURM_JOBID

#--max-db-size 50
