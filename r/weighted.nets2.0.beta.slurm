#!/bin/bash
#SBATCH -N 1
#SBATCH -p LM
#SBATCH --mem=3000GB
#SBATCH -t 48:00:00
#SBATCH --job-name="R.nets"
#SBATCH --output="R.nets.out"
#SBATCH --mail-type=ALL
#SBATCH --mail-user=sleiman.bassim@stonybrook.edu

module load java
module load R

# DONT CHANGE #
scratch=/pylon2/oc4ifip/bassim
home=/home/bassim
tmp=/dev/shm/temp_$SLURM_JOBID
#mkdir -p $tmp

# Files _CHANGE_
_DIR=$scratch/tj
file=final.qiime.local.txt
log=logs

output=$_DIR/networks/$SLURM_JOBID
mkdir -p $output

time=$home/time
jobid=$SLURM_JOBID.networks.R
start=$(date); echo "Job started at: $start" > $time/$jobid.time

cd $output
cp $_DIR/$file $output
cp $home/script*/convertMatrix2graph.R $output
mv $file $log

R CMD BATCH $home/script*/weighted.nets2.0.beta.R

end=$(date); echo "Job ended at: $end" >> $time/$jobid.time

