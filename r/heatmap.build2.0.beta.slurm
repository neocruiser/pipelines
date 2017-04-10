#!/bin/bash
#SBATCH -N 1
#SBATCH -p LM
#SBATCH --mem=3000GB
#SBATCH -t 72:00:00
#SBATCH --job-name="meta"
#SBATCH --output="R.heatmaps.out"
#SBATCH --mail-type=ALL
#SBATCH --mail-user=sleiman.bassim@stonybrook.edu

module load R

# DONT CHANGE #
scratch=/pylon2/oc4ifip/bassim
home=/home/bassim
tmp=/dev/shm/temp_$SLURM_JOBID
#mkdir -p $tmp

# Files _CHANGE_
_DIR=$scratch/tj
file=otus.genus.percentages.txt
log=logs

output=$_DIR/heatmaps/$SLURM_JOBID
mkdir -p $output

time=$home/time
jobid=$SLURM_JOBID.heatmap.R
start=$(date); echo "Job started at: $start" > $time/$jobid.time

cd $output
cp $_DIR/$file $output
mv $file $log

    R CMD BATCH $home/script*/heatmap2.0.beta.R
    rm Rplots.pdf


end=$(date); echo "Job ended at: $end" >> $time/$jobid.time
