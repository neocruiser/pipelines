#!/bin/bash
#SBATCH -N 1
#SBATCH -p LM
#SBATCH --mem=3000GB
#SBATCH -t 01:00:00
#SBATCH --job-name="R.heatmaps"
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
_DIR=$scratch/nodule/qpx-analysis/trinity/trinity_out_dir_A
file=deg.eXpress*/DESeq2.eXpress.tissue.p$p.c$c*
log=$_DIR/$file/diffExpr.P1e-${p}_C${c}.matrix.log2.dat

output=$_DIR/heatmaps/$SLURM_JOBID
mkdir -p $output

time=$home/time
jobid=$SLURM_JOBID.networks.R
start=$(date); echo "Job started at: $start" > $time/$jobid.time

cd $output
cp $log $output


## the for loop will go through all clusters selected from 1 to 32

    R CMD BATCH $home/script*/heatmap2.0.beta.R
    rm Rplots.pdf

done

    end=$(date); echo "Job ended at: $end" >> $time/$jobid.time
