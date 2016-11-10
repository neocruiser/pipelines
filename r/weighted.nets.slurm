#!/bin/bash
#SBATCH -N 1
#SBATCH -p RM
#SBATCH -t 05:00:00
#SBATCH --job-name="R.nets"
#SBATCH --output="R.nets.out"
#SBATCH --mail-type=ALL
#SBATCH --mail-user=sleiman.bassim@stonybrook.edu

module load java
#module load R

# DONT CHANGE #
scratch=/pylon2/oc4ifip/bassim
home=/home/bassim
tmp=/dev/shm/temp_$SLURM_JOBID
#mkdir -p $tmp

# Files _CHANGE_
p=4
c=2
matrix=tissue
_DIR=$scratch/ganglia/trinity/trinity_out_dir_00000
file=deg.eXpress*/DESeq2.eXpress.$matrix.p$p.c$c*
log=$_DIR/$file/diffExpr.P1e-${p}_C${c}.matrix.log2.dat
output=$_DIR/networks/$SLURM_JOBID
mkdir -p $output

time=$home/time
jobid=$SLURM_JOBID.networks.R
start=$(date); echo "Job started at: $start" > $time/$jobid.time

cd $output
cp $log $output
cp $_DIR/*id2description* $output
cp $home/script*/convertMatrix2graph.R $output


R CMD BATCH $home/script*/weighted.nets.R

end=$(date); echo "Job ended at: $end" >> $time/$jobid.time
