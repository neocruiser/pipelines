#!/bin/bash
#SBATCH -N 1
#SBATCH -p LM
#SBATCH --mem=3000GB
#SBATCH -t 24:00:00
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
_DIR=$scratch/nodule/cross-analysis
output=$_DIR/networks/$SLURM_JOBID
mkdir -p $output

# Files _CHANGE_
pa=3
ca=2
matrix_a=tissue
_DIR_a=$scratch/nodule/qpx-analysis/trinity/trinity_out_dir_A
file_a=deg.eXpress*/edgeR.eXpress.$matrix_a.p$pa.c$ca*
log_a=$_DIR_a/$file_a/diffExpr.P1e-${pa}_C${ca}.matrix.log2.dat

pb=3
cb=2
matrix_b=tissue
_DIR_b=$scratch/nodule/clam-analysis/trinity/trinity_out_dir_40p
file_b=deg.eXpress*/edgeR.eXpress.$matrix_b.p$pb.c$cb*
log_b=$_DIR_b/$file_b/diffExpr.P1e-${pb}_C${cb}.matrix.log2.dat


time=$home/time
jobid=$SLURM_JOBID.networks.R
start=$(date); echo "Job started at: $start" > $time/$jobid.time

cd $output
cp $home/script*/convertMatrix2graph.R $output

cp $log_a $output/logs.a
cp $log_b $output/logs.b

# create one id2description file
cp $_DIR_a/qpx-degs/*id2description* $output/qpx.ids
cp $_DIR_b/*id2description* $output/clam.ids
cp qpx.ids qpx.clam.ids
cat clam.ids >> qpx.clam.ids

R CMD BATCH $home/script*/weighted.nets.cross.species.R

end=$(date); echo "Job ended at: $end" >> $time/$jobid.time
