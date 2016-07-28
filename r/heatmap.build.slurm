#!/bin/bash
#SBATCH -N 1
#SBATCH -p RM
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
p=4
c=2
_COR=pearson.average
_DIR=$scratch/ganglia/trinity/trinity_out_dir_00000
file=deg.eXpress*/DESeq2.eXpress.tissue.p$p.c$c*
log=$_DIR/$file/diffExpr.P1e-${p}_C${c}.matrix.log2.dat

output=$_DIR/heatmaps/$SLURM_JOBID
mkdir -p $output

time=$home/time
jobid=$SLURM_JOBID.networks.R
start=$(date); echo "Job started at: $start" > $time/$jobid.time

cd $output
cp $log $output
# get the node ID (contigs) from network or sub-networks
if [ -f $home/heatmaps*/*csv ]; then
    mv $home/heatmaps.data/*csv $output
fi


# get the gene expression logged FPKMs
grep -Fwf <(cat $output/*csv | sed 's/","/ /g' | cut -f7 -d' ' | sed 's/"//g' | awk 'NR>1') $log > $output/fpkms.logs
# get the sample distribution (different for each analysis)
head -n1 $log > $output/samples.logs


R CMD BATCH $home/script*/heatmap.R


end=$(date); echo "Job ended at: $end" >> $time/$jobid.time
