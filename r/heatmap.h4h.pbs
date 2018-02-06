#!/bin/bash
#SBATCH -N 1
#SBATCH -p LM
#SBATCH --mem=3000GB
#SBATCH -t 10:00:00
#SBATCH --job-name="R.heatmaps"
#SBATCH --output="R.heatmaps.out"
#SBATCH --mail-type=ALL
#SBATCH --mail-user=sleiman.bassim@stonybrook.edu

set -eux

module load R

## choose between txt and csv
_exe=dat

# DONT CHANGE #
scratch=/pylon2/oc4ifip/bassim
home=/home/bassim
tmp=/dev/shm/temp_$SLURM_JOBID
#mkdir -p $tmp

# Files _CHANGE_
p=3
c=2
matrix=tissue
_DIR=$scratch/nodule/clam-analysis/trinity/trinity_out_dir_40p
file=deg.eXpress*90*/edgeR.eXpress.$matrix.p$p.c$c*
log=$_DIR/$file/diffExpr.P1e-${p}_C${c}.matrix.log2.dat

output=$_DIR/heatmaps/$SLURM_JOBID
mkdir -p $output

time=$home/time
jobid=$SLURM_JOBID.networks.R
start=$(date); echo "Job started at: $start" > $time/$jobid.time

cd $output
cp $log $output


## the for loop will go through all clusters selected from 1 to 32
for i in 1; do
    SUB_NET=module.30520.$i.$_exe

    # get the node ID (contigs) from network or sub-networks
    if [ -f "$home/heatmaps.data/$SUB_NET" ]; then
        cp $home/heatmaps.data/$SUB_NET $output

        if [ "$_exe" == csv ]; then
            # get the gene expression logged TPMs
            grep -Fwf <(cat $output/$SUB_NET | sed 's/","/ /g' | cut -f6 -d' ' | sed 's/"//g' | awk 'NR>1') $log > $output/fpkms.logs
            # get the sample distribution (different for each analysis)
            head -n1 $log > $output/samples.logs

        elif [ "$_exe" == txt ]; then
            grep -Fwf <(cat $output/$SUB_NET) $log > $output/fpkms.logs
            head -n1 $log > $output/samples.logs

        fi

    elif [ "$_exe" == dat ]; then
            rm $output/$SUB_NET
            cat $log | sed '1d' > $output/fpkms.logs
            head -n1 $log > $output/samples.logs

    else
        echo "!!!ERROR!!! No contig IDs were found. Export $SUB_NET from Cytoscape and place it in Bridges at ~/heatmap.data/. Cheers"
fi



    R CMD BATCH $home/script*/heatmap.R
    rm Rplots.pdf
    mv heatmap.pdf heatmap.$i.pdf
    mv bootstrap.pdf bootstrap.$i.pdf
done

    end=$(date); echo "Job ended at: $end" >> $time/$jobid.time
