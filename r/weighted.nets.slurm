#!/bin/bash
#SBATCH -N 1
#SBATCH -p RM
#SBATCH -t 01:00:00
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
p=4
c=2
#_P=6
#_T=01
_COR=pearson.average
_DIR=$scratch/ganglia/trinity/trinity_out_dir_00000
file=deg.eXpress*/DESeq2.eXpress.tissue.p$p.c$c*
log=$_DIR/$file/diffExpr.P1e-${p}_C${c}.matrix.log2.dat
output=$_DIR/networks/$SLURM_JOBID
mkdir -p $output

time=$home/time
jobid=$SLURM_JOBID.networks.R
start=$(date); echo "Job started at: $start" > $time/$jobid.time

cd $output
cp $log $output
cp $_DIR/contigs*p$p.c$c*id2description.txt $output
cp $home/script*/convertMatrix2graph.R $output


# Extract annotated genes from IPS tsv output
# then combine with Panther annotation.
# The joined annotation is then pushed to R from duplication removal.
# Files contain contig IDs and protein description
# node: $panther_annot is created through automated.scripts.sh
index=contigs.deseq2.p4.c2
ips_annot=$index.prot.fa.tsv
panther_annot=$ips_annot.PANprots.LEN50.EVAL10.txt
graph_annot=$index.tsv_id2description.txt

cat $_DIR/$ips_annot | sed 's/ /./g' | cut -f1,9,13 | awk '{if($2<=0.00001)print$1,$3}' | sed 's/_. / /g' | sort - | uniq -u > $_DIR/$graph_annot
paste <(awk '{print $4}' $_DIR/$panther_annot ) <(awk '{print $3}' $_DIR/$panther_annot ) | sed 's/_.\t/\t/g' | sort -k1 - | uniq -u >> $_DIR/$graph_annot


R CMD BATCH $home/script*/weighted.nets.R

#if [ -f network.graphml ]; then
#    mv network.graphml network.PVAL$p.FOLD$c.POW$_P.COR$_COR.Th$_T.$SLURM_JOBID.graphm#l
#    rm weighted.nets.Rout
#fi

end=$(date); echo "Job ended at: $end" >> $time/$jobid.time
