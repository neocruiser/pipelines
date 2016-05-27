#!/bin/bash
#SBATCH -N 1
#SBATCH -p LM
#SBATCH -t 2:00:00
#SBATCH --job-name="IPS4"
#SBATCH --output="ips.parasite.%j.%N.out"
#SBATCH --mail-type=ALL
#SBATCH --mail-user=sleiman.bassim@stonybrook.edu

module load java

# DONT CHANGE #
scratch=/pylon2/oc4ifip/bassim
home=/home/bassim
#ips=$scratch/db/ips/64/interproscan-5.18-57.0
ips=$scratch/db/ips/64/interproscan-5.17-56.0
tmp=/dev/shm/temp_$SLURM_JOBID

# ChaNGE
#peptides=$(find $scratch -maxdepth 2 -iname "prote*seq*2*fa")
#output=$(dirname $peptides)
mkdir -p $tmp

# Files _CHANGE_
p=3
c=2
transcriptome=4
_DIR=$scratch/ganglia/trinity/trinity_out_dir_parasite.$transcriptome

file=deg*248/DESeq2.eXpress.tissue.p$p.c$c*
log=$_DIR/$file/diffExpr.P1e-${p}_C${c}.matrix.log2.dat
assembled=$_DIR/Trinity-GG.fasta
contigs=$tmp/contigs.fa
peptides=$tmp/peptides-GG.fa
output=$_DIR

time=$home/time
jobid=$SLURM_JOBID.ips
start=$(date); echo "Job started at: $start" > $time/$jobid.time

# Get gene ids (only differentially expressed)
# Get gene sequences (1)
cat $assembled | sed 's/.len.*$//g' | perl -ne 'if(/^>(\S+)/){$c=$i{$1}}$c?print:chomp;$i{$_}=1 if @ARGV' <(cat $log | cut -f1 | grep "^TRINITY" | sort - | uniq) - > $contigs
# Translate
transeq $contigs $peptides -frame=F --clean=yes
# Run interpro Scans
# Dont run HAMAP PRINTS
$ips/interproscan.sh -t p -i $peptides -appl ProDom,TIGRFAM,SUPERFAMILY,PIRSF,Gene3D,Pfam,SMART,PROSITEPROFILES,PROSITEPATTERNS,COILS,PANTHER,HAMAP,PRINTS \
    -iprlookup -goterms --pathways -f TSV --tempdir $tmp -d $output

rm -r $tmp
end=$(date); echo "Job ended at: $end" >> $time/$jobid.time

#    -mode cluster -clusterrunid test \
#    -f TSV, SVG, GFF3, XML, HTML \
#    -appl ProDom,PANTHER,TIGRFAM,SUPERFAMILY,PIRSF,Gene3D,Pfam,SMART,PROSITEPROFILES,PROSITEPATTERNS,COILS \
