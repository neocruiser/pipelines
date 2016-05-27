#!/bin/bash
#SBATCH --partition=LM
#SBATCH --nodes=1
#SBATCH -t 10:00:00
#SBATCH --job-name="degEX1"
#SBATCH --output="degs.%j.%N.out"
#SBATCH --mail-type=ALL
#SBATCH --mail-user=sleiman.bassim@stonybrook.edu

module load R
module load java
module load bowtie
module load samtools

## CHANGE__PROJECT__ID__eXpress__kallisto__salmon__and fasta file
transcriptome=parasite.3
_FASTA=$transcriptome.selected.500.fa
method=eXpress
lib=RF

## DONT__CHANGE
nthreads=64
version=trinityrnaseq-2.2.0
scratch=/pylon2/oc4ifip/bassim
home=/home/bassim
pbs=$SLURM_JOBID
target=trinity_out_dir_$transcriptome
_DIR=$scratch/ganglia/trinity/$target
_FOLD=$_DIR/abundance_${method}
project=$(find $_DIR -name "abundance_*")
reads=$scratch/ganglia/blat/$transcriptome/fastq_genome_bwa
## Analyses
abundance=$home/$version/util/align_and_estimate_abundance.pl
analyze=$home/$version/Analysis/DifferentialExpression/run_DE_analysis.pl
differential=$home/$version/Analysis/DifferentialExpression/analyze_diff_expr.pl
join=$home/$version/util/abundance_estimates_to_matrix.pl
TPM=$home/$version/util/misc/count_matrix_features_given_MIN_TPM_threshold.pl
FPKM=$home/$version/util/misc/count_features_given_MIN_FPKM_threshold.pl
prefix=trans_counts

## Define function that recognize wich alignment type was used to compute contig abundance
# Get the alignment type
function type () {
    if [ ! -z "$project" ]; then
        e=$(grep -oci "express" <(echo $project))
        k=$(grep -oci "kallisto" <(echo $project))
        s=$(grep -oci "salmon" <(echo $project))
        if [ "$e" == 1 ]; then
            method=eXpress
            files=$(find $project -name "results.xprs" | paste -s -d' ')
        elif [ "$k" == 1 ]; then
            method=kallisto
            files=$(find $project -name "abundance.tsv" | paste -s -d' ')
        elif [ "$s" == 1 ]; then
            method=salmon
            files=$(find $project -name "quant.sf" | paste -s -d' ')
        fi
    fi
}


# Run abundance test if eXpress run the first configuration, else (salmon or kallisto)
# Express uses bowtie, so its slow
# Kallisto is fast delivering short summary
cd $_DIR
if [ ! -d "$project" ]; then
    if [ "$method" == eXpress ]; then
        for f in br gg
        do
	          for i in {1..24}
	          do
	              mkdir -p $_FOLD/$f$i
	              perl $abundance --transcripts $_DIR/$_FASTA --SS_lib_type $lib --seqType fq --left $reads/$f.${i}.R1.fastq --right $reads/$f.${i}.R2.fastq --est_method $method --aln_method bowtie --trinity_mode --thread_count $nthreads --output_dir $_FOLD/$f$i --output_prefix $f$i.$method --prep_reference
            done
        done
    else
        for f in br gg
        do
	          for i in {1..24}
	          do
	              mkdir -p $_FOLD/$f$i
	              perl $abundance --transcripts $_DIR/$_FASTA --SS_lib_type $lib --seqType fq --left $reads/$f/$f.${i}.R1.fastq.gz --right $reads/$f/$f.${i}.R2.fastq.gz --est_method $method --trinity_mode --thread_count $nthreads --output_dir $_FOLD/$f$i --output_prefix $f$i.$method --prep_reference
            done
        done
    fi
    type
    touch test2
else
    type
echo $files > test3
fi
