#!/bin/bash
#SBATCH --partition=RM
#SBATCH --nodes=1
#SBATCH -t 1:00:00
#SBATCH --job-name="degEX4"
#SBATCH --output="degs.%j.%N.out"
#SBATCH --mail-type=ALL
#SBATCH --mail-user=sleiman.bassim@stonybrook.edu

module load R
module load java
module load bowtie
module load samtools

## CHANGE__PROJECT__ID__eXpress__kallisto__salmon__and fasta file
transcriptome=parasite.4
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


# Run abundance test if eXpress run the first configuration, else (salmon or kallisto)
# Express uses bowtie, so its slow
# Kallisto is fast delivering short summary
cd $_DIR
if [ ! -d "$project" ]; then
    if [ "$method" == eXpress ]; then
        for f in br gg; do
	          for i in {1..24}; do
	              mkdir -p $_FOLD/$f$i
	              perl $abundance --transcripts $_DIR/$_FASTA --SS_lib_type $lib --seqType fq --left $reads/$f.${i}.R1.fastq --right $reads/$f.${i}.R2.fastq --est_method $method --aln_method bowtie --trinity_mode --thread_count $nthreads --output_dir $_FOLD/$f$i --output_prefix $f$i.$method --prep_reference
            done
        done

    elif [ "$method" == kallisto || "$method" == salmon ]; then
        for f in br gg; do
	          for i in {1..24}; do
	              mkdir -p $_FOLD/$f$i
	              perl $abundance --transcripts $_DIR/$_FASTA --SS_lib_type $lib --seqType fq --left $reads/$f/$f.${i}.R1.fastq.gz --right $reads/$f/$f.${i}.R2.fastq.gz --est_method $method --trinity_mode --thread_count $nthreads --output_dir $_FOLD/$f$i --output_prefix $f$i.$method --prep_reference
            done
        done
    fi
fi

# Get the alignment type that was used to compute contig abundance
project=$(find $_DIR -name "abundance_*")
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


# Join gene counts between samples
cd $project
if [ ! -f $prefix.TPM.not_cross_norm.counts_by_min_TPM_$method ]; then
    perl $join --est_method $method --out_prefix $prefix --name_sample_by_basedir $files
    # merge matrices accross samples to get shared TPM scores
    perl $TPM $prefix.TPM.not_cross_norm > $prefix.TPM.not_cross_norm.counts_by_min_TPM_$method
    # merge matrices accross samples to get shared FPKM scores
    #perl $FPKM $prefix.TPM.not_cross_norm > $prefix.TPM.not_cross_norm.counts_by_min_FPKM_$method
    else
    echo "Matrices have been already compiled"
fi

## Get the length of contigs that contribute to 50% of the assembled trxome (N50)
# Also get the transcript abundance for each Nx
$home/$version/util/misc/contig_ExN50_statistic.pl $project/$prefix.TMM.EXPR.matrix $_DIR/$_FASTA > $_FOLD/ExN50.TRX-$transcriptome.$pbs.stats

# Get discriptive stats of the assembled trxome
$home/$version/util/TrinityStats.pl $_DIR/$_FASTA > $_FOLD/descriptive.TRX-$transcriptome.$pbs.stats

## Choose__matrices [i]
jobid[1]=tissue
jobid[2]=tissue-diet
jobid[3]=tissue-br
jobid[4]=tissue-gg
jobid[5]=tissue-br-females
jobid[6]=tissue-gg-females
jobid[7]=tissue-br-bucephalus

## Get differentially expressed genes
dir=$_DIR/deg.$method.$pbs
mkdir -p $dir
for align in $method
do
    for Rpack in DESeq2 edgeR voom
    do
	      for i in {1..6}
	      do
	          for pval in {1..6}
	          do
		            for cfold in {1..2}
		            do
	    jobid=${jobid[${i}]}
	    matrix=$scratch/ganglia/trinity/matrix/$jobid.txt
	    contrast=$scratch/ganglia/trinity/matrix/contrast.$jobid

	    cd $project
	    $analyze --matrix $project/trans_counts.counts.matrix --method $Rpack --samples_file $matrix --output $dir/$Rpack.$align.$jobid.p$pval.c$cfold.$pbs --contrasts $contrast

	    cd $dir/$Rpack.$align.$jobid.p$pval.c$cfold.$pbs
	    $differential --matrix $project/trans_counts.TMM.EXPR.matrix -P 1e-$pval -C $cfold --samples $matrix
	              done
	          done
	      done
    done
done

## backup abundance folder and cleanup
cd $_DIR
ls *ebwt | xargs rm
mv $project $_DIR/eXpress.${pbs}_abundance
ls *gene_trans_map | xargs rm
ls *bowtie.ok | xargs rm

# Create a table for the number of differentially expressed genes
#tissue tissue-diet tissue-br tissue-gg tissue-br-females tissue-gg-females tissue-br-bucephalus
if [ -d "$dir" ]; then
cd $dir

for m in $method
do
    for i in DESeq2 edgeR voom
    do
        for t in {1..6}
        do
            for p in {1..6}
            do
                for c in {1..2}
                do
                    for f in $i*$m*${jobid[${t}]}.p$p.c$c*
                    do
temp=summary.txt
final=summary.$method.$pbs.txt
summary=$scratch/ganglia/trinity/$target/$final
rm $final

# Get the number of genes per abundance test
cat ${f}/diffExpr*matrix.log2.dat | cut -f 1 >> raw.$m.${jobid[${t}]}.$p.$c
# count number of all and unique differentially expressed genes
all=$(grep "^TRINITY" raw.$m.${jobid[${t}]}.$p.$c | wc -l)
uniq=$(grep "^TRINITY" raw.$m.${jobid[${t}]}.$p.$c | sort - | uniq | wc -l)
paste <(printf "%s\n" "$f") <(printf "%s\n" "$all") <(printf "%s\n" "$uniq") >> $temp
# column names; trandform to tabulated format
cat $temp | sed -e 's/\./\t/g' -e '/\*/d' >> $summary
rm raw.$m.${jobid[${t}]}.$p.$c $temp
                    done
                done
            done
        done
    done
done

else
    echo "A differential expression gene test must be executed first"
    scancel $pbs
fi

