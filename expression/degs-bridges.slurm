#!/bin/bash
#SBATCH --partition=LM
#SBATCH --nodes=1
#SBATCH --mem=3000GB
#SBATCH -t 150:00:00
#SBATCH --job-name="DEGtest"
#SBATCH --output="degs.%j.%N.out"
#SBATCH --mail-type=ALL
#SBATCH --mail-user=sleiman.bassim@stonybrook.edu

module load R
module load java
module load bowtie
module load samtools

## CHANGE__PROJECT__ID__eXpress__kallisto__salmon__and fasta file
transcriptome=40p
_FASTA=clam_raw_trinity.fa
lib=RF

## DONT__CHANGE
nthreads=64
version=trinityrnaseq-2.2.0
scratch=/pylon2/oc4ifip/bassim
home=/home/bassim
pbs=$SLURM_JOBID
target=trinity_out_dir_$transcriptome
### CHANGE ###
reads=$scratch/nodule/raw.reads
_DIR=$scratch/nodule/clam-analysis/trinity/$target
## Analyses
abundance=$home/$version/util/align_and_estimate_abundance.pl
analyze=$home/$version/Analysis/DifferentialExpression/run_DE_analysis.pl
differential=$home/$version/Analysis/DifferentialExpression/analyze_diff_expr.pl
join=$home/$version/util/abundance_estimates_to_matrix.pl
TPM=$home/$version/util/misc/count_matrix_features_given_MIN_TPM_threshold.pl
FPKM=$home/$version/util/misc/count_features_given_MIN_FPKM_threshold.pl
prefix=trans_counts

    method[1]=eXpress
    method[2]=kallisto
    method[3]=salmon

mkdir -p $_DIR
project=$(find $_DIR -name "abundance_*")

# choose one method or more
    for w in 1 ; do
	  method=${method[${w}]}
_FOLD=$_DIR/abundance_${method}

# Run abundance test if eXpress run the first configuration, else (salmon or kallisto)
# Express uses bowtie, so its slow
# Kallisto is fast delivering short summary
# salmon supposedly reduces FDRs
cp $scratch/nodule/clam_trxome/$_FASTA $_DIR
cd $_DIR
if [ ! -d "$project" ]; then

    if [ "$method" == eXpress ]; then
                for f in A B C; do
	              for i in {1..3}; do
	                  mkdir -p $_FOLD/$f$i
	                  perl $abundance --transcripts $_DIR/$_FASTA --SS_lib_type $lib --seqType fq --left $reads/$f${i}R1.fastq.gz --right $reads/$f${i}R2.fastq.gz --est_method $method --aln_method bowtie --trinity_mode --thread_count $nthreads --output_dir $_FOLD/$f$i --output_prefix $f$i.$method --prep_reference
                done
            done

    elif [ "$method" == kallisto ]; then
        for f in A B C; do
            for i in {1..3}; do
  	                  mkdir -p $_FOLD/$f$i
	                  perl $abundance --transcripts $_DIR/$_FASTA --SS_lib_type $lib --seqType fq --left $reads/$f${i}R1.fastq.gz --right $reads/$f${i}R2.fastq.gz --est_method $method --trinity_mode --thread_count $nthreads --output_dir $_FOLD/$f$i --output_prefix $f$i.$method --prep_reference
                done
            done

    elif [ "$method" == salmon ]; then
            for f in A B C; do
	              for i in {1..3}; do
	                  mkdir -p $_FOLD/$f$i
	                  perl $abundance --transcripts $_DIR/$_FASTA --SS_lib_type $lib --seqType fq --left $reads/$f${i}R1.fastq.gz --right $reads/$f${i}R2.fastq.gz --est_method $method --trinity_mode --thread_count $nthreads --output_dir $_FOLD/$f$i --output_prefix $f$i.$method --prep_reference
                done
            done

            else
            echo "Must input eXpress, kallisto or salmon."
fi
fi

touch $_FOLD/$method.abundance.ok

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
    # if TMM is not generated or empty recheck the transcriptome
    # if reads are not those used to assemble the transcriptome
    # edgeR throws an error for TMM
cd $project
if [ ! -f $prefix.TPM.not_cross_norm.counts_by_min_TPM_$method ]; then
    perl $join --est_method $method --cross_sample_norm TMM --out_prefix $prefix --name_sample_by_basedir $files
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
jobid[8]=ab
jobid[9]=ac

## Get differentially expressed genes
dir=$_DIR/deg.$method.$pbs
mkdir -p $dir

for align in $method
do
    for Rpack in edgeR
    do
	      for o in 8 9
	      do
	          for pval in 3
	          do
		            for cfold in 2
		            do
	    jobid=${jobid[${o}]}
	    matrix=$scratch/nodule/matrix/$jobid.txt
	    contrast=$scratch/nodule/matrix/contrast.$jobid

	    cd $project
#	    $analyze --matrix $project/trans_counts.counts.matrix --method $Rpack --samples_file $matrix --output $dir/$Rpack.$align.$jobid.p$pval.c$cfold.$pbs

 	    $analyze --matrix $project/trans_counts.counts.matrix --method $Rpack --samples_file $matrix --output $dir/$Rpack.$align.$jobid.p$pval.c$cfold.$pbs --contrasts $contrast

	    cd $dir/$Rpack.$align.$jobid.p$pval.c$cfold.$pbs
	    $differential --matrix $project/trans_counts.TMM.EXPR.matrix -P 1e-$pval -C $cfold --samples $matrix
	              done
	          done
	      done
    done
done

touch $_FOLD/$method.express.ok

## backup abundance folder and cleanup
cd $_DIR
ls *ebwt | xargs rm
mv $_FOLD $_DIR/$method.${pbs}_abundance
ls *gene_trans_map | xargs rm
ls *idx | xargs rm
ls *bowtie.ok | xargs rm

# Create a table for the number of differentially expressed genes
#tissue tissue-diet tissue-br tissue-gg tissue-br-females tissue-gg-females tissue-br-bucephalus
if [ -d "$dir" ]; then
cd $dir

for m in $method
do
    for g in edgeR
    do
        for t in 8 9
        do
            for p in 3
            do
                for c in 2
                do
                    for f in $g*$m*${jobid[${t}]}.p$p.c$c*
                    do
temp=summary.txt
summary=$dir/summary.$method.$pbs.txt

# Get the number of genes per abundance test
all=$(cat ${f}/diffExpr*matrix.log2.dat | sed '1d' | wc -l)

# count number of all and unique differentially expressed genes
uniq=$(cat ${f}/diffExpr*matrix.log2.dat | sed '1d' | sort - | uniq | wc -l)
paste <(printf "%s\n" "$f") <(printf "%s\n" "$all") <(printf "%s\n" "$uniq") >> $temp

# column names; transform to tabulated format
cat $temp | sed 's/\./\t/g' >> $summary
rm $temp

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

done
