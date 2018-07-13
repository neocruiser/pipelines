#!/bin/bash

home=/cluster/home/sbassim
scratch=/cluster/projects/kridelgroup
project=$scratch/relapse

## adjust the B-statistic
_bval=$1
output=$2

cd $output/summary

## get the list of files that contain the output of Toptable from limma package
_listFiles=$(find $output -maxdepth 2 -iname "*moderated*")


## get gene counts significant or not per contrast after lmfit of limma
## cleaning the annotation database
## mining operation on significantly selected genes
for _index in all fdrAdjPval; do

    touch tmp.$_index

    if [ $_index == "all" ]; then
	## this task requires nothing. All in - All out

	## TASK 1
	## do not select based on bavalue or any other statistic
	## get all genes, significant or not from all contrast comparisons
	for _filename in $_listFiles; do
	    paste <(printf "%s$(cat $_filename | sed -e '1d' -e 's/ /./g' | cut -f1,3,20-25)") \
		  <(printf "%s$(cat $_filename | sed -e '1d' -e 's/ /./g' | cut -f8 | sed 's/.\/\/./\t/g' | cut -f1-3)") | \
		sed "s/^/$(echo $(basename $_filename) | cut -f1-2 -d'.' | sed 's/\./\t/g')\t/g"   >> tmp.$_index
	done
	## make sure no unmerged columns are included in output
	grep "hg.1" $output/summary/tmp.$_index > $output/summary/summary.lmfit.$_index.txt
	rm  tmp.$_index

	## create header
	sed -i "1iComparison\tContrast\tID\tChromosome\tLogFC\tAveExp\tt\tPval\tFDRadjPval\tB\tEnsembl\tSymbol\tFunction" \
	    $output/summary/summary.lmfit.$_index.txt



    elif [ $_index == "fdrAdjPval" ]; then
	## this branch requires a Bstatistics value (already set as second input)
	cd $output/summary
	
	## TASK 2
	## select based on b-value which will assure unbiased selection of adjusted pvalues.
	## get only highly significant genes based on B-statistics
	for _filename in $_listFiles; do
	    paste <(printf "%s$(cat $_filename | sed -e '1d' -e 's/ /./g' | awk -vb=$_bval '{if($25>=b) print $0}' | cut -f1,3,20-25)") \
		  <(printf "%s$(cat $_filename | sed -e '1d' -e 's/ /./g' | awk -vb=$_bval '{if($25>=b) print $0}' | cut -f8 | sed 's/.\/\/./\t/g' | cut -f1-3)") | \
		sed "s/^/$(echo $(basename $_filename) | cut -f1-2 -d'.' | sed 's/\./\t/g')\t/g" \
		    >> tmp.$_index
	done
	## make sure no unmerged columns are included in output
	grep "hg.1" $output/summary/tmp.$_index > $output/summary/summary.lmfit.$_index.txt
	rm  tmp.$_index

	## create header
	sed -i "1iComparison\tContrast\tID\tChromosome\tLogFC\tAveExp\tt\tPval\tFDRadjPval\tB\tEnsembl\tSymbol\tFunction" \
	    $output/summary/summary.lmfit.$_index.txt

    fi

done



## TASK 3
## get a tabulated summary of significant genes and their contrasts
head -n1 $output/summary/summary.full* > $output/summary/summary.lmfit.bval.txt
cat $output/summary/summary.full* | awk '{if($9>1 && $8 == "bval")print$0}' >> $output/summary/summary.lmfit.bval.txt



## TASK 4
## get ABC/GCB known expressed genes
## from lmfit in limma with fold change, Bstats, adj-pval per gene
cd $output
grep -wFf <(cat $output/summary/abc_gcb.genes.counts | cut -f1) $output/summary/summary.lmfit.all.txt  | sort - > $output/summary/abc_gcb.lmFolds.txt
sed -i "1iComparison\tContrast\tID\tChromosome\tLogFC\tAveExp\tt\tPval\tFDRadjPval\tB\tEnsembl\tSymbol\tFunction" $output/summary/abc_gcb.lmFolds.txt

## Get ABC and GCB preselected and preferentially expressed genes
## from RMA log2 quantile normalized scores per sample
grep -wFf <(cat $output/summary/abc_gcb.genes.counts | cut -f1) $output/summary/normalized.subset*  | sort - > $output/summary/abc_gcb.RMA.txt
sed -i "1iID\t$(head -n1 $output/summary/normalized.subset*)" $output/summary/abc_gcb.RMA.txt



## plotting
R CMD BATCH $home/script*/r/affymetrix.pval.distribution.R
