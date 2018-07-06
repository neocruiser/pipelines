#!/bin/bash

scratch=/cluster/projects/kridelgroup
project=$scratch/relapse

## adjust the B-statistic
_bval=$1
output=$2

cd $output/summary

## get the list of files that contain the output of Toptable from limma package
_listFiles=$(find $output -maxdepth 2 -iname "*moderated*")

## create header
echo -e "Comparison\tContrast\tID.hg.1\tChromosome\tLogFC\tAveExp\tt\tPval\tFDRadjPval\tB\tEnsembl\tSymbol\tFunction" > tmp

## mining operation on significantly selected genes
for _filename in $_listFiles; do
    paste <(printf "%s$(cat $_filename | sed -e '1d' -e 's/ /./g' | awk -vb=$_bval '{if($25>=b) print $0}' | cut -f1,3,20-25)") \
	  <(printf "%s$(cat $_filename | sed 's/ /./g' | awk -vb=$_bval '{if($25>=b) print $0}' | cut -f8 | sed 's/.\/\/./\t/g' | cut -f1-3)") | \
	sed "s/^/$(echo $(basename $_filename) | cut -f1-2 -d'.' | sed 's/\./\t/g')\t/g" \
	    >> tmp
done

grep "hg.1" $output/summary/tmp > $output/summary/summary.lmfit.fdrAdjpval.txt
rm  tmp

## get a tabulated summary of significant genes and their contrasts
head -n1 $output/summary/summary.full* > $output/summary/summary.lmfit.bval.txt
cat $output/summary/summary.full* | awk '{if($9>1 && $8 == "bval")print$0}' >> $output/summary/summary.lmfit.bval.txt


