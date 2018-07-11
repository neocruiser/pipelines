#!/bin/bash

## Create a summary of number of genes
## per condition
## and per statistical parameter
# Choose expression parameters: B-stats, adjusted Pval, average expression, logFC
#_B=$1
_aP=$2
_aE=$3
_foldH=$4
_foldL=$5

## server information
pbs=$6
output=$7

# index all files to be summarized
summary=$output/summary/summary.full.$pbs.txt
_TMP=$output/summary/summary.tmp

listFiles=$(find $output -maxdepth 3 -iname "*moderated*")

_BVAL=bval.tmp
_ADJPVAL=adjpval.tmp
_AVGEX=avgex.tmp

touch $_TMP

# Create column names of the summary file
function scaledSummary () {
    ## sort data into columns
    ## data is generated from expression scores and gene lists
    local FILE=$1
    local CATGR=$2
    local BTHRESHOLD=$3
    ## get the percentage
    ## percentage of B-statsitics (pbval) based on the 15% of genes selected
    ## percentange of adjusted Pval and average expression based on the total gene counts w/ only significant B-statistics
    pbval=$(printf "%.2f%%" $(echo "scale=3;($bval/$total_genes)*100" | bc -l))
    padjpval=$(printf "%.2f%%" $(echo "scale=3;($adjpval/$bval)*100" | bc -l))
    pavgex=$(printf "%.2f%%" $(echo "scale=3;($avgex/$bval)*100" | bc -l))

    
    for param in bval adjpval avgex; do
	paste <(printf "%s\n" "$design") \
	      <(printf "%s\n" "$BTHRESHOLD") \
	      <(printf "%s\n" "$_aP") \
	      <(printf "%s\n" "$CATGR")  \
	      <(printf "%s\n" "$param") \
	      <(printf "%s\n" "${!param}") \
	      <(percent=p$param; printf "%s\n" "${!percent}") \
	      <(min=${param}_min; printf "%s\n" "${!min}") \
	      <(max=${param}_max; printf "%s\n" "${!max}") \
	      <(avg=${param}_avg; printf "%s\n" "${!avg}") \
	      <(median=${param}_median; printf "%s\n" "${!median}") \
	      >> $FILE
    done
}


function cleaning () {
    ## organize the coliumns in the summary file
    ## first will add column headers
    ## second will sort depending on the statistical contrasts
    ## finally will discard tmp files
    local OUTPUT=$1
    local FILE=$2
    
    echo -e "Design\tModel\tRpack\tStatistics\tBthreshold\tadjPval\tCategory\tParameter\tTranscripts\tGeneCount\tMinimumB\tMaximumB\tMeanB\tMedianB" > $FILE
    cat $OUTPUT | sort -k2 >> $FILE
    rm *tmp
    rm $OUTPUT
}

# Counting number of genes representing different thresholds of the B-statistics
# higher B-statistics scores higher is the confidence that each gene is significantly expressed
# Up and down regulated genes are also counted
# the R code above will in addition to preprocessing the data and extracting expression scores, create venn diagrams of gene expressions

for _B in -2 -1 0 0.5 1 1.5; do
    for i in total; do
	##    for i in total up down; do
	## bug inside the up/down summarization code
	for files in $listFiles; do
	    design=$(basename $files | sed -e 's/.txt//g' -e 's/\./\t/g')

	    total_genes=$(cat $files | sed '1d' | wc -l)

	    ## Get the number of significant genes based on their B-statisitcs, adjusted P value and average expression
	    cat $files | sed -e '1d' -e 's/ /./g' | awk -vb=$_B '{if($25>=b) print $0}' > $_BVAL
	    cat $files | sed -e '1d' -e 's/ /./g' | awk -vb=$_B '{if($25>=b) print $0}' | awk -vp=$_aP '{if($24<=p) print $0}' > $_ADJPVAL
	    cat $files | sed -e '1d' -e 's/ /./g' | awk -vb=$_B '{if($25>=b) print $0}' | awk -ve=$_aE '{if($21>=e) print $0}' > $_AVGEX
	    
	    
	    ## count the total number of diff expressed significant genes
	    if [ $i == "total" ]; then
		
		bval=$(cat $_BVAL | wc -l)
		adjpval=$(cat $_ADJPVAL | wc -l)
		avgex=$(cat $_AVGEX | wc -l)

		## get the minimum, maximum, average of B-statistic
		## example: B=3 ; value = ( exp(3) / ( exp(3) + 1 ) ) * 100
		## Value is the chance that one gene is differentially expressed
		for alpha in bval adjpval avgex ; do
		    if [ $alpha == "bval" ]; then
			output=$_BVAL
		    elif [ $alpha == "adjpval" ]; then
			output=$_ADJPVAL
		    elif [ $alpha == "avgex" ]; then
			output=$_AVGEX
		    fi

		    _min=$(cat ${output} | awk '{print $25}' | sort -n | sed '2,$d')
		    _max=$(cat ${output} | awk '{print $25}' | sort -nr | sed '2,$d')
		    _avg=$(cat ${output} | awk '{print $25}' | awk '{ sum += $1 } END { print sum/NR }')
		    _median=$(cat ${output} | awk '{print $25}' | sort -n | awk ' { a[i++]=$1; } END { print a[int(i/2)]; }')

		    eval ${alpha}_min=$(printf "%.2f%%" $(echo "(e($_min)/(e($_min)+1))*100" | bc -l))
		    eval ${alpha}_max=$(printf "%.2f%%" $(echo "(e($_max)/(e($_max)+1))*100" | bc -l))
		    eval ${alpha}_avg=$(printf "%.2f%%" $(echo "(e($_avg)/(e($_avg)+1))*100" | bc -l))
		    eval ${alpha}_median=$(printf "%.2f%%" $(echo "(e($_median)/(e($_median)+1))*100" | bc -l))
		done

		scaledSummary $_TMP $i $_B
		
		echo "Done counting the total number of differentially expressed genes for $files"

	    fi
	    
	    ## count the up regulated genes
	    if [ $i == "up" ]; then

		bval=$(cat $_BVAL | awk -vh=$_foldH '{if($20>=h) print $0}' | wc -l)
		adjpval=$(cat $_ADJPVAL | awk -vh=$_foldH '{if($20>=h) print $0}' | wc -l)
		avgex=$(cat $_AVGEX | awk -vh=$_foldH '{if($20>=h) print $0}' | wc -l)

		for alpha in bval adjpval avgex ; do
		    if [ $alpha == "bval" ]; then
			output=$_BVAL
		    elif [ $alpha == "adjpval" ]; then
			output=$_ADJPVAL
		    elif [ $alpha == "avgex" ]; then
			output=$_AVGEX
		    fi

		    _min=$(cat ${output} | awk -vh=$_foldH '{if($20>=h) print $0}' | awk '{print $25}' | sort -n | sed '2,$d')
		    _max=$(cat ${output} | awk -vh=$_foldH '{if($20>=h) print $0}' | awk '{print $25}' | sort -nr | sed '2,$d')
		    _avg=$(cat ${output} | awk -vh=$_foldH '{if($20>=h) print $0}' | awk '{print $25}' | awk '{ sum += $1 } END { print sum/NR }')
		    _median=$(cat ${output} | awk -vh=$_foldH '{if($20>=h) print $0}' | awk '{print $25}' | sort -n | awk ' { a[i++]=$1; } END { print a[int(i/2)]; }')

		    eval ${alpha}_min=$(printf "%.2f%%" $(echo "(e($_min)/(e($_min)+1))*100" | bc -l))
		    eval ${alpha}_max=$(printf "%.2f%%" $(echo "(e($_max)/(e($_max)+1))*100" | bc -l))
		    eval ${alpha}_avg=$(printf "%.2f%%" $(echo "(e($_avg)/(e($_avg)+1))*100" | bc -l))
		    eval ${alpha}_median=$(printf "%.2f%%" $(echo "(e($_median)/(e($_median)+1))*100" | bc -l))
		done

		scaledSummary $_TMP $i $_B

		echo "Done counting the total number of up-regulated genes for $files"
	    fi

	    ## count the down regulated genes
	    if [ $i == "down" ]; then

		bval=$(cat $_BVAL | awk -vl=$_foldL '{if($20<=l) print $0}' | wc -l)
		adjpval=$(cat $_ADJPVAL | awk -vl=$_foldL '{if($20<=l) print $0}' | wc -l)
		avgex=$(cat $_AVGEX | awk -vl=$_foldL '{if($20<=l) print $0}' | wc -l)

		for alpha in bval adjpval avgex ; do
		    if [ $alpha == "bval" ]; then
			output=$_BVAL
		    elif [ $alpha == "adjpval" ]; then
			output=$_ADJPVAL
		    elif [ $alpha == "avgex" ]; then
			output=$_AVGEX
		    fi

		    _min=$(cat ${output} | awk -vl=$_foldL '{if($20<=l) print $0}' | awk '{print $25}' | sort -n | sed '2,$d')
		    _max=$(cat ${output} | awk -vl=$_foldL '{if($20<=l) print $0}' | awk '{print $25}' | sort -nr | sed '2,$d')
		    _avg=$(cat ${output} | awk -vl=$_foldL '{if($20<=l) print $0}' | awk '{print $25}' | awk '{ sum += $1 } END { print sum/NR }')
		    _median=$(cat ${output} | awk -vl=$_foldL '{if($20<=l) print $0}' | awk '{print $25}' | sort -n | awk ' { a[i++]=$1; } END { print a[int(i/2)]; }')

		    eval ${alpha}_min=$(printf "%.2f%%" $(echo "(e($_min)/(e($_min)+1))*100" | bc -l))
		    eval ${alpha}_max=$(printf "%.2f%%" $(echo "(e($_max)/(e($_max)+1))*100" | bc -l))
		    eval ${alpha}_avg=$(printf "%.2f%%" $(echo "(e($_avg)/(e($_avg)+1))*100" | bc -l))
		    eval ${alpha}_median=$(printf "%.2f%%" $(echo "(e($_median)/(e($_median)+1))*100" | bc -l))
		done

		scaledSummary $_TMP $i $_B

		echo "Done counting the total number of down-regulated genes for $files"
	    fi

	done
    done
done

cleaning $_TMP $summary
