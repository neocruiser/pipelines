#!/bin/bash


## cytoscape network table
_CSV=$1
## annotation output from IPS
_TSV=$2


## choose between go, reactome, kegg, or gene expression
printf "Choose term to be extracted from annotation (go|reactome|kegg|deg) -> "
read _KEY
## name of the network in cytoscape
printf "Input name of the network ->"
read NET_NAME


## get the total number of nodes retrieved from cytoscape
## and get the total number of entries annotated through interpro scan/HMM
_NODES=$(cat $_CSV | sed '1d' | wc -l)
_ANNOTATED=$(cat $_TSV | wc -l)


## each if statement will clean the CSV output of cytoscape
## then will select the columns $14 for GO-terms
## or column $15 for kegg and reactome
## Next duplicate lines are merged
## Entries are combined also if the same transcript produced 2 or more annotations
## the end of each if statement will create a separate output file for each term search
## Gene expression retrival requires the output of abundant transcript analysis
## this analysis is done either with eXpress, Salmon or kallisto
## using different p-values, fold change and pairwise comparisons
if [ $_KEY == "go" ]
then
    _out=go.$NET_NAME.ids
    grep -Ff <(cat $_CSV | sed -e 's/","/\t/g' -e 's/"//g' -e '1d'| cut -f7) $_TSV | \
        cut -f1,14 | sed -e 's/ /./g' -e 's/..\t/\t/g' | sort - | uniq | grep "^TRINITY" - |
        awk 'BEGIN{str = ""}{if ( str != $1 ) {if ( NR != 1 ){printf("\n")} {str = $1;printf("%s\t%s",$1,$2)}} else if ( str == $1 ) {printf("%s;",$2)}}END{printf("\n")}' | \
            grep -i "go" | \
            sed -e 's/;//g' -e 's/\./ /g' > $_out

    _COUNT=$(cat $_out | wc -l)
    echo -e "Output saved in $_out"
    echo -e "\n$_COUNT GO-terms found among $_NODES nodes and $_ANNOTATED annotated transcripts"


elif [ $_KEY == "reactome" ]
then
    _out=reactome.$NET_NAME.ids
    grep -Ff <(cat $_CSV | sed -e 's/","/\t/g' -e 's/"//g' -e '1d'| cut -f7) $_TSV | \
        cut -f1,15 | sed -e 's/ /./g' -e 's/..\t/\t/g' | sort - | uniq | grep "^TRINITY" - |
        awk 'BEGIN{str = ""}{if ( str != $1 ) {if ( NR != 1 ){printf("\n")} {str = $1;printf("%s\t%s",$1,$2)}} else if ( str == $1 ) {printf("%s;",$2)}}END{printf("\n")}' | \
            grep -i "reactome" | \
            sed -e 's/;//g' -e 's/\./ /g' > $_out

    _COUNT=$(cat $_out | wc -l)
    echo -e "Output saved in $_out"
    echo -e "\n$_COUNT reactome keys found among $_NODES nodes and $_ANNOTATED annotated transcripts"

elif [ $_KEY == "kegg" ]
then
    _out=kegg.$NET_NAME.ids
    grep -Ff <(cat $_CSV | sed -e 's/","/\t/g' -e 's/"//g' -e '1d'| cut -f7) $_TSV | \
        cut -f1,15 | sed -e 's/ /./g' -e 's/..\t/\t/g' | sort - | uniq | grep "^TRINITY" - |
        awk 'BEGIN{str = ""}{if ( str != $1 ) {if ( NR != 1 ){printf("\n")} {str = $1;printf("%s\t%s",$1,$2)}} else if ( str == $1 ) {printf("%s;",$2)}}END{printf("\n")}' | \
            grep -i "kegg" | \
            sed -e 's/;//g' -e 's/\./ /g' -e's/|Meta.*$//g' > $_out

    _COUNT=$(cat $_out | wc -l)
    echo -e "Output saved in $_out"
    echo -e "\n$_COUNT KEGG found among $_NODES nodes and $_ANNOTATED annotated transcripts"

elif [ $_KEY == "deg" ]
then
    RAW_DEG=trans_counts.counts.matrix.conditionA_vs_conditionB.DESeq2.DE_results
    _out=deg.$NET_NAME.ids

    if [ ! -f "$RAW_DEG" ]
    then
        printf "Input p-value score of transcript abundance [e-1..10] -> "
        read _PVAL
        printf "Input the fold change of transcript log ratios [1-4] -> "
        read _CFOLD
        echo -e "\nRetrieving transcript abundance results from Bridges ..."
        rsync -avP bassim@bridges.psc.xsede.org:/pylon2/oc4ifip/bassim/ganglia/trinity/trinity_out_dir_00000/deg.eXpress.111718/DESeq2.eXpress.tissue.p$_PVAL.c$_CFOLD.111718/trans_counts.counts.matrix.conditionA_vs_conditionB.DESeq2.DE_results .
    fi

    cat $RAW_DEG | sed '2, $d' > $_out
    grep -Fwf <(cat $_CSV| sed -e 's/","/\t/g' -e 's/"//g' -e '1d'| cut -f7) $RAW_DEG >> $_out

    _COUNT=$(cat $_out | sed '1d' | wc -l)
    echo -e "Output saved in $_out"
    echo -e "\n$_COUNT transcript counts and expresion p-values retrieved for $_NODES nodes"


fi
