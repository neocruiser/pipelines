#!/bin/bash

# 1. Translate nucleotide to protein sequences
# 2. Run an interpro scan on Protein queries
# 3. Use IPS unfiltered output with this script

#####################################################################################################################
################################################################################################# DECLARING FUNCTIONS
#####################################################################################################################
function quest () {
    local VAR1=$1
    echo -e "\n---------------------------------------------\n"
    echo "a. Number of proteins found in PANTHER"
    echo "b. Species names"
    echo "c. Most abundant proteins"
    echo "d. Output all alignment info"
    echo "e. Get contig IDs and IPS/PTHR functions"
    echo "f. Get contig IDs and IPS/PTHR/NR functions"
    echo "w. debugging: Extract fasta file"
    echo
    echo "===DANGER==="
    echo "x. Delete the temporary file"
    echo "y. Delete ALL temporary files"
    echo "===EXIT==="
    echo "z. Do nothing and exit"
    echo
    printf "Input only one of the above letters -> "
    read CHOICE
}

## get the location of tsv files (user defined)
function files () {
    local OUTPUT=$1
    local EXE=$2
    if [ "$EXE" == txt ]; then
        echo
	      find $OUTPUT -maxdepth 2 -iname "*blast*$EXE" | nl | tee $EXE.tmp
    elif [ "$EXE" == fa ]; then
        echo
        find $(dirname $(dirname $FILENAME))  -maxdepth 2 -iname "*trinity*fa*" | nl | tee $EXE.tmp
    else
        echo
	      find $OUTPUT -maxdepth 2 -iname "*$EXE" | nl | tee $EXE.tmp
    fi
    echo
}


function summary () {
## alignment length
    local VAR1=$1
## evalue
    local VAR2=$2
## Filename
    local VAR3=$3
    local __out=$VAR3.$VAR1.10-$VAR2.tmp
## create a correct e-value number by repeating zeros
    ZEROS=$(seq -s. "$(echo "${VAR2}+1" | bc)" | tr -d '[:digit:]' | sed 's/./0/g')
    local VARe="0.${ZEROS}1"
## path to panther database
    if [ "$4" == "bridges" ]; then
        VAR4=/pylon1/oc4ifip/bassim/db/panther
    elif [ "$4" == "lired" ]; then
        VAR4=/gpfs/scratch/ballam/db/panther2
    else
        VAR4="$4"
    fi
 ## select panther only proteins by evalue and alignment length
## the awk part will merge PANTHER identified genes and IPS output results using second column of output 1 and third column of output 2
# Panther columns $2=ID $3=family_name $4=subfamily_name $5=GOs
    awk 'NR==FNR {h[$2] = sprintf ("%s\t%s\t%s\t",$1,$2,$3); next} {print h[$3],$1,$3,$4,$5,$6}' <(grep -RFwf <(cat $VAR3 | cut -f1,4,5,7,8,9 | awk -va="$VAR1" -vp="$VARe" '{n=$4-$5?$5-$4:$4-$5; if(n>=a && $6<=p && $2 == "PANTHER") print $0}' | cut -f3  | sed '/^.*\:.*$/d') $VAR4 | sed 's/ /./g' | sort - | cut -f1,3,4,5 | uniq | sed 's/\:SF.*\t/\t/g') <(cat $VAR3 | cut -f1,4,5,7,8,9 | awk -va="$VAR1" -vp="$VARe" '{n=$4-$5?$5-$4:$4-$5; if(n>=a && $6<=p && $2 == "PANTHER") print $0}' | sort -k3 - | sed '/^.*\:.*$/d' | uniq) | sort -k1 - > $__out


}


function extra () {
    local VAR1=$1
    local VAR2=$2
    if [ "$VAR2" == 1 ]; then
## count number of hits
	cat $VAR1 | wc -l
    elif [ "$VAR2" == 2 ]; then
## select species column
	cat $VAR1 | sort -k2 - | cut -f1 | egrep -o "0_.*:" | sed -e 's/0_//g' -e 's/://g' | sort - | uniq -c | sort -nrk1 > $FILENAME.PANspecies.LEN$ALIGNMENT.EVAL$EVAL.tab
    elif [ "$VAR2" == 3 ]; then
## select protein function column
	cat $VAR1 | cut -f3 | sort - | uniq -c | sort -nr | sed 's/\./ /g' > $FILENAME.PANfunctions.LEN$ALIGNMENT.EVAL$EVAL.tab
    elif [ "$VAR2" == 4 ]; then
## get all entries
	cp $VAR1 $FILENAME.PANprots.LEN$ALIGNMENT.EVAL$EVAL.tab
	fi
}


function guidelines () {
    tmp=$FILENAME.$ALIGNMENT.10-$EVAL.tmp

    if [ "$CHOICE" == a ]; then
	echo -e "\nHere is the number of hits found in 13GB of a PANTHER database:"
	extra $tmp 1
	echo

    elif [ "$CHOICE" == b ]; then
	      extra $tmp 2
        _LN=$(grep -c "^" $FILENAME.PANspecies.LEN$ALIGNMENT.EVAL$EVAL.tab)
        _PN=$(head -n 10 $FILENAME.PANspecies.LEN$ALIGNMENT.EVAL$EVAL.tab | grep -c "^")
	echo -e "\nHere is $_PN of the $_LN most abundant species found by aligning your proteins to 13GB of PANTHER sequences:\nwait ..."
	cat $FILENAME.PANspecies.LEN$ALIGNMENT.EVAL$EVAL.tab | column -t | less
	head -n 10 $FILENAME.PANspecies.LEN$ALIGNMENT.EVAL$EVAL.tab
	echo -e "\nI also put the distribution of all other species in:" $FILENAME.PANspecies.LEN$ALIGNMENT.EVAL$EVAL.tab
	echo

    elif [ "$CHOICE" == c ]; then
	      extra $tmp 3
        _LN=$(grep -c "^" $FILENAME.PANfunctions.LEN$ALIGNMENT.EVAL$EVAL.tab)
        _PN=$(head -n 10 $FILENAME.PANfunctions.LEN$ALIGNMENT.EVAL$EVAL.tab | grep -c "^")
	echo -e "\nHere is $_PN of the $_LN most abundant protein functions found among the queried sequences in PANTHER:\nwait ..."
	cat $FILENAME.PANfunctions.LEN$ALIGNMENT.EVAL$EVAL.tab | sed 's/\./ /g' | less
	head -n 10 $FILENAME.PANfunctions.LEN$ALIGNMENT.EVAL$EVAL.tab | sed 's/\./ /g'
	echo -e "\nI also put the full list in:" $FILENAME.PANfunctions.LEN$ALIGNMENT.EVAL$EVAL.tab
	echo

    elif [ "$CHOICE" == d ]; then
	echo -e "\nWriting the protein alignment output from PANTHER to a file. It will contains E-values, protein acc. numbers, contig IDs, and GOs"
	extra $tmp 4
	echo -e "\nI also put the description of all entries found in:" $FILENAME.PANprots.LEN$ALIGNMENT.EVAL$EVAL.tab
	echo -e "\nDone"
	echo

    elif [ "$CHOICE" == e ]; then
# Extract annotated genes from IPS tsv output
        echo -e "\nStep 1: Getting descriptions from InterPro scans ..."
## create a correct e-value number by repeating zeros
    ZEROS=$(seq -s. "$(echo "${EVAL}+1" | bc)" | tr -d '[:digit:]' | sed 's/./0/g')
    local VARe="0.${ZEROS}1"
    cat $FILENAME | sed 's/ /./g' | cut -f1,9,13 | awk -ve="$VARe" '{if($2<=e)print$1,$3}' | sed 's/_. / /g' | sort - | uniq > $FILENAME.id2description.LEN$ALIGNMENT.EVAL$EVAL.tmp

# then combine with Panther annotation: $4 is contig ID and $3 is Panther description
        echo "Step 2: Getting descriptions from PANTHER database ..."
paste <(awk '{print $4}' $tmp ) <(awk '{print $3}' $tmp ) | sed 's/_.\t/\t/g' | sort -k1 - | uniq >> $FILENAME.id2description.LEN$ALIGNMENT.EVAL$EVAL.tmp

        echo "Step 3: Merging descritpions and removing duplicates ..."
# Files contain contig IDs and protein description
## the source of the awk part: http://stackoverflow.com/questions/17832631/combine-rows-in-linux
cat $FILENAME.id2description.LEN$ALIGNMENT.EVAL$EVAL.tmp | sort - | grep "^TRINITY" - | awk 'BEGIN{str = ""}{if ( str != $1 ) {if ( NR != 1 ){printf("\n")} {str = $1;printf("%s\t%s",$1,$2)}} else if ( str == $1 ) {printf("%s;",$2)}}END{printf("\n")}' > $FILENAME.id2description.LEN$ALIGNMENT.EVAL$EVAL.txt
#rm $FILENAME.id2description.LEN$ALIGNMENT.EVAL$EVAL.tmp

_FULL=$(cat $FILENAME.id2description.LEN$ALIGNMENT.EVAL$EVAL.txt | cut -f2 | sed '/^\s*$/d' | wc -l)
_ALL=$(cat $FILENAME.id2description.LEN$ALIGNMENT.EVAL$EVAL.txt | wc -l)
        echo "There is $_FULL annotated proteins found in all databases among $_ALL aligned contigs"
        echo

    elif [ "$CHOICE" == f ]; then
# Check for NCBI database files; download them if they dont exist
        _ftp="ftp://ftp.ncbi.nlm.nih.gov/gene/DATA"

        if [ "$DB" == "bridges" ]; then
            gene_info=/pylon1/oc4ifip/bassim/db/ncbi/gene_info
            gene2accession=/pylon1/oc4ifip/bassim/db/ncbi/gene2accession
            else
            gene_info=/gpfs/scratch/ballam/db/ncbi/gene_info
            gene2accession=/gpfs/scratch/ballam/db/ncbi/gene2accession
        fi

        if [ -f "$gene_info" ]; then
            echo -e "\nStep 1: Checking for gene information from NCBI ... OK"
        else
            echo "!!!ERROR!!! File not found. Downloading gene_info file. Wait ..."
            if [ "$DB" == "bridges" ]; then
                $_path=/pylon1/oc4ifip/bassim/db/ncbi
                mkdir -p $_path
                wget -O $_path/gene_info.gz $_ftp/gene_info.gz
                gunzip -d $_path/gene_info.gz
            elif [ "$DB" == "lired" ]; then
                $_path=/gpfs/scratch/ballam/db/ncbi
                mkdir -p $_path
                wget -O $_path/gene_info.gz $_ftp/gene_info.gz
                gunzip -d $_path/gene_info.gz
            fi
            echo "NCBI gene information has been downloaded"
        fi

        if [ -f "$gene2accession" ]; then
            echo "Step 2: Checking for protein accessions from NCBI ... OK"
        else
            echo "!!!ERROR!!! File not found. Downloading gene2accession file. Wait ..."
            if [ "$DB" == "bridges" ]; then
                $_path=/pylon1/oc4ifip/bassim/db/ncbi
                mkdir -p $_path
                wget -O $_path/gene2accession.gz $_ftp/gene2accession.gz
                gunzip -d $_path/gene2accession.gz
            elif [ "$DB" == "lired" ]; then
                $_path=/gpfs/scratch/ballam/db/ncbi
                mkdir -p $_path
                wget -O $_path/gene2accession.gz $_ftp/gene2accession.gz
                gunzip -d $_path/gene2accession.gz
            fi
            echo "NCBI protein accession numbers have been downloaded"
        fi

## Append NCBIs gene information, gene ID, protein accession nb, and taxid to each contig
            files $FILE__PATH txt
            printf "a. Choose one DIAMOND txt file from the list above (number) -> "
            read _FD
            _diamond=$(awk -vf="$_FD" '{if ($1 == f) print $2}' txt.tmp)
            rm txt.tmp
            printf "b. Choose an E-value for an alignment score [e-0..35] -> "
            read ED
            ZD=$(seq -s. "$(echo "${ED}+1" | bc)" | tr -d '[:digit:]' | sed 's/./0/g')
            local _ev="0.${ZD}1"

            _diamond_summary=NR-PTHR-IPS.LEN.tmp

            echo -e "\nStep 3: Getting descriptions from BLAST (diamond/nr output). Wait ..."
            awk 'NR==FNR {b[$7] = sprintf ("%s\t%s\t%s\t%s\t%s\t%s\t%s\t",$1,$2,$3,$4,$5,$6,$7); next} {print b[$2],$1,$2}' <(awk 'NR==FNR {a[$2] = sprintf ("%s\t%s\t%s\t%s\t%s\t",$1,$2,$3,$4,$5); next} {print a[$1],$1,$2}' <(cat $gene_info | grep -Fwf <(cat $gene2accession | cut -f2,6 | grep -Fwf <(cat $_diamond  | cut -f1,2,11 | awk -ve="$_ev" '{if($3<=e)print$0}' | cut -f1,2 | sed 's/|/\t/g' | cut -f5) - | cut -f1) - | cut -f1,2,3,9,10 | sed 's/ /./g') <(cat $gene2accession | cut -f2,6 | grep -Fwf <(cat $_diamond  | cut -f1,2,11 | awk -ve="$_ev" '{if($3<=e)print$0}' | cut -f1,2 | sed 's/|/\t/g' | cut -f5) -)) <(cat $_diamond  | cut -f1,2,11 | awk -ve="$_ev" '{if($3<=e)print$0}' | cut -f1,2 | sed 's/|/\t/g' | cut -f1,5) | awk '{if ($2 == $6 && $7 == $9) print $1,$2,$3,$4,$5,$7,$8}' | sort - | uniq | cut -f1,4,7 -d ' ' | sort -k2 | uniq  > $_diamond_summary

# restructure the summary
            paste <(cut -f3 -d' ' $_diamond_summary) <(cut -f1,2 -d' ' $_diamond_summary | sed 's/ /:/g') > $FILENAME.id2description.NR-PTHR-IPS.LEN$ALIGNMENT.EVAL$EVAL.tmp
#            rm $_diamond_summary

# merge summary with IPS and Panther summaries
# Extract annotated genes from IPS tsv output
        echo "Step 4: Getting descriptions from InterPro scans ..."
## create a correct e-value number by repeating zeros
    ZEROS=$(seq -s. "$(echo "${EVAL}+1" | bc)" | tr -d '[:digit:]' | sed 's/./0/g')
    local VARe="0.${ZEROS}1"
    cat $FILENAME | sed 's/ /./g' | cut -f1,9,13 | awk -ve="$VARe" '{if($2<=e)print$1,$3}' | sed 's/_. / /g' | sort - | uniq >> $FILENAME.id2description.NR-PTHR-IPS.LEN$ALIGNMENT.EVAL$EVAL.tmp

# then combine with Panther annotation: $4 is contig ID and $3 is Panther description
        echo "Step 5: Getting descriptions from PANTHER database ..."
paste <(awk '{print $4}' $tmp ) <(awk '{print $3}' $tmp ) | sed 's/_.\t/\t/g' | sort -k1 - | uniq >> $FILENAME.id2description.NR-PTHR-IPS.LEN$ALIGNMENT.EVAL$EVAL.tmp

        echo "Step 6: Merging descritpions and removing duplicates ..."
# Files contain contig IDs and protein description
cat $FILENAME.id2description.NR-PTHR-IPS.LEN$ALIGNMENT.EVAL$EVAL.tmp | sort - | grep "^TRINITY" - | awk 'BEGIN{str = ""}{if ( str != $1 ) {if ( NR != 1 ){printf("\n")} {str = $1;printf("%s\t%s",$1,$2)}} else if ( str == $1 ) {printf("%s;",$2)}}END{printf("\n")}' > $FILENAME.id2description.NR-PTHR-IPS.LEN$ALIGNMENT.EVAL$EVAL.txt
#rm $FILENAME.id2description.LEN$ALIGNMENT.EVAL$EVAL.tmp

_FULL=$(cat $FILENAME.id2description.NR-PTHR-IPS.LEN$ALIGNMENT.EVAL$EVAL.txt | cut -f2 | sed '/^\s*$/d' | wc -l)
_ALL=$(cat $FILENAME.id2description.NR-PTHR-IPS.LEN$ALIGNMENT.EVAL$EVAL.txt | wc -l)
        echo "There is $_FULL annotated proteins found in all databases among $_ALL aligned contigs"
        echo

    elif [ "$CHOICE" == x ]; then
	rm $tmp
	echo "$tmp was deleted!"
	exit
    elif [ "$CHOICE" == y ]; then
	find $(dirname $FILENAME) -maxdepth 1 -iname "*tmp" -exec rm {} \;
	echo "All temperory files were deleted!"
	exit
    elif [ "$CHOICE" == z ]; then
	echo -e "Bye!"
	exit
    else
	echo "Stop messing around!"
	echo
    fi
}


## TESTING (ongoing preogress)
function extractFasta () {
##panther output
    if [ -f $FILENAME.PANprots.LEN$ALIGNMENT.EVAL$EVAL.tab ]; then
	printf "String to search for -> "
	read STRING
	egrep -i "$STRING" $FILENAME.PANprots.LEN$ALIGNMENT.EVAL$EVAL.tab >> $FILENAME.PANselected.tsv
	cat $FILENAME.PANselected.tsv | cut -f1 | sed 's/..$//g' | sort - | uniq >> $FILENAME.PANcontigs

	else
	extra $TMP 4
    fi

##all db##
# interpro output
col1="$7-$6"
col2="$8-$7"
for i in $col1 $col2; do
    cat $FILENAME | sed 's/ /./g' | awk '{if($8<=0.0000000001 || $9<=0.0000000001)print $0}' | awk -vc="$i" '{n=c; if(n>=50) print$0}' >> $FILENAME.ALLselected.tsv
done
cat $FILENAME.ALLselected.tsv | cut -f1 | sed 's/..$//g' | sort - | uniq >> $FILENAME.Allcontigs

#cat alp | sort - | uniq | alpp
#rm alp

##compare##
#comm -13 <(cat alpp | sort -| uniq) <(cat pa | sort -| uniq) > pas
#rm pa
#cat pas >> alpp
#rm pas


}


#####################################################################################################################
########################################################################## RUN ANALYSES FOR BLAST AND PANTHER OUTPUTS
#####################################################################################################################
echo
echo -e "\n(p) analyses on PANTHER output, diamond included\n(b) analyses on BLAST output, STRING included\n(s) analyses on Interpro scans"
printf "Choose between a panther, blast or summary analysis (p|b|s) -> "
read ANALYSIS
echo

## user input when executing this script
FILE__PATH=$1

if [ "$ANALYSIS" == s ]; then
#===================== summary analysis, showcase all database entries
#=====================================================================
    while true; do
        files $FILE__PATH tsv
        printf "1. Choose one annotated file from the list above (number) -> "
        read _FN
        FILENAME=$(awk -vf="$_FN" '{if ($1 == f) print $2}' ips.tmp)
        rm ips.tmp
	      printf "2. Choose an E-value for an alignment score [e-0..35] -> "
	      read NUM
	      ZEROS=$(seq -s. "$(echo "$NUM+1" | bc)" | tr -d '[:digit:]' | sed 's/./0/g')
	      EVAL="0.${ZEROS}1"
	      echo "Below is the number of proteins annotated for each database at an E-value of $EVAL"
	      cat $FILENAME | sed 's/ /./g' | cut -f4,9 | awk -ve="$EVAL" '{if($2<=e)print$1}' | sort - | uniq -c | sort -n
        read -n 1 -s
    done

elif [ "$ANALYSIS" == p ]; then
#==================================== analysis on panther output files
#=====================================================================
    files $FILE__PATH tsv
    printf "1. Choose one PANTHER file from the list above (number) -> "
    read _FN
    FILENAME=$(awk -vf="$_FN" '{if ($1 == f) print $2}' tsv.tmp)
    rm tsv.tmp
    printf "2. Choose an acceptable alignment length [20..1000] -> "
    read ALIGNMENT
    printf "3. Choose an E-value for an alignment score [e-0..35] -> "
    read EVAL
    printf "4. Choose a PANTHER database (bridges|lired|...) -> "
    read DB

    COUNTS=$(grep -c "^" $FILENAME)
    TMP=$FILENAME.$ALIGNMENT.10-$EVAL.tmp

    if [ -f "$TMP" ]; then
	while [ -f "$TMP" ]; do
	    quest $COUNTS
	    guidelines
	    read -n 1 -s
	done
    else
	echo -e "\nSearching in 13GB of PANTHER genomes. I found $COUNTS hits from interpro scans ..."
	summary $ALIGNMENT $EVAL $FILENAME $DB
	while [ -f "$TMP" ]; do
	    quest $COUNTS
	    guidelines
	    read -n 1 -s
	done
    fi

elif [ "$ANALYSIS" == b ]; then
#======================================================== BLAST output
#=====================================================================
    while true; do
    files $FILE__PATH txt
    printf "1. Choose one BLAST file from the list above (number) -> "
    read _FN
    FILENAME=$(awk -vf="$_FN" '{if ($1 == f) print $2}' blast.tmp)
    printf "2. Choose an E-value for an alignment score [e-0..35] -> "
    read REP_
## create a correct e-value number by repeating zeros
    ZEROS=$(seq -s. "$(echo "${REP_}+1" | bc)" | tr -d '[:digit:]' | sed 's/./0/g')
    EVAL="0.${ZEROS}1"
    rm blast.tmp
## get the synthax used to identify each assembled contig


    echo -e "\n--------------------------------------------------------------"
    echo "a. Count sequences based on length, eval, identity, mismatches..."
    echo "b. Extract fasta sequences from QUALITY blast"
    echo "c. Extract fasta sequences that do NOT match the quality blast alignment"
    echo "d. Associate STRING connections to contig IDs"
    echo "e. Do nothing and exit"
    echo
    printf "Input only one of the above letters -> "
    read CHOICE_BLAST

    if [ "$CHOICE_BLAST" == a ]; then
## Get distribution of number of hits from ANY blast output
	for i in 2 4 10 11 12 14 15; do
## blast categories found in any blast+ output
	    s[2]="1. Query length [max-min] ->"
	    s[4]="2. Target length [max-min] ->"
	    s[10]="3. Bit Score [max-min] (Higher better) ->"
	    s[11]="4. Alignment length [max-min] ->"
	    s[12]="5. Percentage identity [max-min] ->"
	    s[14]="6. Indentical nt [max-min] ->"
	    s[15]="7. Mismatches [max-min] (Important) ->"
## show both [max and min] intervals for each blast category
	    egrep "^[^#]" $FILENAME | awk -ve="$EVAL" -vg="$i" '{if($9<=e)print$g}' | sort - | uniq -c | sort -k2 -nr | awk -vf="${s[$i]}" 'NR==1{print f,$2"--"}END{print$2}' | sed 'N;s/\n//g'
## get the maximum range of each category
	    _MAX=$(egrep "^[^#]" $FILENAME | awk -ve="$EVAL" -vg="$i" '{if($9<=e)print$g}' | sort - | uniq -c | sort -k2 -nr | awk 'NR==1{print$2}')
## show the number of hits for the different categories
## randomly divide the max number of hits
	    h3=$(echo "$_MAX / ($_MAX/50)" | bc)
	    for f in 2 10 20 $h3; do
		_DIV=$(echo "$_MAX / $f" | bc | awk '{printf "%.0f\n",$1}')
		if [ ! "$i" == 15 ]; then
		   egrep "^[^#]" $FILENAME | awk -ve="$EVAL" -vg="$i" '{if($9<=e)print$g}' | sort - | uniq -c | sort -k2 -nr | awk -vm="$_DIV" '{if($2>=m) print $1}' | awk -vx="$_DIV" '{sum+=$1;print "Nb of hits >= "x,"is "sum}' | tail -n1
		else
		    egrep "^[^#]" $FILENAME | awk -ve="$EVAL" -vg="$i" '{if($9<=e)print$g}' | sort - | uniq -c | sort -k2 -nr | awk -vm="$_DIV" '{if($2<=m) print $1}' | awk -vx="$_DIV" '{sum+=$1;print "Nb of hits <= "x,"is "sum}' | tail -n1
		fi
	    done
	done

	while true; do
## Filter number of hits using the blast categories
	    printf "Select scores from the 7 above blast categories (space separated) -> "
	    IFS=" "
	    read _Q _T _B _A _P _I _M
	    HITS=$(egrep "^[^#]" $FILENAME | awk -ve="$EVAL" -vq="$_Q" -vt="$_T" -vb="$_B" -va="$_A" -vp="$_P" -vi="$_I" -vm="$_M" '{if($9<=e && $2>=q && $4>=t && $10>=b && $11>=a && $12>=p && $14>=i && $15<=m) print $1}' | sort - | uniq | wc -l)
	    echo "Number of unique hits found: $HITS"
	done

	elif [ "$CHOICE_BLAST" == b ]; then
## extract fasta based on quality sequences
    files $FILE__PATH fa
    printf "1. Choose TRANSCRIPTOME (number) -> "
    read _FN
    TRANSCRIPTOME=$(awk -vf="$_FN" '{if ($1 == f) print $2}' fasta.tmp)
    rm fasta.tmp
	  echo -e "\nQuery length, target length, bit score, align length, % identity, identical, mismatches"
	  printf "2. Input 7 scores for sequence quality filtering (space separated) -> "
	  IFS=" "
	  read _Q _T _B _A _P _I _M
	  echo "Wait ..."
	  _FA=$(dirname $TRANSCRIPTOME)/matching.BLAST.q$_Q.t$_T.b$_B.a$_A.p$_P.i$_I.m$_M.EVAL${REP_}.fa
	  cat $TRANSCRIPTOME | sed 's/.len.*$//g' | perl -ne 'if(/^>(\S+)/){$c=$i{$1}}$c?print:chomp;$i{$_}=1 if @ARGV' <(egrep "^[^#]" $FILENAME | awk -ve="$EVAL" -vq="$_Q" -vt="$_T" -vb="$_B" -va="$_A" -vp="$_P" -vi="$_I" -vm="$_M" '{if($9<=e && $2>=q && $4>=t && $10>=b && $11>=a && $12>=p && $14>=i && $15<=m) print $1}' | sort - | uniq) - > $_FA
	  echo "Extracted $(grep -c "^>" $_FA) from $(grep -c "^>" $TRANSCRIPTOME) fasta sequences based on the options you provided at an E-value of 10-$REP_"
	  read -n 1 -s


    elif [ "$CHOICE_BLAST" == c ]; then
## extract fasta that did not align by blast
    files $FILE__PATH fa
    printf "1. Choose TRANSCRIPTOME (number) -> "
    read _FN
    TRANSCRIPTOME=$(awk -vf="$_FN" '{if ($1 == f) print $2}' fasta.tmp)
    rm fasta.tmp
	  echo -e "\nQuery length, target length, bit score, align length, % identity, identical, mismatches"
	  printf "2. Input 7 scores for sequence quality filtering (space separated) -> "
	  IFS=" "
	  read _Q _T _B _A _P _I _M
	  echo "Wait ..."

	  _FA=$(dirname $TRANSCRIPTOME)/nonmatching.BLAST.q$_Q.t$_T.b$_B.a$_A.p$_P.i$_I.m$_M.EVAL${REP_}.fa
    cat $TRANSCRIPTOME  | sed 's/.len.*$//g' | perl -ne 'if(/^>(\S+)/){$c=$i{$1}}$c?print:chomp;$i{$_}=1 if @ARGV' <(comm -13 <(grep "^>" <(cat $TRANSCRIPTOME | sed 's/.len.*$//g' | perl -ne 'if(/^>(\S+)/){$c=$i{$1}}$c?print:chomp;$i{$_}=1 if @ARGV' <(egrep "^[^#]" $FILENAME | awk -ve="$EVAL" -vq="$_Q" -vt="$_T" -vb="$_B" -va="$_A" -vp="$_P" -vi="$_I" -vm="$_M" '{if($9<=e && $2>=q && $4>=t && $10>=b && $11>=a && $12>=p && $14>=i && $15<=m) print $1}' | sort - | uniq) -) | sort -) <(grep "^>" $TRANSCRIPTOME | cut -f1 -d ' ' | sort -) | sed 's/>//g') - > $_FA
	  echo "Extracted $(grep -c "^>" $_FA) from $(grep -c "^>" $TRANSCRIPTOME) fasta sequences based on the options you provided at an E-value of 10-$REP_"
	  read -n 1 -s


	elif [ "$CHOICE_BLAST" == d ]; then
## extract STRING connections
	  echo -e "\nQuery length, target length, bit score, align length, % identity, identical, mismatches"
	  printf "1. Input 7 scores for sequence quality filtering (space separated) -> "
	  IFS=" "
	  read _Q _T _B _A _P _I _M
    printf "2. Choose a STRING database (bridges|lired|...) -> "
    read DB

	  _TX=$(dirname $FILENAME)/string2string.q$_Q.t$_T.b$_B.a$_A.p$_P.i$_I.m$_M.EVAL${REP_}

## set the string DB path
    if [ "$DB" == "bridges" ]; then
        STRING_DB=/pylon1/oc4ifip/bassim/db/string/protein.links.full.v10.txt
    elif [ "$DB" == "lired" ]; then
        STRING_DB=/gpfs/scratch/ballam/db/string/protein.links.full.v10.txt
    else
        STRING_DB="$DB"
    fi

## get the number of hits found by blast in string database
    _strings=$(egrep "^[^#]" $FILENAME | awk -ve="$EVAL" -vq="$_Q" -vt="$_T" -vb="$_B" -va="$_A" -vp="$_P" -vi="$_I" -vm="$_M" '{if($9<=e && $2>=q && $4>=t && $10>=b && $11>=a && $12>=p && $14>=i && $15<=m) print $3}' | sort - | uniq | wc -l)
    _contigs=$(egrep "^[^#]" $FILENAME | awk -ve="$EVAL" -vq="$_Q" -vt="$_T" -vb="$_B" -va="$_A" -vp="$_P" -vi="$_I" -vm="$_M" '{if($9<=e && $2>=q && $4>=t && $10>=b && $11>=a && $12>=p && $14>=i && $15<=m) print $1}' | sort - | uniq | wc -l)

    echo -e "\n10-20 minutes to retrieve STRING-STRING connections for $_strings STRING IDs that matched to $_contigs genes\nWait..."

## get connections
	  grep -Fwf <(egrep "^[^#]" $FILENAME | awk -ve="$EVAL" -vq="$_Q" -vt="$_T" -vb="$_B" -va="$_A" -vp="$_P" -vi="$_I" -vm="$_M" '{if($9<=e && $2>=q && $4>=t && $10>=b && $11>=a && $12>=p && $14>=i && $15<=m) print $3}' | sort - | uniq) $STRING_DB > $_TX.tmp

## merge string/blast results (input 1) with string connections (input 2)
    awk 'NR==FNR {h[$3] = sprintf ("%s\t%s\t%s\t",$1,$3,$9); next} {print h[$1],$0}' <(egrep "^[^#]" $FILENAME | awk -ve="$EVAL" -vq="$_Q" -vt="$_T" -vb="$_B" -va="$_A" -vp="$_P" -vi="$_I" -vm="$_M" '{if($9<=e && $2>=q && $4>=t && $10>=b && $11>=a && $12>=p && $14>=i && $15<=m) print $0}') <(cat $_TX.tmp) | grep "^TRINITY" > $_TX.txt


	  echo "Extracted $(cat $_TX.txt | wc -l) STRING connections for $_contigs unique genes at an E-value of 10-$REP_"
    rm $_TX.tmp
	  read -n 1 -s


    else
	      exit
    fi
    done
## from blast of STRING
#grep -i "^trinity" selected.ips.allDB.pval.10-10.fa.string.blastx.67494.txt | awk '{if($9<=0.00001 && $10>=50 && $11>=100 && $12>=30 && $14>=50 && $15<=50) print $3}' | sort - | uniq | wc -l

fi



#cat $transcriptome  | sed 's/.len.*$//g' | perl -ne 'if(/^>(\S+)/){$c=$i{$1}}$c?print:chomp;$i{$_}=1 if @ARGV' <(comm -13 <(grep "^>" $selected | sort -) <(grep "^>" $transcriptome | cut -f1 -d ' ' | sort -) | sed 's/>//g') - > $output
