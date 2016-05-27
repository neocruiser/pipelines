#! /usr/bin/bash

genome[1]=diatome/Fracy1_goinfo_FilteredModels1.tab
genome[2]=diatome/Psemu1_GeneCatalog_proteins_20111011_GO.tab
genome[3]=diatome/Thaps3_chromosomes_goinfo_FilteredModels2.tab
genome[4]=diatome/Phatr2_chromosomes_goinfo_FilteredModels2.tab
genome[5]=laby/Aplke1_GeneCatalog_proteins_20121220_GO.tab
genome[6]=laby/Aurli1_GeneCatalog_proteins_20120618_GO.tab
genome[7]=laby/Schag1_GeneCatalog_proteins_20121220_GO.tab
genome[8]=oomycete/Phyca11_filtered_proteins_GO.tab
genome[9]=oomycete/Phyci1_GeneCatalog_proteins_20120612_GO.tab
genome[10]=oomycete/Physo3_GeneCatalog_proteins_20110401_GO.tab
genome[11]=oomycete/pramor1_GeneCatalog_proteins_2011_GO.tab

## choose GO ID after enrichment analysis
go="GO:0005737"
cc="cytoplasm"
## Choose either biological_process or molecular_function
term="molecular_function"
te="MF"
## minimum number of genes per GO term
Z=1
## output file
output=cc2cytoplasm_gene1_mf
touch $output

for i in 1 2 3 4 5 6 7 8 9 10 11
do
    genome=${genome[${i}]}
## get the 5 first letters of each filename as genome names
    name=$(echo $genome | perl -pe 's/^.*\/(.....).*[0-9]\_.*$/\1/g')
## get the number of genes per GO term
    cat $genome | \
        grep $go | \
        awk '{print $1}' | \
        grep -Fwf - $genome | \
        grep $term | \
        cut -f 3 | \
        sort - | \
        uniq -c | \
        awk -vZ="$Z" '{if ($1 >= Z) print $0}' | \
        sort -n | \
        perl -pe 's/([0-9]) /\1\t/g' | \
        awk -vgo="$go" -vcc="$cc" -vte="$te" -vname="$name" \
        '{print $0,"\t",te,"\t",go,"\t",cc,"\t",name,"\t","1"}' \
        >> $output
done

