#!/bin/bash

out=summary.txt
final=summary.txt
rm $final

for d in sample.{1..13}.*; do
    for f in $d/kraken.mpa.report; do
	abbrv=$(echo $d | cut -f2 -d '.')
	sed -e 's/$/\t'$abbrv'/g' <(cat $f) >> $out
    done
done


## Run in R
## combine the results of all samples by phylum
#require(tidyr)
#require(dplyr)
#dat <- read.table("./summary.txt", header = F) %>%
#    spread(V1,V2)
#write.table(t(dat),"summary.txt",quote=F,sep="\t")


