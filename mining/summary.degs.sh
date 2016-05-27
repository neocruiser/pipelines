#!/bin/bash

## Choose__matrices [i]
jobid[1]=tissue
jobid[2]=tissue-diet
jobid[3]=tissue-br
jobid[4]=tissue-gg
jobid[5]=tissue-br-females
jobid[6]=tissue-gg-females
jobid[7]=tissue-br-bucephalus

out=summary.txt
final=summary.br.txt
rm $final

## tissue tissue.diet tissue.gg tissue.br
for m in eXpress
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
# Get the number of genes per abundance test
cat ${f}/diffExpr*matrix.log2.dat | cut -f 1 >> raw.$m.${jobid[${t}]}.$p.$c
# count number of raw and unique differentially expressed genes
all=$(grep "^TRINITY" raw.$m.${jobid[${t}]}.$p.$c | wc -l)
raw=$(grep "^TRINITY" raw.$m.${jobid[${t}]}.$p.$c | sort - | uniq | wc -l)
paste <(printf "%s\n" "$f") <(printf "%s\n" "$all") <(printf "%s\n" "$all") >> $out
# column names; trandform to tabulated format
#cat $out | sed -e 's/.br/-br/g' -e 's/.gg/-gg/g' -e 's/.diet/-diet/g' -e 's/\./\t/g' >> $final
cat $out | sed -e 's/\./\t/g' >> $final
rm raw.$m.${jobid[${t}]}.$p.$c $out
		                done
		            done
	          done
	      done
    done
done
