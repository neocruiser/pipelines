#!/bin/bash
#SBATCH --partition=RM
#SBATCH --nodes=1
#SBATCH -t 08:00:00
#SBATCH --job-name="annot.SUM"
#SBATCH --output="annot.summary.%j.%N.out"
#SBATCH --mail-type=ALL
#SBATCH --mail-user=sleiman.bassim@stonybrook.edu

## Append NCBIs gene information, gene ID, protein accession nb, and taxid to each contig
scratch=/pylon2/oc4ifip/bassim
home=/home/bassim

gene_info=/pylon1/oc4ifip/bassim/db/ncbi/gene_info
gene2accession=/pylon1/oc4ifip/bassim/db/ncbi/gene2accession

_DIR=$scratch/ganglia/trinity/trinity_out_dir_00000/lectins
_diamond=$_DIR/diamond.contigIDs-lectins.p4.c2.tissue.nr.blastx.00000.txt

## e-value for blast search
ED=5
ZD=$(seq -s. "$(echo "${ED}+1" | bc)" | tr -d '[:digit:]' | sed 's/./0/g')
_ev="0.${ZD}1"

_addon=$(echo $(basename $_diamond))
_diamond_summary=$_DIR/NR-PTHR-IPS.$_addon.tmp

awk 'NR==FNR {b[$7] = sprintf ("%s\t%s\t%s\t%s\t%s\t%s\t%s\t",$1,$2,$3,$4,$5,$6,$7); next} {print b[$2],$1,$2}' <(awk 'NR==FNR {a[$2] = sprintf ("%s\t%s\t%s\t%s\t%s\t",$1,$2,$3,$4,$5); next} {print a[$1],$1,$2}' <(cat $gene_info | grep -Fwf <(cat $gene2accession | cut -f2,6 | grep -Fwf <(cat $_diamond  | cut -f1,2,11 | awk -ve="$_ev" '{if($3<=e)print$0}' | cut -f1,2 | sed 's/|/\t/g' | cut -f5) - | cut -f1) - | cut -f1,2,3,9,10 | sed 's/ /./g') <(cat $gene2accession | cut -f2,6 | grep -Fwf <(cat $_diamond  | cut -f1,2,11 | awk -ve="$_ev" '{if($3<=e)print$0}' | cut -f1,2 | sed 's/|/\t/g' | cut -f5) -)) <(cat $_diamond  | cut -f1,2,11 | awk -ve="$_ev" '{if($3<=e)print$0}' | cut -f1,2 | sed 's/|/\t/g' | cut -f1,5) | awk '{if ($2 == $6 && $7 == $9) print $1,$2,$3,$4,$5,$7,$8}' | sort - | uniq | cut -f1,4,7 -d ' ' | sort -k2 | uniq  > $_diamond_summary
