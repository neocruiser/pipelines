#!/bin/bash
#PBS -q vsmp
#PBS -l nodes=1:ppn=256:vsmp
#PBS -l walltime=48:00:00
#PBS -N trinity_all_btfy
#PBS -e trinity_all.err
#PBS -A sun108
#PBS -M sleiman.bassim@stonybrook.edu
#PBS -m abe
#PBS -V

# set stack to unlimited
# because of large datasets
ulimit -s unlimited
# echo stdout to output file
set -x
# xsede directories
oasis=/oasis/projects/nsf/sun108
scratch=/oasis/scratch/silo/temp_project
home=/home/silo
##################################
# output directories for trinity #
##################################
jobid=all
fileid=raw.all

########################
# Dont change anything #
########################
workdir=${scratch}/ganglia/trinity/trinity_out_dir_${jobid}/
mkdir -p ${workdir}
cd ${workdir}
# fastq raw files, reads
sense=${scratch}/ganglia/data/$fileid.R1.fq.gz
antisense=${scratch}/ganglia/data/$fileid.R2.fq.gz
# maximum memory for jellyfish to use
JM=900G

###########
# TRINITY #
###########
## Phase 1: construct Kmers
#--no_run_chrysalis
## Phase 2: cluster reads, DeBruijn, and map to kmers
${home}/trinity/Trinity --seqType fq \
--SS_lib_type RF \
--left ${sense} \
--right ${antisense} \
--quality_trimming_params "ILLUMINACLIP:adapters.fa:2:30:10 LEADING:5 TRAILING:5 SLIDINGWINDOW:4:15 MINLEN:36" \
--normalize_max_read_cov 50 \
--min_contig_length 200 \
--output ${workdir} \
--max_memory ${JM} \
--CPU 64 \
--bflyCPU 64 \
--bflyGCThreads 30 \
--no_run_chrysalis \
>& ${home}/trinity.${jobid}_output.log