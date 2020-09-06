#!/bin/bash

usage="
USAGE: This script process the BAM file and generate exon and gene level RPKM.
sh 01_bam_to_rpkm.sh <sampleId> <sample.bam> <mapped_reads> <genes.bed> <exons.bed>
"

if [ $# -ne 5 ] ; then printf "Error: No arguments\n${usage}" >&2 ; exit 1 ; fi

sampleId=$1			## sample ID which will be used as output file prefix
file_bam=$2			## bam file
mappedReads=$3		## number of mapped reads
file_genes=$4		## bed file for gene
file_exons=$5		## bed file for exons


## Get exon coverage
bedtools2-2.19.0/bin/coverageBed -abam ${file_bam} -b ${file_exons} > ${sampleId}.cov.txt

## Get exonwise rpkm
cat ${sampleId}.cov.txt | awk -v readCount="$mappedReads" '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t",($5*10^9)/(($3-$2) * readCount)}' > ${sampleId}.exons.rpkm.txt

## Get genewise rpkm
perl 01_exon_to_gene_rpkm.pl ${file_genes} ${sampleId}.exons.rpkm.txt >> ${sampleId}.genes.rpkm.txt
awk -v sample="$sampleId" 'BEGIN{print "chr\tstart\tend\tgene\t"sample}{print $1"\t"$2"\t"$3"\t"$4"\t"$5}' ${sampleId}.genes.rpkm.txt >> ${sampleId}.txt

