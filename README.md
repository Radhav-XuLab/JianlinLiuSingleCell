# Single cell exome sequencing data analysis scripts
This repository has scripts related to single cell exome sequencing project. The work is published ...

### 01_raw_data_process.sh
This workflow was used to align the raw fastq data against reference genome using `bwa` followed by GATK Best practice for variant detection.

### 02_bam_to_rpkm.sh
This script generates genewise rpkm from bam files  
**Usage:** `sh 01_bam_to_rpkm.sh <sampleId> <sample.bam> <mapped_reads> <genes.bed> <exons.bed>`  
**Description:**  
- ``<sampleId>``: BED file with gene coordinates
- ``<sample.bam>``: BAM file processed using GATK Best practice workflow for variant detection.
- ``<mapped_reads>``: number of mapped reads for RPKM calculation
- ``<genes.bed>``: BED file with gene coordinates
- ``<exons.bed>``: BED file with exon coordinates  
**Output**: sample.cov.txt, sample.rpkm.txt



