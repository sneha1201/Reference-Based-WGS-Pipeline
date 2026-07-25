# Reference-Based-WGS-Pipeline
Reference-Based-WGS-Pipeline
# Reference-Based WGS Variant Calling Pipeline

## Overview

This pipeline performs end-to-end Whole Genome Sequencing (WGS) analysis from paired-end FASTQ files to filtered variant calls.

## Workflow

FASTQ

↓

FastQC

↓

Fastp

↓

FastQC

↓

SeqKit Statistics

↓

BWA-MEM Alignment

↓

Samtools Sort

↓

Picard MarkDuplicates

↓

Samtools Flagstat

↓

Mosdepth Coverage

↓

GATK HaplotypeCaller

↓

SNP & INDEL Separation

↓

Variant Filtering

↓

Merge PASS Variants

↓

MultiQC Report

## Software

- FastQC
- Fastp
- SeqKit
- BWA
- Samtools
- Picard
- GATK 4
- Mosdepth
- MultiQC


## Usage

python pipeline/wgs_pipeline.py \
-i fastq/ \
-r reference/hg38 \
-o output \
-t 16

## Output

01_Raw_QC/

02_Clean_QC/

03_Seqkit_Stats/

04_Trimmed_Reads/

05_Alignment/

06_Coverage/

07_Variant_Calling/

08_Filtered_Variants/

09_Final_Variants/

10_MultiQC_Report/

11_Flagstat/

## Author

Sneha Goel
