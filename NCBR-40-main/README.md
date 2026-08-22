<div align="center">
   
  <h1>NCBR-40 🔬</h1>
  
  **_Data and script for processing NCBR-40_**

</div>

# Overview

This repository contains data and scripts for analyzing NCBR-40. 

Any data accompanying this project can be stored in the `data/` directory, while any scripts used to process the data can be stored in the `scripts/` directory.

> [!NOTE]  
> _**By default**, any data or files added to the `data/` directory are ignored by git due to our `.gitignore`, so you can store large files or data here without worrying about them being uploaded to Github._ If you would to upload a small data file (<5MB) to Github, you can stage the file using force, `-f`, option to the `git add`:  
> ```bash
> # Stage file for commit
> git add -f data/counts.tsv
> # Commit the file to history
> git commit -m "Adding small counts matrix"
> ```

## Installation

To install the repository locally, you can use the following command:

```bash
# Clone the github repository
# and change your working directory
git clone https://github.com/OpenOmics/NCBR-40.git
cd NCBR-40/

# Download the per-sample RSEM
# isoform counts from Box,
# for H-T samples
rclone copy -v \
  "Box:/NCBR-40/NCBR-40_IsoformSwitchAnalyzeR_2024/IsoformSwitchAnalyzeR/Inputs/RSEM/H-T/Isoforms/" \
  data/H-T/rsem_isoforms/
# for the U-T samples
rclone copy -v \
  "Box:/NCBR-40/NCBR-40_IsoformSwitchAnalyzeR_2024/IsoformSwitchAnalyzeR/Inputs/RSEM/U-T/Isoforms/" \
  data/U-T/rsem_isoforms/
# Download the per-sample RSEM
# isoform counts from Box, for
# the healthy/patient samples
# from NCBR-390. These samples
# use the same H-T transcriptome.
rclone copy -v \
  "Box:/NCBR-40/NCBR-40_IsoformSwitchAnalyzeR_2024/IsoformSwitchAnalyzeR/Inputs/RSEM/NCBR-390/H-T/Isoforms/" \
  data/NCBR-390/H-T/rsem_isoforms/

# Download the GTF and transcriptomic
# fasta files from Box for CDS/ORF 
# prediction
# GTF file created from H-T samples
rclone copy -v \
  "Box:/NCBR-40/NCBR-40_IsoformSwitchAnalyzeR_2024/IsoformSwitchAnalyzeR/Inputs/refs/GRCh38_Gencode_CHR_v40_plus_H-T_Isoseq.gtf" \
  data/
# Transcriptomic FASTA file created
# from H-T samples 
rclone copy -v \
  "Box:/NCBR-40/NCBR-40_IsoformSwitchAnalyzeR_2024/IsoformSwitchAnalyzeR/Inputs/refs/transcripts_H-T.fa" \
  data/
# GTF file created from U-T samples
rclone copy -v \
  "Box:/NCBR-40/NCBR-40_IsoformSwitchAnalyzeR_2024/IsoformSwitchAnalyzeR/Inputs/refs/GRCh38_Gencode_CHR_v40_plus_U-T_Isoseq.gtf" \
  data/
# Transcriptomic FASTA file created
# from U-T samples
rclone copy -v \
  "Box:/NCBR-40/NCBR-40_IsoformSwitchAnalyzeR_2024/IsoformSwitchAnalyzeR/Inputs/refs/transcripts_U-T.fa" \
  data/
```

## Setup your environment

To setup your environment and download any missing packages, you can use the following command:

```bash
# Install any missing or 
# required python packages
# in a virtual environment
python -m venv .venv
source .venv/bin/activate
pip install -U pip
pip install -r requirements.txt

# Install any missing R packages
./packages.R
```

## Reproduce the analyses

This is where you can add any steps to reproduce the analyses. For example, you can add the following command to run the script:

```bash
# Run IsoformSwitchAnalyzeR to characterize
# any differential isoform switching and 
# alternative splicing events using the
# H-T and U-T sample transcriptome
module load R/4.4
# Results are filtered via:
#   • adjusted p-value <= 0.05
#   • abs(dIF) >= 0.1
#   • IF >= 0.01
./scripts/IsoformSwitchAnalyzeR.R \
  -i data/H-T/rsem_isoforms/ \
  -o results/H-T/IsoformSwitchAnalyzeR/ \
  -g data/GRCh38_Gencode_CHR_v40_plus_H-T_Isoseq.gtf \
  -s data/sample_sheet_isa_H-T.tsv \
  -c1 "H_minus" -c2 "H_plus" \
  -t data/transcripts_H-T.fa
./scripts/IsoformSwitchAnalyzeR.R \
  -i data/H-T/rsem_isoforms/ \
  -o results/H-T/IsoformSwitchAnalyzeR/ \
  -g data/GRCh38_Gencode_CHR_v40_plus_H-T_Isoseq.gtf \
  -s data/sample_sheet_isa_H-T.tsv \
  -c1 "T_minus" -c2 "T_plus" \
  -t data/transcripts_H-T.fa
# Using the U-T sample transcriptome
./scripts/IsoformSwitchAnalyzeR.R \
  -i data/U-T/rsem_isoforms/ \
  -o results/U-T/IsoformSwitchAnalyzeR/ \
  -g data/GRCh38_Gencode_CHR_v40_plus_U-T_Isoseq.gtf \
  -s data/sample_sheet_isa_U-T.tsv \
  -c1 "U_minus" -c2 "U_plus" \
  -t data/transcripts_U-T.fa
./scripts/IsoformSwitchAnalyzeR.R \
  -i data/U-T/rsem_isoforms/ \
  -o results/U-T/IsoformSwitchAnalyzeR/ \
  -g data/GRCh38_Gencode_CHR_v40_plus_U-T_Isoseq.gtf \
  -s data/sample_sheet_isa_U-T.tsv \
  -c1 "T_minus" -c2 "T_plus" \
  -t data/transcripts_U-T.fa

# Get unfiltered results, i.e:
#   • adjusted p-value <= 1.0
#   • abs(dIF) >= 0.0
#   • IF >= 0.0
./scripts/IsoformSwitchAnalyzeR.R \
  -i data/H-T/rsem_isoforms/ \
  -o results/H-T/IsoformSwitchAnalyzeR_unfiltered/ \
  -g data/GRCh38_Gencode_CHR_v40_plus_H-T_Isoseq.gtf \
  -s data/sample_sheet_isa_H-T.tsv \
  -c1 "H_minus" -c2 "H_plus" \
  -t data/transcripts_H-T.fa \
  -p 1.0 -f 0.0 -u 0.0 
./scripts/IsoformSwitchAnalyzeR.R \
  -i data/H-T/rsem_isoforms/ \
  -o results/H-T/IsoformSwitchAnalyzeR_unfiltered/ \
  -g data/GRCh38_Gencode_CHR_v40_plus_H-T_Isoseq.gtf \
  -s data/sample_sheet_isa_H-T.tsv \
  -c1 "T_minus" -c2 "T_plus" \
  -t data/transcripts_H-T.fa \
  -p 1.0 -f 0.0 -u 0.0
# Using the U-T sample transcriptome
./scripts/IsoformSwitchAnalyzeR.R \
  -i data/U-T/rsem_isoforms/ \
  -o results/U-T/IsoformSwitchAnalyzeR_unfiltered/ \
  -g data/GRCh38_Gencode_CHR_v40_plus_U-T_Isoseq.gtf \
  -s data/sample_sheet_isa_U-T.tsv \
  -c1 "U_minus" -c2 "U_plus" \
  -t data/transcripts_U-T.fa \
  -p 1.0 -f 0.0 -u 0.0
./scripts/IsoformSwitchAnalyzeR.R \
  -i data/U-T/rsem_isoforms/ \
  -o results/U-T/IsoformSwitchAnalyzeR_unfiltered/ \
  -g data/GRCh38_Gencode_CHR_v40_plus_U-T_Isoseq.gtf \
  -s data/sample_sheet_isa_U-T.tsv \
  -c1 "T_minus" -c2 "T_plus" \
  -t data/transcripts_U-T.fa \
  -p 1.0 -f 0.0 -u 0.0
```

**NCBR-390 ISA comparisons**

Using the same H-T transcriptome, IsoformSwitchAnalyzeR was on NCBR-390 patient/heathly samples.

```bash
# Run IsoformSwitchAnalyzeR to characterize
# any differential isoform switching and 
# alternative splicing events using the
# H-T sample transcriptome with samples
# from NCBR-390
module load R/4.4
# Results are filtered via:
#   • adjusted p-value <= 0.05
#   • abs(dIF) >= 0.1
#   • IF >= 0.01
./scripts/IsoformSwitchAnalyzeR.R \
  -i data/NCBR-390/H-T/rsem_isoforms/ \
  -o results/NCBR-390/H-T/IsoformSwitchAnalyzeR/ \
  -g data/GRCh38_Gencode_CHR_v40_plus_H-T_Isoseq.gtf \
  -s data/sample_sheet_isa_NCBR-390_H-T.tsv \
  -c1 "grpA" -c2 "grpB" \
  -t data/transcripts_H-T.fa
./scripts/IsoformSwitchAnalyzeR.R \
  -i data/NCBR-390/H-T/rsem_isoforms/ \
  -o results/NCBR-390/H-T/IsoformSwitchAnalyzeR/ \
  -g data/GRCh38_Gencode_CHR_v40_plus_H-T_Isoseq.gtf \
  -s data/sample_sheet_isa_NCBR-390_H-T.tsv \
  -c1 "grpC" -c2 "grpD" \
  -t data/transcripts_H-T.fa
./scripts/IsoformSwitchAnalyzeR.R \
  -i data/NCBR-390/H-T/rsem_isoforms/ \
  -o results/NCBR-390/H-T/IsoformSwitchAnalyzeR/ \
  -g data/GRCh38_Gencode_CHR_v40_plus_H-T_Isoseq.gtf \
  -s data/sample_sheet_isa_NCBR-390_H-T.tsv \
  -c1 "grpE" -c2 "grpF" \
  -t data/transcripts_H-T.fa
./scripts/IsoformSwitchAnalyzeR.R \
  -i data/NCBR-390/H-T/rsem_isoforms/ \
  -o results/NCBR-390/H-T/IsoformSwitchAnalyzeR/ \
  -g data/GRCh38_Gencode_CHR_v40_plus_H-T_Isoseq.gtf \
  -s data/sample_sheet_isa_NCBR-390_H-T.tsv \
  -c1 "grpG" -c2 "grpH" \
  -t data/transcripts_H-T.fa
./scripts/IsoformSwitchAnalyzeR.R \
  -i data/NCBR-390/H-T/rsem_isoforms/ \
  -o results/NCBR-390/H-T/IsoformSwitchAnalyzeR/ \
  -g data/GRCh38_Gencode_CHR_v40_plus_H-T_Isoseq.gtf \
  -s data/sample_sheet_isa_NCBR-390_H-T.tsv \
  -c1 "grpI" -c2 "grpJ" \
  -t data/transcripts_H-T.fa
./scripts/IsoformSwitchAnalyzeR.R \
  -i data/NCBR-390/H-T/rsem_isoforms/ \
  -o results/NCBR-390/H-T/IsoformSwitchAnalyzeR/ \
  -g data/GRCh38_Gencode_CHR_v40_plus_H-T_Isoseq.gtf \
  -s data/sample_sheet_isa_NCBR-390_H-T.tsv \
  -c1 "grpK" -c2 "grpL" \
  -t data/transcripts_H-T.fa
./scripts/IsoformSwitchAnalyzeR.R \
  -i data/NCBR-390/H-T/rsem_isoforms/ \
  -o results/NCBR-390/H-T/IsoformSwitchAnalyzeR/ \
  -g data/GRCh38_Gencode_CHR_v40_plus_H-T_Isoseq.gtf \
  -s data/sample_sheet_isa_NCBR-390_H-T.tsv \
  -c1 "grpM" -c2 "grpN" \
  -t data/transcripts_H-T.fa
./scripts/IsoformSwitchAnalyzeR.R \
  -i data/NCBR-390/H-T/rsem_isoforms/ \
  -o results/NCBR-390/H-T/IsoformSwitchAnalyzeR/ \
  -g data/GRCh38_Gencode_CHR_v40_plus_H-T_Isoseq.gtf \
  -s data/sample_sheet_isa_NCBR-390_H-T.tsv \
  -c1 "grpO" -c2 "grpP" \
  -t data/transcripts_H-T.fa

# Get unfiltered results, i.e:
#   • adjusted p-value <= 1.0
#   • abs(dIF) >= 0.0
#   • IF >= 0.0
./scripts/IsoformSwitchAnalyzeR.R \
  -i data/NCBR-390/H-T/rsem_isoforms/ \
  -o results/NCBR-390/H-T/IsoformSwitchAnalyzeR_unfiltered/ \
  -g data/GRCh38_Gencode_CHR_v40_plus_H-T_Isoseq.gtf \
  -s data/sample_sheet_isa_NCBR-390_H-T.tsv \
  -c1 "grpA" -c2 "grpB" \
  -t data/transcripts_H-T.fa \
  -p 1.0 -f 0.0 -u 0.0
./scripts/IsoformSwitchAnalyzeR.R \
  -i data/NCBR-390/H-T/rsem_isoforms/ \
  -o results/NCBR-390/H-T/IsoformSwitchAnalyzeR_unfiltered/ \
  -g data/GRCh38_Gencode_CHR_v40_plus_H-T_Isoseq.gtf \
  -s data/sample_sheet_isa_NCBR-390_H-T.tsv \
  -c1 "grpC" -c2 "grpD" \
  -t data/transcripts_H-T.fa \
  -p 1.0 -f 0.0 -u 0.0
./scripts/IsoformSwitchAnalyzeR.R \
  -i data/NCBR-390/H-T/rsem_isoforms/ \
  -o results/NCBR-390/H-T/IsoformSwitchAnalyzeR_unfiltered/ \
  -g data/GRCh38_Gencode_CHR_v40_plus_H-T_Isoseq.gtf \
  -s data/sample_sheet_isa_NCBR-390_H-T.tsv \
  -c1 "grpE" -c2 "grpF" \
  -t data/transcripts_H-T.fa \
  -p 1.0 -f 0.0 -u 0.0
./scripts/IsoformSwitchAnalyzeR.R \
  -i data/NCBR-390/H-T/rsem_isoforms/ \
  -o results/NCBR-390/H-T/IsoformSwitchAnalyzeR_unfiltered/ \
  -g data/GRCh38_Gencode_CHR_v40_plus_H-T_Isoseq.gtf \
  -s data/sample_sheet_isa_NCBR-390_H-T.tsv \
  -c1 "grpG" -c2 "grpH" \
  -t data/transcripts_H-T.fa \
  -p 1.0 -f 0.0 -u 0.0
./scripts/IsoformSwitchAnalyzeR.R \
  -i data/NCBR-390/H-T/rsem_isoforms/ \
  -o results/NCBR-390/H-T/IsoformSwitchAnalyzeR_unfiltered/ \
  -g data/GRCh38_Gencode_CHR_v40_plus_H-T_Isoseq.gtf \
  -s data/sample_sheet_isa_NCBR-390_H-T.tsv \
  -c1 "grpI" -c2 "grpJ" \
  -t data/transcripts_H-T.fa \
  -p 1.0 -f 0.0 -u 0.0
./scripts/IsoformSwitchAnalyzeR.R \
  -i data/NCBR-390/H-T/rsem_isoforms/ \
  -o results/NCBR-390/H-T/IsoformSwitchAnalyzeR_unfiltered/ \
  -g data/GRCh38_Gencode_CHR_v40_plus_H-T_Isoseq.gtf \
  -s data/sample_sheet_isa_NCBR-390_H-T.tsv \
  -c1 "grpK" -c2 "grpL" \
  -t data/transcripts_H-T.fa \
  -p 1.0 -f 0.0 -u 0.0
./scripts/IsoformSwitchAnalyzeR.R \
  -i data/NCBR-390/H-T/rsem_isoforms/ \
  -o results/NCBR-390/H-T/IsoformSwitchAnalyzeR_unfiltered/ \
  -g data/GRCh38_Gencode_CHR_v40_plus_H-T_Isoseq.gtf \
  -s data/sample_sheet_isa_NCBR-390_H-T.tsv \
  -c1 "grpM" -c2 "grpN" \
  -t data/transcripts_H-T.fa \
  -p 1.0 -f 0.0 -u 0.0
./scripts/IsoformSwitchAnalyzeR.R \
  -i data/NCBR-390/H-T/rsem_isoforms/ \
  -o results/NCBR-390/H-T/IsoformSwitchAnalyzeR_unfiltered/ \
  -g data/GRCh38_Gencode_CHR_v40_plus_H-T_Isoseq.gtf \
  -s data/sample_sheet_isa_NCBR-390_H-T.tsv \
  -c1 "grpO" -c2 "grpP" \
  -t data/transcripts_H-T.fa \
  -p 1.0 -f 0.0 -u 0.0
```
