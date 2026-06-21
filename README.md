# NGS Variant Calling Pipeline

An end-to-end NGS pipeline for germline variant detection, validated on the NA12878/HG001 GIAB gold standard sample (chr20).

## Tools Used
- FastQC / fastp — quality control and trimming
- BWA-MEM — read alignment to hg38
- SAMtools — BAM sorting and indexing
- Picard — duplicate marking
- GATK HaplotypeCaller — variant calling
- bcftools — VCF filtering and stats

## Pipeline Steps
1. Quality control (FastQC → fastp)
2. Alignment to hg38 reference (BWA-MEM)
3. Sort and index BAM (SAMtools)
4. Mark duplicates (Picard)
5. Variant calling (GATK HaplotypeCaller)
6. VCF filtering (bcftools)

## Validation Results (NA12878 / HG001 — chr20)
| Metric | Value |
|--------|-------|
| Total reads processed | 4,172,718 |
| Mapped reads | 241,854 (5.78% — chr20 only) |
| Duplicate reads | 4,815 |
| SNPs detected | 9,953 |
| Indels detected | 1,912 |

> Note: Low mapping rate is expected — reads are whole-genome but aligned to chr20 only.

## How to Run
```bash
python scripts/run_pipeline.py
```

## Author
Taruna Gupta | MSc Bioinformatics | Mumbai University
