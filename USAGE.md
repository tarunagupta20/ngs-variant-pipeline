# How to Run This Pipeline

## Requirements
- Python 3.8+
- BWA, SAMtools, Picard, GATK 4.x, bcftools
- Minimum 16GB RAM, 100GB disk space

## Setup
Clone this repository:
```bash
git clone https://github.com/tarunagupta20/ngs-variant-pipeline.git
cd ngs-variant-pipeline
```

## Input
- Paired-end FASTQ files (R1 and R2)
- hg38 reference genome

## Run
```bash
python scripts/run_pipeline.py \
  --r1 data/raw/sample_R1.fastq \
  --r2 data/raw/sample_R2.fastq \
  --ref data/reference/hg38.fa \
  --out results/
```

## Output
- `results/variants/` — filtered VCF file
- `results/qc/` — FastQC and fastp reports
- `results/aligned/` — sorted, deduplicated BAM
