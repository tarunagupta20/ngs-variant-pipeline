# Changelog

## [1.1.0] - 2026-06-21
### Changed
- Removed hardcoded server-specific GATK path; now requires --gatk argument

## [1.0.0] - 2026-05-23
### Added
- End-to-end pipeline: FASTQ to filtered VCF
- 9-step workflow: FastQC, fastp, BWA-MEM, SAMtools, Picard, GATK, bcftools
- Validated on NA12878/HG001 GIAB gold standard (chr20)
- Logging to stdout and pipeline.log
