#!/usr/bin/env python3
"""
NGS Variant Calling Pipeline
FASTQ → QC → Trimming → Alignment → Variant Calling → Filtered VCF
Author: Taruna Gupta
GitHub: github.com/tarunagupta20
"""

import os
import sys
import subprocess
import argparse
import logging

# ── Logging setup ──────────────────────────────────────────────────────────────
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.StreamHandler(sys.stdout),
        logging.FileHandler("pipeline.log")
    ]
)
log = logging.getLogger(__name__)

# ── Helper function ────────────────────────────────────────────────────────────
def run(cmd, step_name):
    log.info(f"Starting: {step_name}")
    log.info(f"Command:  {cmd}")
    result = subprocess.run(cmd, shell=True)
    if result.returncode != 0:
        log.error(f"FAILED at step: {step_name}")
        log.error("Check the error above, fix it, and re-run.")
        sys.exit(1)
    log.info(f"Completed: {step_name}\n")

# ── Pipeline ───────────────────────────────────────────────────────────────────
def main():
    parser = argparse.ArgumentParser(
        description="End-to-end NGS Variant Calling Pipeline",
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument("--r1",      required=True,       help="Path to Read 1 FASTQ file")
    parser.add_argument("--r2",      required=True,       help="Path to Read 2 FASTQ file")
    parser.add_argument("--ref",     required=True,       help="Path to reference genome FASTA")
    parser.add_argument("--sample",  required=True,       help="Sample name (used for output files)")
    parser.add_argument("--outdir",  default="results",   help="Output directory (default: results/)")
    parser.add_argument("--threads", default=4, type=int, help="Number of threads (default: 4)")
    parser.add_argument("--picard",  default="/usr/share/java/picard.jar", help="Path to picard.jar")
    parser.add_argument("--gatk",   default="/home/taruna/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar", help="Path to GATK jar")
    args = parser.parse_args()

    # Output subdirectories
    qc_dir      = f"{args.outdir}/qc"
    trim_dir    = f"{args.outdir}/trimmed"
    align_dir   = f"{args.outdir}/aligned"
    variant_dir = f"{args.outdir}/variants"

    for d in [qc_dir, trim_dir, align_dir, variant_dir]:
        os.makedirs(d, exist_ok=True)

    s   = args.sample
    ref = args.ref
    t   = args.threads
    pic = args.picard
    gatk = args.gatk

    log.info("=" * 60)
    log.info(f"  NGS Variant Calling Pipeline")
    log.info(f"  Sample:    {s}")
    log.info(f"  Reference: {ref}")
    log.info(f"  Threads:   {t}")
    log.info("=" * 60 + "\n")

    # ── STEP 1: Quality Control ──────────────────────────────────────────────
    run(
        f"fastqc {args.r1} {args.r2} -o {qc_dir} -t {t}",
        "Step 1: FastQC — Quality Control"
    )

    # ── STEP 2: Trimming ─────────────────────────────────────────────────────
    r1_trim   = f"{trim_dir}/{s}_R1_trimmed.fastq.gz"
    r2_trim   = f"{trim_dir}/{s}_R2_trimmed.fastq.gz"

    run(
        f"fastp "
        f"-i {args.r1} -I {args.r2} "
        f"-o {r1_trim} -O {r2_trim} "        
        f"--thread {t} "
        f"--detect_adapter_for_pe "
        f"--qualified_quality_phred 15 "
        f"--length_required 36 "
        f"--json {qc_dir}/{s}_fastp.json "
        f"--html {qc_dir}/{s}_fastp.html",
        "Step 2: fastp — Adapter Trimming & QC"
    )
    # ── STEP 3: Alignment ────────────────────────────────────────────────────
    raw_bam    = f"{align_dir}/{s}_raw.bam"
    sorted_bam = f"{align_dir}/{s}_sorted.bam"

    run(
        f"bwa mem -t {t} "
        f"-R '@RG\\tID:{s}\\tSM:{s}\\tPL:ILLUMINA\\tLB:lib1\\tPU:unit1' "
        f"{ref} {r1_trim} {r2_trim} | "
        f"samtools view -Sb - > {raw_bam}",
        "Step 3: BWA-MEM — Read Alignment"
    )

    # ── STEP 4: Sort & Index BAM ─────────────────────────────────────────────
    run(
        f"samtools sort -@ {t} -o {sorted_bam} {raw_bam} && "
        f"samtools index {sorted_bam}",
        "Step 4: SAMtools — Sort & Index BAM"
    )

    # Clean up unsorted BAM to save space
    os.remove(raw_bam)
    log.info(f"Removed intermediate file: {raw_bam}\n")

    # ── STEP 5: Mark Duplicates ──────────────────────────────────────────────
    dedup_bam = f"{align_dir}/{s}_dedup.bam"
    metrics   = f"{align_dir}/{s}_dup_metrics.txt"

    run(
        f"java -jar {pic} MarkDuplicates "
        f"I={sorted_bam} O={dedup_bam} M={metrics} "
        f"REMOVE_DUPLICATES=false CREATE_INDEX=true "
        f"VALIDATION_STRINGENCY=SILENT",
        "Step 5: Picard — Mark Duplicates"
    )

    # ── STEP 6: Alignment Stats ──────────────────────────────────────────────
    run(
        f"samtools flagstat {dedup_bam} > {align_dir}/{s}_flagstat.txt",
        "Step 6: SAMtools — Alignment Statistics"
    )

    # ── STEP 7: Variant Calling ──────────────────────────────────────────────
    raw_vcf = f"{variant_dir}/{s}_raw.vcf.gz"

    run(
        f"java -jar {gatk} HaplotypeCaller "
        f"-R {ref} "
        f"-I {dedup_bam} "
        f"-O {raw_vcf} "
        f"--sample-name {s}",
        "Step 7: GATK HaplotypeCaller — Variant Calling"
    )

    # ── STEP 8: Filter Variants ──────────────────────────────────────────────
    filtered_vcf = f"{variant_dir}/{s}_filtered.vcf.gz"

    run(
        f"java -jar {gatk} VariantFiltration "
        f"-R {ref} "
        f"-V {raw_vcf} "
        f"--filter-expression 'QD < 2.0'   --filter-name 'QD2' "
        f"--filter-expression 'FS > 60.0'  --filter-name 'FS60' "
        f"--filter-expression 'MQ < 40.0'  --filter-name 'MQ40' "
        f"-O {filtered_vcf}",
        "Step 8: GATK VariantFiltration — Hard Filtering"
    )

    # ── STEP 9: Summary Stats ────────────────────────────────────────────────
    run(
        f"bcftools stats {filtered_vcf} > {variant_dir}/{s}_stats.txt",
        "Step 9: BCFtools — Variant Statistics"
    )

    # ── Print final summary ──────────────────────────────────────────────────
    log.info("=" * 60)
    log.info("  PIPELINE COMPLETE")
    log.info(f"  Final VCF:        {filtered_vcf}")
    log.info(f"  Variant stats:    {variant_dir}/{s}_stats.txt")
    log.info(f"  Alignment stats:  {align_dir}/{s}_flagstat.txt")
    log.info(f"  Dup metrics:      {metrics}")
    log.info(f"  FastQC reports:   {qc_dir}/")
    log.info("=" * 60)

if __name__ == "__main__":
    main()
