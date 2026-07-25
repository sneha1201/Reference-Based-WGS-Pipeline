#!/usr/bin/env python3

import os
import glob
import subprocess
from pathlib import Path

# =====================================================
# CONFIG
# =====================================================

snpEff = "/apps/snpEff5.0/snpEff.jar"

snpEff_config = "/apps/snpEff5.0/snpEff.config"

organism = "Bacillus_altitudinis_BSA5"

base_dir = Path("/mnt/scratch/sneha/PFY2526S01WGS2777/merged")

finalvcf_dir = base_dir / "09_Final_Variants"

filtered_dir = base_dir / "08_Filtered_Variants"

annotation_dir = base_dir / "12_snpEff_Annotation"

log_dir = base_dir / "Logs"

annotation_dir.mkdir(parents=True, exist_ok=True)

# =====================================================
# GET SAMPLE NAMES
# =====================================================

samples = []

for file in glob.glob(str(finalvcf_dir / "*.final.vcf.gz")):

    sample = Path(file).name.replace(".final.vcf.gz", "")

    samples.append(sample)

print("\nSamples detected:")
print(samples)

# =====================================================
# ANNOTATE FINAL MERGED VCFs
# =====================================================

print("\nAnnotating FINAL merged VCFs")

for sample in samples:

    input_vcf = finalvcf_dir / f"{sample}.final.vcf.gz"

    output_vcf = annotation_dir / f"{sample}.final.annotated.vcf"

    csv_stats = annotation_dir / f"{sample}.final.stats.csv"

    html_stats = annotation_dir / f"{sample}.final.summary.html"

    log = log_dir / f"{sample}_snpeff_final.log"

    cmd = [
        "java", "-Xmx8g",
        "-jar", snpEff,
        "-c", snpEff_config,
        "-v", organism,
        "-noLog",
        "-no-downstream",
        "-no-upstream",
        "-no-utr",
        "-o", "vcf",
        str(input_vcf),
        "-csvStats", str(csv_stats),
        "-htmlStats", str(html_stats)
    ]

    with open(log, "w") as logfile, open(output_vcf, "w") as vcf_out:

        process = subprocess.Popen(
            cmd,
            stdout=vcf_out,
            stderr=subprocess.PIPE,
            text=True
        )

        for line in process.stderr:
            print(line, end="")
            logfile.write(line)

        process.wait()

    if process.returncode != 0:

        raise RuntimeError(
            f"snpEff annotation failed for {sample}"
        )

    print(f"Annotated FINAL VCF : {output_vcf}")

# =====================================================
# ANNOTATE PASS SNPs
# =====================================================

print("\nAnnotating PASS SNPs")

for sample in samples:

    input_vcf = filtered_dir / f"{sample}.snps.pass.vcf.gz"

    output_vcf = annotation_dir / f"{sample}.snps.annotated.vcf"

    csv_stats = annotation_dir / f"{sample}.snps.stats.csv"

    html_stats = annotation_dir / f"{sample}.snps.summary.html"

    log = log_dir / f"{sample}_snpeff_snps.log"

    cmd = [
        "java", "-Xmx8g",
        "-jar", snpEff,
        "-c", snpEff_config,
        "-v", organism,
        "-noLog",
        "-no-downstream",
        "-no-upstream",
        "-no-utr",
        "-o", "vcf",
        str(input_vcf),
        "-csvStats", str(csv_stats),
        "-htmlStats", str(html_stats)
    ]

    with open(log, "w") as logfile, open(output_vcf, "w") as vcf_out:

        process = subprocess.Popen(
            cmd,
            stdout=vcf_out,
            stderr=subprocess.PIPE,
            text=True
        )

        for line in process.stderr:
            print(line, end="")
            logfile.write(line)

        process.wait()

    if process.returncode != 0:

        raise RuntimeError(
            f"SNP annotation failed for {sample}"
        )

    print(f"Annotated SNP VCF : {output_vcf}")

# =====================================================
# ANNOTATE PASS INDELs
# =====================================================

print("\nAnnotating PASS INDELs")

for sample in samples:

    input_vcf = filtered_dir / f"{sample}.indels.pass.vcf.gz"

    output_vcf = annotation_dir / f"{sample}.indels.annotated.vcf"

    csv_stats = annotation_dir / f"{sample}.indels.stats.csv"

    html_stats = annotation_dir / f"{sample}.indels.summary.html"

    log = log_dir / f"{sample}_snpeff_indels.log"

    cmd = [
        "java", "-Xmx8g",
        "-jar", snpEff,
        "-c", snpEff_config,
        "-v", organism,
        "-noLog",
        "-no-downstream",
        "-no-upstream",
        "-no-utr",
        "-o", "vcf",
        str(input_vcf),
        "-csvStats", str(csv_stats),
        "-htmlStats", str(html_stats)
    ]

    with open(log, "w") as logfile, open(output_vcf, "w") as vcf_out:

        process = subprocess.Popen(
            cmd,
            stdout=vcf_out,
            stderr=subprocess.PIPE,
            text=True
        )

        for line in process.stderr:
            print(line, end="")
            logfile.write(line)

        process.wait()

    if process.returncode != 0:

        raise RuntimeError(
            f"INDEL annotation failed for {sample}"
        )

    print(f"Annotated INDEL VCF : {output_vcf}")

print("\nALL snpEff ANNOTATIONS COMPLETED")
