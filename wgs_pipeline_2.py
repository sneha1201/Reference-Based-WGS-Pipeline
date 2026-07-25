#!/usr/bin/env python3

import os
import glob
import argparse
import subprocess
import datetime
from pathlib import Path


# =========================================================
#        REFERENCE BASED WGS PIPELINE
# =========================================================

class WGSPipeline:

    def __init__(self, input_dir, reference_prefix, output_dir, threads):

        self.input_dir = Path(input_dir).resolve()

        # Reference prefix WITHOUT extension
        self.ref_prefix = Path(reference_prefix).resolve()

        self.reference = self.detect_reference()

        self.output_dir = Path(output_dir).resolve()

        self.threads = threads

        self.picard = "/apps/picard/build/libs/picard.jar"
        self.gatk = "/apps/gatk-4.2.6.1/gatk"

        self.paths = self.make_directories()

        self.samples = self.find_samples()

    # =====================================================
    # DETECT REFERENCE FASTA
    # =====================================================

    def detect_reference(self):

        extensions = [".fa", ".fasta", ".fna"]

        for ext in extensions:

            fasta = str(self.ref_prefix) + ext

            if os.path.exists(fasta):
                return Path(fasta)

        raise FileNotFoundError(
            f"Reference genome not found for prefix: {self.ref_prefix}"
        )

    # =====================================================
    # CREATE OUTPUT DIRECTORIES
    # =====================================================

    def make_directories(self):

        folders = {

            "rawqc": self.output_dir / "01_Raw_QC",

            "cleanqc": self.output_dir / "02_Clean_QC",

            "seqkit": self.output_dir / "03_Seqkit_Stats",

            "trimmed": self.output_dir / "04_Trimmed_Reads",

            "alignment": self.output_dir / "05_Alignment",

            "coverage": self.output_dir / "06_Coverage",

            "variants": self.output_dir / "07_Variant_Calling",

            "filtered": self.output_dir / "08_Filtered_Variants",

            "finalvcf": self.output_dir / "09_Final_Variants",

            "multiqc": self.output_dir / "10_MultiQC_Report",

            "flagstat": self.output_dir / "11_Flagstat",

            "logs": self.output_dir / "Logs"
        }

        for path in folders.values():
            path.mkdir(parents=True, exist_ok=True)

        return folders

    # =====================================================
    # SAMPLE IDENTIFICATION
    # =====================================================

    def find_samples(self):

        files = glob.glob(
            str(self.input_dir / "*_R1_001.fastq.gz")
        )

        samples = []

        for f in files:

            sample = Path(f).name.replace(
                "_R1_001.fastq.gz",
                ""
            )

            samples.append(sample)

        return sorted(samples)

    # =====================================================
    # COMMAND RUNNER
    # =====================================================

    def run_command(self, cmd, logfile):

        with open(logfile, "w") as log:

            process = subprocess.Popen(
                cmd,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                text=True
            )

            for line in process.stdout:
                print(line, end="")
                log.write(line)

            process.wait()

        if process.returncode != 0:
            raise RuntimeError(f"FAILED : {logfile}")

    # =====================================================
    # FETCH FASTQ FILES
    # =====================================================

    def reads(self, sample):

        r1 = self.input_dir / f"{sample}_R1_001.fastq.gz"

        r2 = self.input_dir / f"{sample}_R2_001.fastq.gz"

        return r1, r2

    # =====================================================
    # STEP 1 : RAW FASTQC
    # =====================================================

    def raw_fastqc(self):

        print("\nSTEP 1 : RAW FASTQC")

        for sample in self.samples:

            r1, r2 = self.reads(sample)

            cmd = [
                "fastqc",
                "-t", str(self.threads),
                str(r1),
                str(r2),
                "-o",
                str(self.paths["rawqc"])
            ]

            self.run_command(
                cmd,
                self.paths["logs"] / f"{sample}_raw_fastqc.log"
            )

    # =====================================================
    # STEP 2 : RAW SEQKIT STATS
    # =====================================================

    def raw_seqkit(self):

        print("\nSTEP 2 : RAW READ STATISTICS")

        raw_reads = glob.glob(
            str(self.input_dir / "*.fastq.gz")
        )

        cmd = f"""
        seqkit stats -a -j {self.threads} {' '.join(raw_reads)} \
        > {self.paths['seqkit']}/raw_read_statistics.tsv
        """

        self.run_command(
            ["bash", "-c", cmd],
            self.paths["logs"] / "raw_seqkit.log"
        )

    # =====================================================
    # STEP 3 : FASTP TRIMMING
    # =====================================================

    def trimming(self):

        print("\nSTEP 3 : FASTP TRIMMING")

        for sample in self.samples:

            r1, r2 = self.reads(sample)

            clean_r1 = self.paths["trimmed"] / f"{sample}_R1.clean.fq.gz"

            clean_r2 = self.paths["trimmed"] / f"{sample}_R2.clean.fq.gz"

            html = self.paths["trimmed"] / f"{sample}.html"

            json = self.paths["trimmed"] / f"{sample}.json"

            cmd = [
                "fastp",

                "-i", str(r1),

                "-I", str(r2),

                "-o", str(clean_r1),

                "-O", str(clean_r2),

                "-w", str(self.threads),

                "-h", str(html),

                "-j", str(json),

                "--detect_adapter_for_pe",

                "-c"
            ]

            self.run_command(
                cmd,
                self.paths["logs"] / f"{sample}_fastp.log"
            )

    # =====================================================
    # STEP 4 : CLEAN FASTQC
    # =====================================================

    def clean_fastqc(self):

        print("\nSTEP 4 : CLEAN FASTQC")

        for sample in self.samples:

            r1 = self.paths["trimmed"] / f"{sample}_R1.clean.fq.gz"

            r2 = self.paths["trimmed"] / f"{sample}_R2.clean.fq.gz"

            cmd = [
                "fastqc",
                "-t", str(self.threads),
                str(r1),
                str(r2),
                "-o",
                str(self.paths["cleanqc"])
            ]

            self.run_command(
                cmd,
                self.paths["logs"] / f"{sample}_clean_fastqc.log"
            )

    # =====================================================
    # STEP 5 : CLEAN SEQKIT STATS
    # =====================================================

    def clean_seqkit(self):

        print("\nSTEP 5 : CLEAN READ STATISTICS")

        clean_reads = glob.glob(
            str(self.paths["trimmed"] / "*.fq.gz")
        )

        cmd = f"""
        seqkit stats -a -j {self.threads} {' '.join(clean_reads)} \
        > {self.paths['seqkit']}/clean_read_statistics.tsv
        """

        self.run_command(
            ["bash", "-c", cmd],
            self.paths["logs"] / "clean_seqkit.log"
        )

    # =====================================================
    # STEP 6 : CHECK BWA INDEX
    # =====================================================

    def bwa_index_check(self):

        print("\nSTEP 6 : REFERENCE INDEX CHECK")

        required = [".amb", ".ann", ".bwt", ".pac", ".sa"]

        missing = []

        for ext in required:

            idx = str(self.reference) + ext

            if not os.path.exists(idx):
                missing.append(idx)

        if len(missing) == 0:

            print("BWA index files already present")

        else:

            print("Missing BWA index files")
            print("Generating BWA index")

            self.run_command(
                ["bwa", "index", str(self.reference)],
                self.paths["logs"] / "bwa_index.log"
            )

    # =====================================================
    # STEP 7 : ALIGNMENT
    # =====================================================

    def alignment(self):

        print("\nSTEP 7 : BWA ALIGNMENT")

        for sample in self.samples:

            r1 = self.paths["trimmed"] / f"{sample}_R1.clean.fq.gz"

            r2 = self.paths["trimmed"] / f"{sample}_R2.clean.fq.gz"

            bam = self.paths["alignment"] / f"{sample}.bam"

            rg = f"@RG\\tID:{sample}\\tSM:{sample}\\tPL:ILLUMINA"

            bwa_cmd = [
                "bwa", "mem",
                "-t", str(self.threads),
                "-M",
                "-R", rg,
                str(self.reference),
                str(r1),
                str(r2)
            ]

            sort_cmd = [
                "samtools", "sort",
                "-@", str(self.threads),
                "-o", str(bam)
            ]

            logfile = self.paths["logs"] / f"{sample}_alignment.log"

            with open(logfile, "w") as log:

                p1 = subprocess.Popen(
                    bwa_cmd,
                    stdout=subprocess.PIPE,
                    stderr=log
                )

                p2 = subprocess.Popen(
                    sort_cmd,
                    stdin=p1.stdout,
                    stdout=log,
                    stderr=log
                )

                p1.stdout.close()

                p2.communicate()

            if p2.returncode != 0:

                raise RuntimeError(
                    if"Alignment failed for {sample}"
                )

             BAM INDEX
            self.run_command(
                [
                    "samtools",
                    "index",
                    str(bam)
                ],
                self.paths["logs"] / f"{sample}_bam_index.log"
            )

    # =====================================================
    # STEP 8 : MARK DUPLICATES
    # =====================================================

    def mark_duplicates(self):

        print("\nSTEP 8 : MARK DUPLICATES")

        for sample in self.samples:

            bam = self.paths["alignment"] / f"{sample}.bam"

            dedup = self.paths["alignment"] / f"{sample}.dedup.bam"

            metrics = self.paths["alignment"] / f"{sample}.metrics.txt"

            cmd = [
                "java",
                "-Xmx8g",
                "-jar",
                self.picard,

                "MarkDuplicates",

                f"I={bam}",

                f"O={dedup}",

                f"M={metrics}",

                "CREATE_INDEX=true"
            ]

            self.run_command(
                cmd,
                self.paths["logs"] / f"{sample}_markdup.log"
            )

    # =====================================================
    # STEP 9 : SAMTOOLS FLAGSTAT
    # =====================================================

    def flagstat_summary(self):

        print("\nSTEP 9 : SAMTOOLS FLAGSTAT")

        summary_file = (
            self.paths["flagstat"] /
            "primary_alignment_summary.tsv"
        )

        with open(summary_file, "w") as summary:

            summary.write(
                "Sample\tPrimary_Alignment_Percent\n"
            )

            for sample in self.samples:

                bam = (
                    self.paths["alignment"] /
                    f"{sample}.dedup.bam"
                )

                flagstat_out = (
                    self.paths["flagstat"] /
                    f"{sample}.flagstat.txt"
                )

                with open(flagstat_out, "w") as fout:

                    process = subprocess.run(
                        [
                            "samtools",
                            "flagstat",
                            str(bam)
                        ],
                        stdout=fout,
                        stderr=subprocess.PIPE,
                        text=True
                    )

                if process.returncode != 0:

                    raise RuntimeError(
                        f"Flagstat failed for {sample}"
                    )

                primary_percent = "NA"

                with open(flagstat_out) as fin:

                    for line in fin:

                        if "primary mapped" in line:

                            start = line.find("(") + 1

                            end = line.find("%")

                            primary_percent = line[start:end]

                            break

                summary.write(
                    f"{sample}\t{primary_percent}\n"
                )


    # =====================================================
    # STEP 10 : MOSDEPTH COVERAGE
    # =====================================================

    def mosdepth_coverage(self):

        print("\nSTEP 10 : MOSDEPTH COVERAGE")

        for sample in self.samples:

            bam = self.paths["alignment"] / f"{sample}.dedup.bam"

            prefix = self.paths["coverage"] / sample

            cmd = [
                "mosdepth",
                "-t", str(self.threads),

                "-n",

                "--by", "1000",

                "--thresholds", "1,5,10,20,30",

                str(prefix),

                str(bam)
            ]

            self.run_command(
                cmd,
                self.paths["logs"] / f"{sample}_mosdepth.log"
            )

    # =====================================================
    # STEP 11 : PREPARE REFERENCE FOR GATK
    # =====================================================

    def prepare_gatk_reference(self):

        print("\nSTEP 11 : GATK REFERENCE PREPARATION")

        fai = str(self.reference) + ".fai"

        if not os.path.exists(fai):

            self.run_command(
                ["samtools", "faidx", str(self.reference)],
                self.paths["logs"] / "samtools_faidx.log"
            )

        dict_file = self.reference.with_suffix(".dict")

        if not os.path.exists(dict_file):

            self.run_command([
                "java",
                "-jar",
                self.picard,
                "CreateSequenceDictionary",
                f"R={self.reference}",
                f"O={dict_file}"

            ], self.paths["logs"] / "picard_dict.log")

    # =====================================================
    # STEP 12 : HAPLOTYPECALLER
    # =====================================================

    def variant_calling(self):

        print("\nSTEP 12 : VARIANT CALLING")

        for sample in self.samples:

            bam = self.paths["alignment"] / f"{sample}.dedup.bam"

            vcf = self.paths["variants"] / f"{sample}.vcf.gz"

            cmd = [
                self.gatk,
                "HaplotypeCaller",

                "-R", str(self.reference),

                "-I", str(bam),

                "-O", str(vcf)
            ]

            self.run_command(
                cmd,
                self.paths["logs"] / f"{sample}_haplotypecaller.log"
            )

    # =====================================================
    # STEP 13 : SPLIT SNPs AND INDELs
    # =====================================================

    def split_variants(self):

        print("\nSTEP 13 : SPLITTING SNPs AND INDELs")

        for sample in self.samples:

            raw_vcf = self.paths["variants"] / f"{sample}.vcf.gz"

            snp_vcf = self.paths["filtered"] / f"{sample}.snps.vcf.gz"

            indel_vcf = self.paths["filtered"] / f"{sample}.indels.vcf.gz"

            self.run_command([
                self.gatk,
                "SelectVariants",

                "-V", str(raw_vcf),

                "--select-type-to-include", "SNP",

                "-O", str(snp_vcf)

            ], self.paths["logs"] / f"{sample}_snps.log")

            self.run_command([
                self.gatk,
                "SelectVariants",

                "-V", str(raw_vcf),

                "--select-type-to-include", "INDEL",

                "-O", str(indel_vcf)

            ], self.paths["logs"] / f"{sample}_indels.log")

    # =====================================================
    # STEP 14 : SNP FILTERING
    # =====================================================

    def snp_filtering(self):

        print("\nSTEP 14 : SNP FILTERING")

        for sample in self.samples:

            input_vcf = self.paths["filtered"] / f"{sample}.snps.vcf.gz"

            filtered_vcf = self.paths["filtered"] / f"{sample}.snps.filtered.vcf.gz"

            pass_vcf = self.paths["filtered"] / f"{sample}.snps.pass.vcf.gz"

            self.run_command([
                self.gatk,
                "VariantFiltration",

                "-V", str(input_vcf),

                "--filter-expression", "QD < 2.0",
                "--filter-name", "LowQD",

                "--filter-expression", "QUAL < 30.0",
                "--filter-name", "LowQUAL",

                "--filter-expression", "FS > 60.0",
                "--filter-name", "HighFS",

                "--filter-expression", "MQ < 40.0",
                "--filter-name", "LowMQ",

                "-O", str(filtered_vcf)

            ], self.paths["logs"] / f"{sample}_snp_filter.log")

            self.run_command([
                self.gatk,
                "SelectVariants",

                "--exclude-filtered",

                "-V", str(filtered_vcf),

                "-O", str(pass_vcf)

            ], self.paths["logs"] / f"{sample}_pass_snps.log")

    # =====================================================
    # STEP 15 : INDEL FILTERING
    # =====================================================

    def indel_filtering(self):

        print("\nSTEP 15 : INDEL FILTERING")

        for sample in self.samples:

            input_vcf = self.paths["filtered"] / f"{sample}.indels.vcf.gz"

            filtered_vcf = self.paths["filtered"] / f"{sample}.indels.filtered.vcf.gz"

            pass_vcf = self.paths["filtered"] / f"{sample}.indels.pass.vcf.gz"

            self.run_command([
                self.gatk,
                "VariantFiltration",

                "-V", str(input_vcf),

                "--filter-expression", "QD < 2.0",
                "--filter-name", "LowQD",

                "--filter-expression", "FS > 200.0",
                "--filter-name", "HighFS",

                "--filter-expression", "QUAL < 30.0",
                "--filter-name", "LowQUAL",

                "-O", str(filtered_vcf)

            ], self.paths["logs"] / f"{sample}_indel_filter.log")

            self.run_command([
                self.gatk,
                "SelectVariants",

                "--exclude-filtered",

                "-V", str(filtered_vcf),

                "-O", str(pass_vcf)

            ], self.paths["logs"] / f"{sample}_pass_indels.log")

    # =====================================================
    # STEP 16 : MERGE PASS VARIANTS
    # =====================================================

    def merge_variants(self):

        print("\nSTEP 16 : MERGING PASS VARIANTS")

        for sample in self.samples:

            snps = self.paths["filtered"] / f"{sample}.snps.pass.vcf.gz"

            indels = self.paths["filtered"] / f"{sample}.indels.pass.vcf.gz"

            merged = self.paths["finalvcf"] / f"{sample}.final.vcf.gz"

            cmd = [
                "java",
                "-jar",
                self.picard,

                "MergeVcfs",

                f"I={snps}",

                f"I={indels}",

                f"O={merged}"
            ]

            self.run_command(
                cmd,
                self.paths["logs"] / f"{sample}_merge.log"
            )

    # =====================================================
    # STEP 17 : MULTIQC
    # =====================================================

    def multiqc_report(self):

        print("\nSTEP 17 : MULTIQC REPORT")

        cmd = [
            "multiqc",
            str(self.output_dir),
            "-o",
            str(self.paths["multiqc"])
        ]

        self.run_command(
            cmd,
            self.paths["logs"] / "multiqc.log"
        )

    # =====================================================
    # RUN COMPLETE PIPELINE
    # =====================================================

    def run_pipeline(self):

        print("\nPIPELINE STARTED :", datetime.datetime.now())

        self.raw_fastqc()

        self.raw_seqkit()

        self.trimming()

        self.clean_fastqc()

        self.clean_seqkit()

        self.bwa_index_check()

        self.alignment()

        self.mark_duplicates()

        self.flagstat_summary()

        self.mosdepth_coverage()

        self.prepare_gatk_reference()

        self.variant_calling()

        self.split_variants()

        self.snp_filtering()

        self.indel_filtering()

        self.merge_variants()

        self.multiqc_report()

        print("\nPIPELINE COMPLETED :", datetime.datetime.now())


# =========================================================
# MAIN
# =========================================================

if __name__ == "__main__":

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "-i",
        "--input",
        required=True,
        help="Directory containing FASTQ files"
    )

    parser.add_argument(
        "-r",
        "--reference",
        required=True,
        help="Reference genome prefix WITHOUT extension"
    )

    parser.add_argument(
        "-o",
        "--output",
        required=True,
        help="Output directory"
    )

    parser.add_argument(
        "-t",
        "--threads",
        type=int,
        default=8
    )

    args = parser.parse_args()

    pipeline = WGSPipeline(
        args.input,
        args.reference,
        args.output,
        args.threads
    )

    pipeline.run_pipeline()
