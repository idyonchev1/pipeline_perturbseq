import glob
import os
from pathlib import Path

# Get list of R1 files and derive R2 files
R1_FILES = sorted([f for f in glob.glob("*.fastq.gz") if "_R1_" in f])
R2_FILES = [f.replace("_R1_", "_R2_") for f in R1_FILES]

# Verify all R2 files exist
for r1, r2 in zip(R1_FILES, R2_FILES):
    if not os.path.exists(r2):
        raise ValueError(f"Missing R2 file for {r1}: {r2}")

rule all:
    input:
        expand("half_reads/{r1}", r1=[os.path.basename(f) for f in R1_FILES]),
        expand("half_reads/{r2}", r2=[os.path.basename(f) for f in R2_FILES])

rule subsample_paired_fastq:
    input:
        r1="{prefix}_R1_{suffix}",
        r2="{prefix}_R2_{suffix}"
    output:
        r1="half_reads/{prefix}_R1_{suffix}",
        r2="half_reads/{prefix}_R2_{suffix}"
    params:
        seed=42,
        fraction=0.5
    conda:
        "subsample"
    threads: 4
    resources:
        mem_mb=8000,
        time="24:00:00"
    shell:
        """
        # Use same seed for both files to maintain pairs
        seqtk sample -s {params.seed} {input.r1} {params.fraction} | gzip > {output.r1}
        seqtk sample -s {params.seed} {input.r2} {params.fraction} | gzip > {output.r2}
        """
