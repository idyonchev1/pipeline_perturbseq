#To generate bigwig files of every pertubation and export a Non-Targeting control bam for target enrichment probe design
#run this script after generating a "cell-target-identities" text file containing two columns of cell identities and sgRNA target. If you analyzed data on an aggr object,
#you can rename each GEM identitiy (e.g. -2) to -1 and export separate text files. The bams will be merged down the line. 

import os
import pandas as pd

samples_table = pd.read_table("samples.txt")
samples=list(samples_table.Samples.unique())

num_targets = len(target)

rule all:
    input:[expand("logs/{samples}_sinto_done.txt",samples=samples),expand("logs/{samples}_index_done.txt",samples=samples)]

rule extract_targets_sinto:
    input:
        sample="{samples}/outs/possorted_genome_bam.bam",
        targets="{samples}_target_identities.tsv"
    output:touch("logs/{samples}_sinto_done.txt")
    log:
        "logs/{samples}_extract_targets_sinto.log"
    threads:32
    resources:
        mem_mb=100000,
        time="2-00:00:00"
    conda:
        "sinto"
    shell:"sinto filterbarcodes -b {input.sample} -c {input.targets} --outdir demultiplex_cells/{wildcards.samples} -p 32 &> {log}"

rule index:
    input:"logs/{samples}_sinto_done.txt"
    output:touch("logs/{samples}_index_done.txt")
    threads:1
    resources:
        mem_mb=4000,
        time="24:00:00"
    shell:
        """
        for bamfile in demultiplex_cells/{wildcards.samples}/*.bam; do
            samtools index $bamfile
        done
        """
