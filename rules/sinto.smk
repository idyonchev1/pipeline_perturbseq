#To generate bigwig files of every pertubation and export a Non-Targeting control bam for target enrichment probe design
#run this script after generating a "cell-target-identities" text file containing two columns of cell identities and sgRNA target. If you analyzed data on an aggr object,
#you can rename each GEM identitiy (e.g. -2) to -1 and export separate text files. The bams will be merged down the line. 

import os
import pandas as pd

samples_table = pd.read_table("samples.txt")
samples=list(samples_table.Samples.unique())

rule all:
    input:[expand("logs/{samples}_sinto_done.txt",samples=samples)]

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
    shell:
        """
        sinto filterbarcodes -b {input.sample} -c {input.targets} --outdir demultiplex_cells/{wildcards.samples} -p 32 &> {log}
        rm demultiplex_cells/{wildcards.samples}/*_*
        samtools index -M demultiplex_cells/{wildcards.samples}/*bam
        """
