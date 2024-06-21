#To generate bigwig files of every pertubation and export a Non-Targeting control bam for target enrichment probe design
#run this script after generating a "cell-target-identities" text file containing two columns of cell identities and sgRNA target. If you analyzed data on an aggr object,
#you can rename each GEM identitiy to end in -1 and export separate text files to be run on each GEM, prior to merging back together. 

import os
import pandas as pd

samples_table = pd.read_table("samples.txt")
samples=list(samples_table.Samples.unique())

targets_table = pd.read_table("targets.txt")
target=list(targets_table.Target.unique())

rule all:
    input:[expand("logs/{samples}_sinto_done.txt",samples=samples),expand("demultiplex_cells/{target}_deduplicated_merged.bam",target=target)]

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

rule checkpoint:
    input:"logs/{sample}_sinto_done.txt"
    output:temp(expand("demultiplex_cells/{sample}/{target}.bam",sample="{sample}",target=target))
    threads:1
    resources:
        mem_mb=2000,
        time="1:00:00"
    shell:"touch {input}"

rule deduplicate:
    input:
        dummy="logs/{samples}_sinto_done.txt",
        sample="demultiplex_cells/{samples}/{target}.bam"
    output:temp("demultiplex_cells/{samples}/{target}_deduplicated.bam")
    log:
        "logs/deduplicate/{samples}_{target}_deduplicate.log"
    threads:4
    resources:
        mem_mb=10000,
        time="2-00:00:00"
    conda:
        "cgat-apps"
    shell:
        """
        umi_tools dedup --umi-tag=UR --cell-tag=CR --extract-umi-method=tag --per-cell --paired --chimeric-pairs=discard --unpaired-reads=discard --stdin={input.sample} --stdout={output} --mapping-quality=255 &> {log}
        samtools index -M {output}
        """

rule merge:
    input:
        lambda wildcards: expand("demultiplex_cells/{samples}/{target}_deduplicated.bam", samples=samples, target=wildcards.target)
    output:"demultiplex_cells/{target}_deduplicated_merged.bam"
    threads:4
    resources:
        mem_mb=10000,
        time="24:00:00"
    shell:
        """
        samtools merge -@{threads} -o {output} {input}
        samtools index -M {output}
        """
