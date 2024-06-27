#To generate bigwig files of every pertubation and export a Non-Targeting control bam for target enrichment probe design
#run this script after generating a "cell-target-identities" text file containing two columns of cell identities and sgRNA target. If you analyzed data on an aggr object,
#you can rename each GEM identitiy to end in -1 and export separate text files to be run on each GEM, prior to merging back together. 
#RPKM normalized tracks output for IGV
#For obtaining accurate 5' start sites - remove read 2 and use the first base of read 1. --samFlagExclude 128 --Offset 1 in export_cpm_bigwig

import os
import pandas as pd

samples_table = pd.read_table("samples.txt")
samples=list(samples_table.Samples.unique())

targets_table = pd.read_table("targets.txt")
target=list(targets_table.Target.unique())

rule all:
    input:[expand("logs/{samples}_sinto_done.txt",samples=samples),expand("demultiplex_cells/{target}_deduplicated_merged.bam",target=target),expand("demultiplex_cells/tracks/{target}_fwd.bw",target=target),expand("demultiplex_cells/tracks/{target}_rev.bw",target=target),"probe_design/fastq/Non-Targeting.1.fastq.gz","probe_design/fastq/Non-Targeting.2.fastq.gz","probe_design/Non-Targeting.bedgraph"]

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
    input:
        "logs/{sample}_sinto_done.txt"
    output:
        temp(expand("demultiplex_cells/{sample}/{target}.bam",sample="{sample}",target=target))
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
    output:
        "demultiplex_cells/{target}_deduplicated_merged.bam"
    threads:4
    resources:
        mem_mb=10000,
        time="24:00:00"
    shell:
        """
        samtools merge -@{threads} -o {output} {input}
        samtools index -M {output}
        """

rule extract_control_reads:
    input:
        "demultiplex_cells/Non-Targeting_deduplicated_merged.bam"
    threads:4
    resources:
        mem_mb=10000,
        time="24:00:00"
    log:
        "logs/extract_control_reads.log"
    output:
        read1_temp=temp("probe_design/fastq/Non-Targeting-untrimmed.1.fastq"),
        read2_un=temp("probe_design/fastq/Non-Targeting.2.fastq"),
        read1_un=temp("probe_design/fastq/Non-Targeting.1.fastq"),
        read1="probe_design/fastq/Non-Targeting.1.fastq.gz",
        read2="probe_design/fastq/Non-Targeting.2.fastq.gz"
    conda:
        "cgat-apps"
    shell:
        """
        samtools collate -u -O {input} -@4 | samtools fastq -f 3 -1 {output.read1_temp} -2 {output.read2_un} -@4 &> {log}
        cutadapt -u 13 {output.read1_temp} -o {output.read1_un} &> {log}
        gzip < {output.read1_un} > {output.read1}
        gzip < {output.read2_un} > {output.read2}
        """

rule map_control_reads_transcriptome:
    input:
        sample=["probe_design/fastq/Non-Targeting.1.fastq.gz", "probe_design/fastq/Non-Targeting.2.fastq.gz"]
    output:
        out1=temp("probe_design/Non-Targeting_unsorted.bam"),
        out2="probe_design/Non-Targeting_sorted_deduplicated.bam"
    log:
        "logs/map_control_reads_transcriptome.log",
    threads:8
    resources:
        mem_mb=100000,
        time="3-00:00:00"
    shell:
        """
        bowtie2 -x bowtie2_transcriptome_index/refdata-gex-GRCh38-2024-A-dCas9Zim3 -1 {input[0]} -2 {input[1]} -X 600 --very-sensitive-local -p 8 --no-mixed --no-discordant -k 100 2>>{log} | samtools view -f 3 -bS - > {output.out1} 
        samtools sort {output.out1} -@ 8 -o {output.out2}
        samtools index {output.out2}        
        """

rule export_bedgraph:
    input:"probe_design/{sample}_sorted_deduplicated.bam"
    output:"probe_design/{sample}.bedgraph"
    threads: 1
    resources:
        mem_mb=8000,
        time="24:00:00"
    shell:
        "genomeCoverageBed -bga -ibam {input} | sort -k1,1 -k2,2n > {output}"

rule export_stranded_bigwigs:
    input:
        "demultiplex_cells/{target}_deduplicated_merged.bam"
    output:
        fwdbam1=temp("demultiplex_cells/tracks/{target}_fwd1.bam"),
        fwdbam2=temp("demultiplex_cells/tracks/{target}_fwd2.bam"),
        fwdbam=temp("demultiplex_cells/tracks/{target}_fwd.bam"),
        revbam1=temp("demultiplex_cells/tracks/{target}_rev1.bam"),
        revbam2=temp("demultiplex_cells/tracks/{target}_rev2.bam"),
        revbam=temp("demultiplex_cells/tracks/{target}_rev.bam"),
        fwdbw="demultiplex_cells/tracks/{target}_fwd.bw",
        revbw="demultiplex_cells/tracks/{target}_rev.bw"
    log:
        "logs/tracks/deeptools_{target}.log"
    threads: 4
    resources:
        mem_mb=10000,
        time="24:00:00"
    conda:
        "cgat-apps"
    shell:
        """
        samtools view -b -f 128 -F 16 {input} > {output.revbam1}
        samtools view -b -f 80 {input} > {output.revbam2}
        samtools merge -f {output.revbam} {output.revbam1} {output.revbam2}
        samtools index {output.revbam}
        samtools view -b -f 144 {input} > {output.fwdbam1}
        samtools view -b -f 64 -F 16 {input} > {output.fwdbam2}
        samtools merge -f {output.fwdbam} {output.fwdbam1} {output.fwdbam2}
        samtools index {output.fwdbam}
        bamCoverage --bam {output.revbam} -o {output.revbw} -p 4 --normalizeUsing RPKM --exactScaling --binSize 10 --effectiveGenomeSize 2913022398 &> {log} 
        bamCoverage --bam {output.fwdbam} -o {output.fwdbw} -p 4 --normalizeUsing RPKM --exactScaling --binSize 10 --effectiveGenomeSize 2913022398 &> {log} 
        rm demultiplex_cells/tracks/{wildcards.target}*.bai
        """
