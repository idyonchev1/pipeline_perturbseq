rule cellranger:
    output: "{sample}/outs/web_summary.html"
    log: "log/{sample}/cellranger_count.log"
    resources:
        mem_mb=128000,
        disk_mb=200000,
        time="72:00:00"
    threads: 16
    run:
        cmd1 ="rm -r {wildcards.sample}"
        cmd2 = "cellranger count --id={wildcards.sample} --libraries=libraries_{wildcards.sample}.csv --transcriptome=refdata-gex-GRCh38-2024-A-dCas9Zim3 --feature-ref=feature_ref.csv --create-bam=false --include-introns=false --localcores=16 --localmem=115 --jobmode=local &> {log}"
        shell(cmd1)
        shell(cmd2)
