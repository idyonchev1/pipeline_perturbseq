FASTQ, = glob_wildcards("input/{fastq}_R1_001.fastq.gz")

rule fastqc:
    input:
        ["input/{fastq}_R1_001.fastq.gz", "input/{fastq}_R2_001.fastq.gz"]
    output:
        ["data/fastqc/{fastq}_R1_001_fastqc.html","data/fastqc/{fastq}_R2_001_fastqc.html"]
    params:
    threads: 8
    resources:
        mem_mb=16000,
        time="48:00:00"
    shell:
        "fastqc -o data/fastqc/ -t 1 {input}"

rule fastq_screen:
    input:
        ["input/{fastq}_R1_001.fastq.gz","input/{fastq}_R2_001.fastq.gz"]
    output:
        txt="data/fastqc/{fastq}.fastq_screen.txt",
        png="data/fastqc/{fastq}.fastq_screen.png"
    params:
        fastq_screen_config="config/fastq_screen.conf",
        subset=100000,
        aligner='bowtie2'
    threads: 8
    resources:
        mem_mb=16000,
        time="48:00:00"
    conda:
        "wrapper"
    wrapper:
        "v2.0.0/bio/fastq_screen"

rule multiqc:
    input:
        expand("data/fastqc/{fastq}_R1_001_fastqc.html", fastq=FASTQ),
        expand("data/fastqc/{fastq}_R2_001_fastqc.html",fastq=FASTQ),
        expand("data/fastqc/{fastq}.fastq_screen.txt",fastq=FASTQ)
    output:
        "data/multiqc/multiqc_report.html"
    threads: 1
    resources:
        mem_mb=4000,
        time="24:00:00"
    shell:
        "multiqc data/fastqc -o data/multiqc"
