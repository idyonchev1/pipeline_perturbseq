Snakemake pipeline for processing Perturb-seq data. 
Designed to run multiqc, cellranger count and cellranger aggr to aggregate output from multiple wells.

To execute on SLURM cluster: snakemake --cluster "sbatch -t {resources.time} --mem={resources.mem_mb} -c {threads}" -j 30 --use-conda 

To build rulegraph:snakemake --rulegraph | dot -Tpdf > rulegraph.pdf

Requirements:

fastq_screen_index : folder with bowtie2 indices for fastq screen. Example indices in config/fastq_screen.conf
samples.txt : column Sample
libraries_{Sample}.txt : columns 'fastqs,sample,library_type'
aggr.csv: columns 'sample_id,molecule_h5,experiment'
feature_ref.csv : columns 'id,name,read,pattern,sequence,feature_type,target_gene_id,target_gene_name'
refdata-gex-GRCh38-2024-A : Genome+Transcriptome to align to.
chrNameLength.txt - chromosome sizes

Tools:
Snakemake version: 7.28.2
CellRanger version: 8.0.1
MultiQC v1.21
FastQC 0.12.1

