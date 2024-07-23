#Script to use UMI-tools-deduplicated 5' scRNA-seq bam file for target enrichment probe design. A gtf containing transcripts of interest is passed to featurecounts. The total read counts per gene and per transcript are separately calculated, after which transcripts containing >80% of the read counts are selected for probe design.
#5' scRNA-seq is forward stranded, change accordingly for 3'
#!/bin/bash

# Check if the correct number of arguments is provided
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <gtf_file> <bam_file> <output_prefix>"
    exit 1
fi

# Assign input arguments to variables
gtf_file="$1"
bam_file="$2"
output_prefix="$3"

# Run featureCounts for gene-level counts
featureCounts -p -s 1 -t exon -g gene_id -M -O -a "$gtf_file" -o "${output_prefix}_gene.txt" "$bam_file"

# Run featureCounts for transcript-level counts
featureCounts -p -s 1 -t exon -g transcript_id -M -O -a "$gtf_file" -o "${output_prefix}_transcript.txt" "$bam_file"

echo "FeatureCounts analysis complete. Output files:"
echo "${output_prefix}_gene.txt"
echo "${output_prefix}_transcript.txt"
