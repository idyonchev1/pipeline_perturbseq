#5' scRNA-seq sometimes shows spurious 3' peaks where the polymerase fails to elongate, which in the case of very long 3' UTRs can dominate the transcript counts. To be confident in the identity of the transcript we are quantifying, we want to target the 5' ends of the transcript. To prevent the alghoritm selecting final exons and UTRs for probe design, this script removes the last exon in transcripts that are >2 exons long.
#!/bin/bash

if [ "$#" -ne 2 ]; then
    echo "Usage: $0 input.gtf output.gtf"
    exit 1
fi

INPUT_FILE=$1
OUTPUT_FILE=$2

if [ ! -f "$INPUT_FILE" ]; then
    echo "Error: Input file $INPUT_FILE does not exist."
    exit 1
fi

awk '
BEGIN {FS="\t"; OFS="\t"}
{
    if (tolower($3) ~ /exon/) {
        print "Processing exon line: " $0 > "/dev/stderr"
        if (match($0, /transcript_id "?([^";]+)"?/, tid) && match($0, /exon_number "?([0-9]+)"?/, enum)) {
            transcript = tid[1]
            exon_num = enum[1]
            
            if (!(transcript in max_exon) || exon_num > max_exon[transcript]) {
                max_exon[transcript] = exon_num
            }
            
            exons[transcript][exon_num] = $0
            exon_count++
        } else {
            print "Warning: Could not extract transcript_id or exon_number from line: " $0 > "/dev/stderr"
        }
    } else {
        print $0
    }
}
END {
    print "Total exons processed: " exon_count > "/dev/stderr"
    for (transcript in exons) {
        print "Processing transcript: " transcript ", max exon: " max_exon[transcript] > "/dev/stderr"
        for (exon_num in exons[transcript]) {
            if (max_exon[transcript] <= 2 || (max_exon[transcript] > 2 && exon_num != max_exon[transcript])) {
                print exons[transcript][exon_num]
                kept_exons++
            }
        }
    }
    print "Exons kept: " kept_exons > "/dev/stderr"
}' "$INPUT_FILE" > "$OUTPUT_FILE"

echo "Processing complete. Output written to $OUTPUT_FILE"

INPUT_EXON_COUNT=$(grep -iP '\texon\t' "$INPUT_FILE" | wc -l)
OUTPUT_EXON_COUNT=$(grep -iP '\texon\t' "$OUTPUT_FILE" | wc -l)

echo "Number of exons in input file: $INPUT_EXON_COUNT"
echo "Number of exons in output file: $OUTPUT_EXON_COUNT"
echo "Exons removed: $((INPUT_EXON_COUNT - OUTPUT_EXON_COUNT))"
