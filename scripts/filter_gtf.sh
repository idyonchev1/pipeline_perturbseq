#!/bin/bash

if [ "$#" -ne 3 ]; then
    echo "Usage: $0 input.gtf enst_list.txt output.gtf"
    exit 1
fi

INPUT_GTF=$1
ENST_LIST=$2
OUTPUT_GTF=$3

if [ ! -f "$INPUT_GTF" ]; then
    echo "Error: Input GTF file $INPUT_GTF does not exist."
    exit 1
fi

if [ ! -f "$ENST_LIST" ]; then
    echo "Error: ENST list file $ENST_LIST does not exist."
    exit 1
fi

# Create a temporary file with ENST IDs in a grep-friendly format
TMP_ENST=$(mktemp)
sed 's/^/transcript_id "/; s/$/"/' "$ENST_LIST" > "$TMP_ENST"

# Filter the GTF file
grep -F -f "$TMP_ENST" "$INPUT_GTF" > "$OUTPUT_GTF"

# Clean up
rm "$TMP_ENST"

# Count entries
INPUT_COUNT=$(wc -l < "$INPUT_GTF")
OUTPUT_COUNT=$(wc -l < "$OUTPUT_GTF")
RETAINED_COUNT=$OUTPUT_COUNT
REMOVED_COUNT=$((INPUT_COUNT - OUTPUT_COUNT))

echo "Processing complete. Output written to $OUTPUT_GTF"
echo "Total entries in input file: $INPUT_COUNT"
echo "Entries retained: $RETAINED_COUNT"
echo "Entries removed: $REMOVED_COUNT"
