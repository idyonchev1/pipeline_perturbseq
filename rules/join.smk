import glob
import os
from pathlib import Path
import re

# Function to get the base pattern without S number
def get_base_pattern(filename):
    # Extract the common pattern and R number
    pattern = re.match(r'(PS0004_17111FL-65-01-\d+)_S\d+_(L007)_(R[12])_(001\.fastq\.gz)', filename)
    if pattern:
        return pattern.group(1), pattern.group(2), pattern.group(3), pattern.group(4)
    return None

# Get all files
input_dir = "input/raw"
ALL_FILES = sorted(glob.glob(f"{input_dir}/*.fastq.gz"))

# Group files by their base pattern
file_groups = {}
for filepath in ALL_FILES:
    filename = os.path.basename(filepath)
    pattern = get_base_pattern(filename)
    if pattern:
        base, lane, read, suffix = pattern
        key = (base, lane, read, suffix)
        if key not in file_groups:
            file_groups[key] = []
        file_groups[key].append(filepath)

# Create list of expected output files
OUTPUT_FILES = []
for (base, lane, read, suffix) in file_groups:
    # Get S numbers for this group
    s_numbers = sorted([re.search('_S(\d+)_', os.path.basename(f)).group(1) 
                       for f in file_groups[(base, lane, read, suffix)]])
    
    # Create new filename with mRNA-GEM format
    new_base = base.replace("17111FL-65-01-0", "mRNA-GEM")
    output_name = f"{new_base}_S{'-S'.join(s_numbers)}_{lane}_{read}_{suffix}"
    OUTPUT_FILES.append(output_name)

wildcard_constraints:
    prefix="PS0004_mRNA-GEM\d+",
    suffix="001\.fastq\.gz"

rule all:
    input:
        expand("input/{output}", output=[os.path.basename(f) for f in OUTPUT_FILES])

rule combine_and_rename_fastq:
    input:
        lambda wildcards: file_groups[(
            wildcards.prefix.replace("mRNA-GEM", "PS0004_17111FL-65-01-0"),
            wildcards.lane,
            wildcards.read,
            wildcards.suffix
        )]
    output:
        "input/{prefix}_S{s_numbers}_{lane}_{read}_{suffix}"
    threads: 1
    resources:
        mem_mb=4000,
        time="24:00:00"
    log:
        "logs/combine_fastq/{prefix}_S{s_numbers}_{lane}_{read}_{suffix}.log"
    shell:
        """
        cat {input} > {output} 2> {log}
        """

# Print some debug information
print("\nFile groups:")
for key, files in file_groups.items():
    print(f"Group: {key}")
    for f in files:
        print(f"  {f}")
print("\nExpected outputs:")
for f in OUTPUT_FILES:
    print(f"  {f}")
