#!/bin/bash

# Create a folder named "Reads" if it doesn't exist
mkdir -p Reads

# Array of SRR accession numbers
SRR_LIST=(
)

# Loop through each accession number
for SRR in "${SRR_LIST[@]}"; do
    echo "Downloading and extracting reads for $SRR..."

    # Extract 1,000,000 reads from the SRA file
    fastq-dump --split-files -X 100000 $SRR -O Reads
done

echo "Download and extraction complete."
