#!/bin/bash

# define base directory
BASE_DIR="/root/Alix"
# Define reference genome
REFERENCE="GCF_000011805.1_ASM1180v1_genomic.fna"
# Directory containing de novo assemblies
ASSEMBLIES_DIR="${BASE_DIR}/output/asm/"
# Output directory for SAM/BAM files
OUTPUT_DIR="${BASE_DIR}/output/msa/"

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# define samples
samples=(
    "1561"
    "1561ER"
    "1561TB"
    "1562"
    "1563"
    "1563A"
    "1563B"
    "1563C"
    "1563D"
    "1564"
)

# Loop through each sample and align to the reference
for sample in "${samples[@]}"; do
    # Define the input and output files
    ASSEMBLY="${ASSEMBLIES_DIR}/${sample}/filtered_scaffolds.fasta"
    OUTPUT_SAM="${OUTPUT_DIR}/${sample}.sam"

    # Align assembly to reference using Minimap2
    minimap2 -x asm5 -a "$REFERENCE" "$ASSEMBLY" > "$OUTPUT_SAM"

    # Optional: Convert SAM to BAM and sort
    samtools view -bS "$OUTPUT_SAM" | samtools sort -o "${OUTPUT_DIR}/${sample}.sorted.bam"
    samtools index "${OUTPUT_DIR}/${sample}.sorted.bam"

    # Clean up intermediate SAM file 
    rm "$OUTPUT_SAM"

    # extract aligned sequence as fasta
    samtools fasta "${OUTPUT_DIR}/${sample}.sorted.bam" > "${OUTPUT_DIR}/${sample}.sorted.fasta"
done

# Combine all FASTA files into a single file
cat $OUTPUT_DIR*.sorted.fasta > "${OUTPUT_DIR}/combined_sequences.fasta"