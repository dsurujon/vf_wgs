#!/bin/bash

# Array of sample identifiers
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

for sample in "${samples[@]}"; do
    spades.py -o output/asm/${sample}/ \
        -1 output/fastq_dedup/${sample}_R1_dedup.fastq.gz \
        -2 output/fastq_dedup/${sample}_R2_dedup.fastq.gz \
        --isolate
    seqtk seq -L 500 output/asm/${sample}/scaffolds.fasta > output/asm/${sample}/filtered_scaffolds.fasta
    poetry run quast.py output/asm/${sample}/filtered_scaffolds.fasta -o output/asm/${sample}/quast

    poetry run bakta output/asm/${sample}/filtered_scaffolds.fasta --output output/asm/${sample}/bakta/ --force --skip-trna --skip-crispr
done