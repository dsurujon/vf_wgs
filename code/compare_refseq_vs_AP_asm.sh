wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/011/805/GCF_000011805.1_ASM1180v1/GCF_000011805.1_ASM1180v1_genomic.fna.gz
gunzip GCF_000011805.1_ASM1180v1_genomic.fna.gz

fastANI -k 20 -q ./drive-download-20241111T172643Z-001/ES114_annot_AP.fasta \
    -r ./GCF_000011805.1_ASM1180v1_genomic.fna \
    -o ./output/fastani/ASM1180v1_vs_newAP_fastANI.tsv \
    --visualize
    
poetry run python ./code/fastani_visualize.py \
    ./drive-download-20241111T172643Z-001/ES114_annot_AP.fasta \
    ./GCF_000011805.1_ASM1180v1_genomic.fna \
    ./output/fastani/ASM1180v1_vs_newAP_fastANI.tsv.visual \
    ./output/fastani/ASM1180v1_vs_newAP_fastANI.png