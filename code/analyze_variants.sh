#transfer vcf files into one location 
python ./code/copy_gd.py -e ./drive-download-20241111T172643Z-001/QUO1002548-20221114/Variant\ Calling/ -o ./output/variants/vcf/ -t vcf

cd ./output/variants/vcf/

# convert vcf files to zipped bcf file
bcftools view -O b  -o VF_WT_1561.bcf.gz VF_WT_1561.vcf
bcftools view -O b  -o VF_WT_1561ER.bcf.gz VF_WT_1561ER.vcf
bcftools view -O b  -o VF_WT_1561TB.bcf.gz VF_WT_1561TB.vcf
bcftools view -O b  -o VF_WT_1562.bcf.gz VF_WT_1562.vcf
bcftools view -O b  -o VF_WT_1563.bcf.gz VF_WT_1563.vcf
bcftools view -O b  -o VF_WT_1563A.bcf.gz VF_WT_1563A.vcf
bcftools view -O b  -o VF_WT_1563B.bcf.gz VF_WT_1563B.vcf
bcftools view -O b  -o VF_WT_1563C.bcf.gz VF_WT_1563C.vcf
bcftools view -O b  -o VF_WT_1563D.bcf.gz VF_WT_1563D.vcf
bcftools view -O b  -o VF_WT_1564.bcf.gz VF_WT_1564.vcf

# generate index for each sample
tabix  VF_WT_1561.bcf.gz 
tabix  VF_WT_1561ER.bcf.gz 
tabix  VF_WT_1561TB.bcf.gz 
tabix  VF_WT_1562.bcf.gz 
tabix  VF_WT_1563.bcf.gz 
tabix  VF_WT_1563A.bcf.gz 
tabix  VF_WT_1563B.bcf.gz 
tabix  VF_WT_1563C.bcf.gz 
tabix  VF_WT_1563D.bcf.gz 
tabix  VF_WT_1564.bcf.gz 

# merge each group into a single vcf
bcftools merge -o ../merged_ACS.vcf  -O v  VF_WT_1561.bcf.gz VF_WT_1562.bcf.gz VF_WT_1563.bcf.gz VF_WT_1563A.bcf.gz 
bcftools merge -o ../merged_WT.vcf  -O v  VF_WT_1561ER.bcf.gz VF_WT_1561TB.bcf.gz VF_WT_1563D.bcf.gz VF_WT_1564.bcf.gz 

cd ..
# normalize
bcftools norm -m -both -o ACS_normalized.vcf -O v merged_ACS.vcf 
bcftools norm -m -both -o WT_normalized.vcf -O v merged_WT.vcf 
# sort
bcftools sort -o ACS_sorted.vcf ACS_normalized.vcf
bcftools sort -o WT_sorted.vcf WT_normalized.vcf
# covert to bcf.gz
bcftools view -O b  -o ACS_sorted.bcf.gz ACS_sorted.vcf
bcftools view -O b  -o WT_sorted.bcf.gz WT_sorted.vcf
# make index
tabix ACS_sorted.bcf.gz
tabix WT_sorted.bcf.gz

# mutations exclusive to ACS
bcftools isec -C -o exclusive_to_ACS.vcf -O v ACS_sorted.bcf.gz WT_sorted.bcf.gz
# mutations exclusive to WT
bcftools isec -C -o exclusive_to_WT.vcf -O v WT_sorted.bcf.gz ACS_sorted.bcf.gz

