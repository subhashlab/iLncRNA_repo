#verify common entries between noncode and gencode

bedtools intersect \
-a noncode_curated_lincRNA_16_1_gene_level.bed \
-b gencode_curated_lincRNA_gene_level.bed \
-wa -wb | head

bedtools intersect \
-a noncode_curated_lincRNA_16_1_gene_level.bed \
-b gencode_curated_lincRNA_gene_level.bed \
-u | wc -l

#verify common entries between rna central and gencode 
bedtools intersect \
-a rnacentral_curated_lincRNA_v2_9_1_gene_level.bed \
-b gencode_curated_lincRNA_gene_level.bed \
-wa -wb | head

bedtools intersect \
-a rnacentral_curated_lincRNA_v2_9_1_gene_level.bed \
-b gencode_curated_lincRNA_gene_level.bed \
-u | wc -l

# verify common entries between noncode and rna central 

bedtools intersect \
-a rnacentral_curated_lincRNA_v2_9_1_gene_level.bed \
-b noncode_curated_lincRNA_16_1_gene_level.bed \
-wa -wb | head

bedtools intersect \
-a rnacentral_curated_lincRNA_v2_9_1_gene_level.bed \
-b noncode_curated_lincRNA_16_1_gene_level.bed \
-u | wc -l