#total transcript and genes in master_lincRNA_hg38.gtf
awk -F'\t' '$3=="transcript"{nt++} $3=="gene"{ng++} END{print ng+0" genes, "nt+0" transcripts"}' master_lincRNA_hg38.gtf

#number of transcript in master_lincRNA_chr_patch_removed.gtf (this gives the number of redundant transcripts prior applying exclusion criteria)

awk -F'\t' '$3=="transcript"{n++} END{print n+0}' master_lincRNA_chr_patch_removed.gtf

#intergenic lncRNAs (transcripts) map to genes originate from NONCODE v6 and RNA Central v26

awk -F'\t' '$2!~/^(HAVANA|ENSEMBL)$/{if($3=="gene")g++; if($3=="transcript")t++} END{print g+0" genes, "t+0" transcripts"}' master_lincRNA_hg38.gtf


