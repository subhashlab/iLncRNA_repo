#strip versions following gene_name to make them non-versioned gene_names(GENCODE standard)

sed -E 's/(gene_name "[A-Za-z0-9_.-]+)\.[0-9]+"/\1"/g' master_lincRNA_chr_patch_removed.gtf > master_lincRNA_nonversioned.gtf