#check unique entries in three gtf
for gtf in \
    gencode_curated_lincRNA.gtf \
    noncode_curated_lincRNA_16_1.gtf \
    rnacentral_curated_lincRNA_v2_9_1.gtf
do
    echo "===== $gtf ====="
    grep -v "^#" "$gtf" | cut -f3 | sort | uniq -c
done