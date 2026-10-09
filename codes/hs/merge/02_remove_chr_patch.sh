# 1. Before removal
echo "=== Before removal ==="
awk -F'\t' '{print $3}' master_lincRNA.gtf | sort | uniq -c
echo "Total lines: $(wc -l < master_lincRNA.gtf)"

# 2. Check what chr patches exist
echo "=== Chr patches ==="
awk -F'\t' '{print $1}' master_lincRNA.gtf | sort -u | grep -v "^chr[0-9XYM]*$" > chr_patch.txt

# 3. Remove chr patches
awk -F'\t' '$1 ~ /^chr([0-9]+|X|Y|M)$/' master_lincRNA.gtf > master_lincRNA_chr_patch_removed.gtf

# 4. After removal
echo "=== After removal ==="
awk -F'\t' '{print $3}' master_lincRNA_1.gtf | sort | uniq -c
echo "Total lines: $(wc -l < master_lincRNA_1.gtf)"
