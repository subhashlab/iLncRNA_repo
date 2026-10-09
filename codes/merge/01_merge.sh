#worflow same as merge_v2.sh, using gencode_curated_v2.gtf now
#STEP1: 
cat gencode_curated_lincRNA_v2.gtf \
    noncode_curated_lincRNA_16_1.gtf \
    rnacentral_curated_lincRNA_12.gtf | \
grep -v "^#" > merged_temp.gtf
echo "Step 1 done: $(wc -l < merged_temp.gtf) lines"

#STEP2: 
awk -F'\t' '$3=="gene"' merged_temp.gtf | \
sort -k1,1V -k4,4n -k5,5n | \
awk -F'\t' '{match($9, /gene_id "([^"]+)"/, a); print a[1]}' > gene_order.txt
echo "Step 2 done: $(wc -l < gene_order.txt) genes"

#STEP3:
awk -F'\t' '
NR==FNR {
    rank[$1] = NR
    next
}
{
    match($9, /gene_id "([^"]+)"/, a); gid = a[1]
    match($9, /transcript_id "([^"]+)"/, b); tid = b[1]
    if ($3=="gene") {
        gene[gid] = $0
    } else if ($3=="transcript") {
        trans[gid SUBSEP tid] = $0
        t_order[gid] = (t_order[gid] == "" ? tid : t_order[gid] "\t" tid)
    } else if ($3=="exon") {
        exon[gid SUBSEP tid] = exon[gid SUBSEP tid] $0 "\n"
    }
}
END {
    while ((getline gid < "gene_order.txt") > 0) {
        print gene[gid]
        n = split(t_order[gid], tids, "\t")
        for (i=1; i<=n; i++) {
            tid = tids[i]
            print trans[gid SUBSEP tid]
            printf "%s", exon[gid SUBSEP tid]
        }
    }
}
' gene_order.txt merged_temp.gtf > master_lincRNA.gtf
echo "Done: $(wc -l < master_lincRNA.gtf) lines written"