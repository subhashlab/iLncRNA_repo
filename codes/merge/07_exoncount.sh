awk -F'\t' '
$3=="exon" {
    match($9, /transcript_id "[^"]+"/)
    t = substr($9, RSTART+15, RLENGTH-16)
    n[t]++
}
END {
    for (t in n) {
        x = n[t]
        if (x > 5) b = ">5"
        else       b = x
        c[b]++
    }
    m = split("1,2,3,4,5,>5", o, ",")
    print "Exon count\tTranscripts"
    for (i = 1; i <= m; i++) print o[i] "\t" c[o[i]]+0
}' master_lincRNA_hg38.gtf

#same for gencode_curated_lincRNA.sh and gencode_curated_lincRNA_v2.sh