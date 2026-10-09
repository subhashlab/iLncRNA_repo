#TRANSCRIPTS

#verify number of transcripts with seq length distribution

awk -F'\t' '
$3=="exon" {
    match($9, /transcript_id "[^"]+"/)
    t = substr($9, RSTART+15, RLENGTH-16)
    L[t] += $5 - $4 + 1
}
END {
    for (t in L) {
        l = L[t]
        if      (l < 200)    b = "<200"
        else if (l <= 500)   b = "200-500"
        else if (l <= 1000)  b = "500-1000"
        else if (l <= 2000)  b = "1000-2000"
        else if (l <= 5000)  b = "2000-5000"
        else if (l <= 10000) b = "5000-10000"
        else                 b = ">10000"
        c[b]++
    }
    n = split("<200,200-500,500-1000,1000-2000,2000-5000,5000-10000,>10000", o, ",")
    print "Seq Length(nt)\tTranscripts"
    for (i = 1; i <= n; i++) print o[i] "\t" c[o[i]]+0
}' master_lincRNA_hg38.gtf


#GENES

#same for gencode_curated_lincRNA.gtf and gencode_curated_lincRNA_v2.gtf

awk -F'\t' '
$3=="exon" {
    match($9, /gene_id "[^"]+"/)
    g = substr($9, RSTART+9, RLENGTH-10)
    if (!(g in s) || $4 < s[g]) s[g] = $4
    if ($5 > e[g]) e[g] = $5
}
END {
    for (g in s) {
        l = e[g] - s[g] + 1
        if      (l < 200)     b = "<200"
        else if (l <= 500)    b = "200-500"
        else if (l <= 1000)   b = "500-1000"
        else if (l <= 2000)   b = "1000-2000"
        else if (l <= 5000)   b = "2000-5000"
        else if (l <= 10000)  b = "5000-10000"
        else if (l <= 50000)  b = "10000-50000"
        else if (l <= 100000) b = "50000-100000"
        else                  b = ">100000"
        c[b]++
    }
    n = split("<200,200-500,500-1000,1000-2000,2000-5000,5000-10000,10000-50000,50000-100000,>100000", o, ",")
    print "Seq Length(nt)\tGenes"
    for (i = 1; i <= n; i++) print o[i] "\t" c[o[i]]+0
}' master_lincRNA_hg38.gtf


#same for gencode_curated_lincRNA.gtf and gencode_curated_lincRNA_v2.gtf


#EXONS

awk -F'\t' '
$3=="exon" {
    l = $5 - $4 + 1
    if      (l < 50)    b = "<50"
    else if (l <= 100)  b = "50-100"
    else if (l <= 200)  b = "100-200"
    else if (l <= 300)  b = "200-300"
    else if (l <= 500)  b = "300-500"
    else if (l <= 1000) b = "500-1000"
    else                b = ">1000"
    c[b]++
}
END {
    n = split("<50,50-100,100-200,200-300,300-500,500-1000,>1000", o, ",")
    print "Seq Length(nt)\tExons"
    for (i = 1; i <= n; i++) print o[i] "\t" c[o[i]]+0
}' master_lincRNA_hg38.gtf

#same for gencode_curated_lincRNA.gtf and gencode_curated_lincRNA_v2.gtf
