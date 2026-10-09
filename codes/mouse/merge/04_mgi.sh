#!/bin/bash
gawk -F'\t' 'BEGIN{OFS="\t"}
$2!="HAVANA" && $2!="ENSEMBL" {
    match($9, /gene_id "([^"]+)"/, g)
    gid = g[1]
    if ($9 ~ /mgi_id "[^"]*"/) {
        sub(/mgi_id "[^"]*"/, "mgi_id \"MGI:" gid "\"", $9)
    } else {
        $9 = $9 " mgi_id \"MGI:" gid "\";"
    }
}
{print}' master_lincRNA_1.1.gtf > mm_master_iLncRNA_annotation_v1.0.gtf
