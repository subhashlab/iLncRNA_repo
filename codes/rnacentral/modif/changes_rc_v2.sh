#GENE

#remove ID attribute
awk 'BEGIN{FS="\t"; OFS="\t"} $3=="gene" {gsub(/ID "[^"]*"; /, "", $9)} {print}' rnacentral_curated_lincRNA_v2.gtf > rnacentral_curated_lincRNA_v2_1.gtf

#after gene_id,  add 'gene_type "lncRNA";',  rename 'Name "gene_id";' to 'gene_name "gene_id";' , add 'level 0;', 'hgnc_id "gene_id", 'havana_gene "gene_id"'  
awk 'BEGIN{FS="\t"; OFS="\t"}
$3=="gene" {
    match($9, /gene_id "([^"]+)"/, g)
    gid = g[1]
    match($9, /Name "([^"]+)"/, n)
    gname = n[1]
    # Remove gene_id, Name, and trailing attributes we will reorder
    rest = $9
    gsub(/gene_id "[^"]*"; /, "", rest)
    gsub(/Name "[^"]*"; /, "", rest)
    $9 = "gene_id \"" gid "\"; gene_type \"lncRNA\"; gene_name \"" gname "\"; level 0; hgnc_id \"" gid "\"; havana_gene \"" gid "\"; " rest
}
{print}' rnacentral_curated_lincRNA_v2_1.gtf > rnacentral_curated_lincRNA_v2_2.gtf

#remove remaining column 
awk 'BEGIN{FS="\t"; OFS="\t"}
$3=="gene" {
    match($9, /gene_id "([^"]+)"/, g)
    gid = g[1]
    match($9, /Name "([^"]+)"/, n)
    gname = n[1]
    $9 = "gene_id \"" gid "\"; gene_type \"lncRNA\"; gene_name \"" gname "\"; level 0; hgnc_id \"" gid "\"; havana_gene \"" gid "\";"
}
{print}' rnacentral_curated_lincRNA_v2_1.gtf > rnacentral_curated_lincRNA_v2_2_1.gtf

#TRANSCRIPT

#remove ID attribute after transcript id

awk 'BEGIN{FS="\t"; OFS="\t"} $3=="transcript" {gsub(/ID "[^"]*"; /, "", $9)} {print}' rnacentral_curated_lincRNA_v2_2_1.gtf > rnacentral_curated_lincRNA_v2_3.gtf

#after transcript_id, add 'gene_type "lncRNA";' , 'gene_name "gene_id";' , 'transcript_type "lncRNA";' , 'transcript_name "transcript_id";' , 'level 0;' , 'hgnc_id "gene_id";' , 'tag "NONCODE_v6";' , 'havana_gene "gene_id"; ' and remaining attributes after these attributes

awk 'BEGIN{FS="\t"; OFS="\t"}
$3=="transcript" {
    match($9, /gene_id "([^"]+)"/, g)
    gid = g[1]
    match($9, /transcript_id "([^"]+)"/, t)
    tid = t[1]
    rest = $9
    gsub(/gene_id "[^"]*"; /, "", rest)
    gsub(/transcript_id "[^"]*"; /, "", rest)
    $9 = "gene_id \"" gid "\"; transcript_id \"" tid "\"; gene_type \"lncRNA\"; gene_name \"" gid "\"; transcript_type \"lncRNA\"; transcript_name \"" tid "\"; level 0; hgnc_id \"" gid "\"; tag \"RNA_Central_v26\"; havana_gene \"" gid "\"; " rest
}
{print}' rnacentral_curated_lincRNA_v2_3.gtf > rnacentral_curated_lincRNA_v2_4.gtf

#remove remaining attributes 
awk 'BEGIN{FS="\t"; OFS="\t"}
$3=="transcript" {
    match($9, /gene_id "([^"]+)"/, g)
    gid = g[1]
    match($9, /transcript_id "([^"]+)"/, t)
    tid = t[1]
    $9 = "gene_id \"" gid "\"; transcript_id \"" tid "\"; gene_type \"lncRNA\"; gene_name \"" gid "\"; transcript_type \"lncRNA\"; transcript_name \"" tid "\"; level 0; hgnc_id \"" gid "\"; tag \"RNA_Central_v26\"; havana_gene \"" gid "\";"
}
{print}' rnacentral_curated_lincRNA_v2_3.gtf > rnacentral_curated_lincRNA_v2_4_1.gtf

#EXON
#1. keep gene_id "RNACG..."; transcript_id "URS..."; gene_type "lncRNA"; gene_name "RNACG..."; transcript_type "lncRNA"; transcript_name "URS..."; 2. keep remaining attributes 
awk 'BEGIN{FS="\t"; OFS="\t"}
$3=="exon" {
    match($9, /gene_id "([^"]+)"/, g)
    gid = g[1]
    match($9, /transcript_id "([^"]+)"/, t)
    tid = t[1]
    rest = $9
    gsub(/gene_id "[^"]*"; /, "", rest)
    gsub(/transcript_id "[^"]*"; /, "", rest)
    $9 = "gene_id \"" gid "\"; transcript_id \"" tid "\"; gene_type \"lncRNA\"; gene_name \"" gid "\"; transcript_type \"lncRNA\"; transcript_name \"" tid "\"; " rest
}
{print}' rnacentral_curated_lincRNA_v2_4_1.gtf > rnacentral_curated_lincRNA_v2_5.gtf

#prepend 'exon_number number;' after transcript_name
awk 'BEGIN{FS="\t"; OFS="\t"}
$3=="exon" {
    match($9, /gene_id "([^"]+)"/, g)
    gid = g[1]
    match($9, /transcript_id "([^"]+)"/, t)
    tid = t[1]
    match($9, /ID "[^"]*:ncRNA_exon([0-9]+)"/, e)
    enum = e[1]
    rest = $9
    gsub(/gene_id "[^"]*"; /, "", rest)
    gsub(/transcript_id "[^"]*"; /, "", rest)
    gsub(/gene_type "[^"]*"; /, "", rest)
    gsub(/gene_name "[^"]*"; /, "", rest)
    gsub(/transcript_type "[^"]*"; /, "", rest)
    gsub(/transcript_name "[^"]*"; /, "", rest)
    $9 = "gene_id \"" gid "\"; transcript_id \"" tid "\"; gene_type \"lncRNA\"; gene_name \"" gid "\"; transcript_type \"lncRNA\"; transcript_name \"" tid "\"; exon_number " enum "; " rest
}
{print}' rnacentral_curated_lincRNA_v2_5.gtf > rnacentral_curated_lincRNA_v2_6.gtf
#add 'exon_id "exon_id";' after 'exon_number number' from ID "URS0002853A2B_9606.1:ncRNA_exon2" = ID "URS0002853A2B_9606.12

awk 'BEGIN{FS="\t"; OFS="\t"}
$3=="exon" {
    match($9, /gene_id "([^"]+)"/, g)
    gid = g[1]
    match($9, /transcript_id "([^"]+)"/, t)
    tid = t[1]
    match($9, /ID "([^:]+):ncRNA_exon([0-9]+)"/, e)
    eid = e[1] e[2]
    enum = e[2]
    rest = $9
    gsub(/gene_id "[^"]*"; /, "", rest)
    gsub(/transcript_id "[^"]*"; /, "", rest)
    gsub(/gene_type "[^"]*"; /, "", rest)
    gsub(/gene_name "[^"]*"; /, "", rest)
    gsub(/transcript_type "[^"]*"; /, "", rest)
    gsub(/transcript_name "[^"]*"; /, "", rest)
    gsub(/exon_number [0-9]+; /, "", rest)
    gsub(/exon_id "[^"]*"; /, "", rest)
    $9 = "gene_id \"" gid "\"; transcript_id \"" tid "\"; gene_type \"lncRNA\"; gene_name \"" gid "\"; transcript_type \"lncRNA\"; transcript_name \"" tid "\"; exon_number " enum "; exon_id \"" eid "\"; " rest
}
{print}' rnacentral_curated_lincRNA_v2_6.gtf > rnacentral_curated_lincRNA_v2_7.gtf

#remove ID
awk 'BEGIN{FS="\t"; OFS="\t"} $3=="exon" {gsub(/ID "[^"]*"; /, "", $9)} {print}' rnacentral_curated_lincRNA_v2_7.gtf > rnacentral_curated_lincRNA_v2_8.gtf
#After exon_id , add 'level 0; hgnc_id "gene_id"; tag "RNA_Central_v26"; havane_gene "gene";' and keep remaining columns
awk 'BEGIN{FS="\t"; OFS="\t"}
$3=="exon" {
    match($9, /gene_id "([^"]+)"/, g)
    gid = g[1]
    match($9, /transcript_id "([^"]+)"/, t)
    tid = t[1]
    match($9, /exon_number ([0-9]+)/, e)
    enum = e[1]
    match($9, /exon_id "([^"]+)"/, x)
    eid = x[1]
    rest = $9
    gsub(/gene_id "[^"]*"; /, "", rest)
    gsub(/transcript_id "[^"]*"; /, "", rest)
    gsub(/gene_type "[^"]*"; /, "", rest)
    gsub(/gene_name "[^"]*"; /, "", rest)
    gsub(/transcript_type "[^"]*"; /, "", rest)
    gsub(/transcript_name "[^"]*"; /, "", rest)
    gsub(/exon_number [0-9]+; /, "", rest)
    gsub(/exon_id "[^"]*"; /, "", rest)
    $9 = "gene_id \"" gid "\"; transcript_id \"" tid "\"; gene_type \"lncRNA\"; gene_name \"" gid "\"; transcript_type \"lncRNA\"; transcript_name \"" tid "\"; exon_number " enum "; exon_id \"" eid "\"; level 0; hgnc_id \"" gid "\"; tag \"RNA_Central_v26\"; havana_gene \"" gid "\"; " rest
}
{print}' rnacentral_curated_lincRNA_v2_8.gtf > rnacentral_curated_lincRNA_v2_9.gtf
#remove remaining columns after havana_gene in exon entries
awk 'BEGIN{FS="\t"; OFS="\t"}
$3=="exon" {
    match($9, /gene_id "([^"]+)"/, g)
    gid = g[1]
    match($9, /transcript_id "([^"]+)"/, t)
    tid = t[1]
    match($9, /exon_number ([0-9]+)/, e)
    enum = e[1]
    match($9, /exon_id "([^"]+)"/, x)
    eid = x[1]
    $9 = "gene_id \"" gid "\"; transcript_id \"" tid "\"; gene_type \"lncRNA\"; gene_name \"" gid "\"; transcript_type \"lncRNA\"; transcript_name \"" tid "\"; exon_number " enum "; exon_id \"" eid "\"; level 0; hgnc_id \"" gid "\"; tag \"RNA_Central_v26\"; havana_gene \"" gid "\";"
}
{print}' rnacentral_curated_lincRNA_v2_9.gtf > rnacentral_curated_lincRNA_v2_9_1.gtf
