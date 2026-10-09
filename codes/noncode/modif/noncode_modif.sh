#GENE

#append 'gene_type "lncRNA";' attribute to col9 in all gene entries
awk 'BEGIN{OFS="\t"} $3=="gene" {sub(/gene_id "[^"]+";/, "& gene_type \"lncRNA\";"); print; next} {print}' noncode_curated_lincRNA.gtf > noncode_curated_lincRNA_1.gtf

#append 'gene_name "gene_id";' attribute to col9 in all gene entries
awk -F'\t' 'BEGIN{OFS="\t"} $3=="gene" {
    match($9, /gene_id "([^"]+)"/, a);
    $9=$9" gene_name \""a[1]"\";";
    print
} $3!="gene" {print}' noncode_curated_lincRNA_1.gtf > noncode_curated_lincRNA_2.gtf

#append 'level 0;' 
awk -F'\t' 'BEGIN{OFS="\t"} $3=="gene" {$9=$9" level 0;"; print} $3!="gene" {print}' noncode_curated_lincRNA_2.gtf > noncode_curated_lincRNA_3.gtf

#append  'hgnc_id "gene_id";'
awk -F'\t' 'BEGIN{OFS="\t"} $3=="gene" {match($9,/gene_id "([^"]+)"/,a); $9=$9" hgnc_id \""a[1]"\";"; print} $3!="gene" {print}' noncode_curated_lincRNA_3.gtf > noncode_curated_lincRNA_4.gtf

#append 'havana_gene "gene_id";' 
awk -F'\t' 'BEGIN{OFS="\t"} $3=="gene" {match($9,/gene_id "([^"]+)"/,a); $9=$9" havana_gene \""a[1]"\";"; print} $3!="gene" {print}' noncode_curated_lincRNA_4.gtf > noncode_curated_lincRNA_5.gtf

#TRANSCRIPT

# add 'gene_type "lncRNA";' after transcript_id 
awk -F'\t' 'BEGIN{OFS="\t"} $3=="transcript" {sub(/transcript_id "[^"]+";/, "& gene_type \"lncRNA\";"); print} $3!="transcript" {print}' noncode_curated_lincRNA_5.gtf > noncode_curated_lincRNA_6.gtf

# add 'gene_name "gene_id";' after 'gene_type "lncRNA";'
awk -F'\t' 'BEGIN{OFS="\t"} $3=="transcript" {match($9,/gene_id "([^"]+)"/,a); sub(/gene_type "lncRNA";/, "& gene_name \""a[1]"\";"); print} $3!="transcript" {print}' noncode_curated_lincRNA_6.gtf > noncode_curated_lincRNA_7.gtf

# add 'transcript_type "lncRNA";' after 'gene_name "gene_id";
awk -F'\t' 'BEGIN{OFS="\t"} $3=="transcript" {sub(/gene_name "[^"]+";/, "& transcript_type \"lncRNA\";"); print} $3!="transcript" {print}' noncode_curated_lincRNA_7.gtf > noncode_curated_lincRNA_8.gtf

# add 'transcript_name "transcript_id";' after 'transcript_type "lncRNA";' 

awk -F'\t' 'BEGIN{OFS="\t"} $3=="transcript" {match($9,/transcript_id "([^"]+)"/,a); sub(/transcript_type "lncRNA";/, "& transcript_name \""a[1]"\";"); print} $3!="transcript" {print}' noncode_curated_lincRNA_8.gtf > noncode_curated_lincRNA_9.gtf

# add 'level 0;' after 'transcript_name "transcript_id";'

awk -F'\t' 'BEGIN{OFS="\t"} $3=="transcript" {match($9,/transcript_name "([^"]+)"/,a); sub(/transcript_name "[^"]+";/, "& level 0;"); print} $3!="transcript" {print}' noncode_curated_lincRNA_9.gtf > noncode_curated_lincRNA_10.gtf

# add 'hgnc_id "gene_id";' after 'level 0;'
awk -F'\t' 'BEGIN{OFS="\t"} $3=="transcript" {match($9,/gene_id "([^"]+)"/,a); sub(/level 0;/, "& hgnc_id \""a[1]"\";"); print} $3!="transcript" {print}' noncode_curated_lincRNA_10.gtf > noncode_curated_lincRNA_11.gtf

# add 'tag "NONCODE_v6";' after 'hgnc_id "gene_id";'

awk -F'\t' 'BEGIN{OFS="\t"} $3=="transcript" {sub(/hgnc_id "[^"]+";/, "& tag \"NONCODE_v6\";"); print} $3!="transcript" {print}' noncode_curated_lincRNA_11.gtf > noncode_curated_lincRNA_12.gtf

# add 'havana_gene "gene_id" after tag "NONCODE_v6";' 
awk -F'\t' 'BEGIN{OFS="\t"} $3=="transcript" {match($9,/gene_id "([^"]+)"/,a); sub(/tag "NONCODE_v6";/, "& havana_gene \""a[1]"\";"); print} $3!="transcript" {print}' noncode_curated_lincRNA_12.gtf > noncode_curated_lincRNA_13.gtf

#remove 'FPKM "0"; exon_number "";' from transcript entries 
awk -F'\t' 'BEGIN{OFS="\t"} $3=="transcript" {sub(/ FPKM "0"; exon_number "[^"]*";/, ""); print} $3!="transcript" {print}' noncode_curated_lincRNA_13.gtf > noncode_curated_lincRNA_13_1.gtf

#EXON

# add 'gene_type "lncRNA"; 'gene_name "gene_id"; 'transcript_type "lncRNA"; 'transcript_name "transcript_id";'  after transcript_id

awk -F'\t' 'BEGIN{OFS="\t"} $3=="exon" {match($9,/gene_id "([^"]+)"/,a); match($9,/transcript_id "([^"]+)"/,b); sub(/transcript_id "[^"]+";/, "& gene_type \"lncRNA\"; gene_name \""a[1]"\"; transcript_type \"lncRNA\"; transcript_name \""b[1]"\";"); print} $3!="exon" {print}' noncode_curated_lincRNA_13_1.gtf > noncode_curated_lincRNA_14.gtf

#place 'exon_number "number";' after 'transcript_name "transcript_id";' and remove "" in 'exon_number "number"

awk -F'\t' 'BEGIN{OFS="\t"} $3=="exon" {
    match($9,/exon_number "([^"]+)"/,a);
    sub(/ exon_number "[^"]*";/, "");
    sub(/transcript_name "[^"]+";/, "& exon_number "a[1]";");
    print
} $3!="exon" {print}' noncode_curated_lincRNA_14.gtf > noncode_curated_lincRNA_15.gtf


#After 'exon_number "number";' , add 'exon_id "exon_id"; level 0; hgnc_id "gene_id"; tag "NONCODE_v6"; havane_gene "gene";' wherein exon_id is derived by replacing NONHSAT with NONHSAE in transcript_id

awk -F'\t' 'BEGIN{OFS="\t"} $3=="exon" {
    match($9,/transcript_id "([^"]+)"/,t);
    match($9,/gene_id "([^"]+)"/,g);
    exon_id=t[1]; gsub(/NONHSAT/,"NONHSAE",exon_id);
    sub(/exon_number [0-9]+;/, "& exon_id \""exon_id"\"; level 0; hgnc_id \""g[1]"\"; tag \"NONCODE_v6\"; havana_gene \""g[1]"\";");
    print
} $3!="exon" {print}' noncode_curated_lincRNA_15.gtf > noncode_curated_lincRNA_16.gtf


#remove FPKM
awk -F'\t' 'BEGIN{OFS="\t"} $3=="exon" {
    match($9,/transcript_id "([^"]+)"/,t);
    match($9,/gene_id "([^"]+)"/,g);
    match($9,/exon_number ([0-9]+)/,e);
    exon_id=t[1]; gsub(/NONHSAT/,"NONHSAE",exon_id);
    sub(/ FPKM "0";/, "");
    sub(/exon_number [0-9]+/, "& exon_id \""exon_id"\"; level 0; hgnc_id \""g[1]"\"; tag \"NONCODE_v6\"; havana_gene \""g[1]"\"");
    print
} $3!="exon" {print}' noncode_curated_lincRNA_15.gtf > noncode_curated_lincRNA_16_1.gtf

