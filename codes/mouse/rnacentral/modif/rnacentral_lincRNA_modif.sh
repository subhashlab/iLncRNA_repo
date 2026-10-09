#GENE


# replace gene_id with URS ID with 'G' before '_' , derived from Name
#remove ID attribute
#awk 'BEGIN{FS="\t"; OFS="\t"} $3=="gene" {gsub(/ID "[^"]*"; /, "", $9)} {print}' rnacentral_curated_lincRNA.gtf > rnacentral_curated_lincRNA_1.gtf
#after gene_id,  add 'gene_type "lncRNA";',  rename 'Name "gene_id";' to 'gene_name "gene_id";' , add 'level 0;', 'mgi_id "gene_id", 'havana_gene "gene_id"'  
#remove remaining column 

awk -F'\t' 'BEGIN{OFS="\t"} $3=="gene"{match($9,/Name "([^"]+)"/,n); gname=n[1]; gid=gensub(/_/,"G_",1,gname); $9="gene_id \""gid"\"; gene_type \"lncRNA\"; gene_name \""gname"\"; level 0; mgi_id \""gid"\"; havana_gene \""gid"\";"} 1' rnacentral_curated_lincRNA.gtf > rnacentral_curated_lincRNA_1.gtf



#TRANSCRIPT
#remove ID attribute after transcript id
# replace gene_id with URS ID with 'G' before '_' , derived from Name
#after transcript_id, add 'gene_type "lncRNA";' , 'gene_name "gene_id";' , 'transcript_type "lncRNA";' , 'transcript_name "transcript_id";' , 'level 0;' , ', mgi_id "gene_id";' , 'tag "NONCODE_v6";' , 'havana_gene "gene_id"; ' and remaining attributes after these attributes
#remove remaining attributes 

awk -F'\t' 'BEGIN{OFS="\t"} $3=="transcript"{match($9,/Name "([^"]+)"/,n); gid=gensub(/_/,"G_",1,n[1]); match($9,/transcript_id "([^"]+)"/,t); tid=t[1]; $9="gene_id \""gid"\"; transcript_id \""tid"\"; gene_type \"lncRNA\"; gene_name \""gid"\"; transcript_type \"lncRNA\"; transcript_name \""tid"\"; level 0; mgi_id \""gid"\"; tag \"RNA_Central_v17\"; havana_gene \""gid"\";"} 1' rnacentral_curated_lincRNA_1.gtf > rnacentral_curated_lincRNA_2.gtf


#EXON

# replace gene_id with URS ID with 'G' before '_' , derived from Name

# keep gene_id "URS..."; transcript_id "URS..."; gene_type "lncRNA"; gene_name "RNACG..."; transcript_type "lncRNA"; transcript_name "URS..."; 2. keep remaining attributes 

#prepend 'exon_number number;' after transcript_name


#remove ID

#After exon_id , add 'level 0; mgi_id "gene_id"; tag "RNA_Central_v17"; havane_gene "gene";' and keep remaining columns

#remove remaining columns after havana_gene in exon entries

awk -F'\t' 'BEGIN{OFS="\t"} $3=="exon"{match($9,/Name "([^"]+)"/,n); gid=gensub(/_/,"G_",1,n[1]); match($9,/transcript_id "([^"]+)"/,t); tid=t[1]; match($9,/ID "([^:]+):ncRNA_exon([0-9]+)"/,e); enum=e[2]; eid=e[1] e[2]; $9="gene_id \""gid"\"; transcript_id \""tid"\"; gene_type \"lncRNA\"; gene_name \""gid"\"; transcript_type \"lncRNA\"; transcript_name \""tid"\"; exon_number "enum"; exon_id \""eid"\"; level 0; mgi_id \""gid"\"; tag \"RNA_Central_v17\"; havana_gene \""gid"\";"} 1' rnacentral_curated_lincRNA_2.gtf > rnacentral_curated_lincRNA_3.gtf