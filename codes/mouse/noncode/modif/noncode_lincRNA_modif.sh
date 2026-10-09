
#GENE

#append 'gene_type "lncRNA";' attribute to col9 in all gene entries
#append 'gene_name "gene_id";' attribute to col9 in all gene entries
#insert 'level 0;' 
#insert  'mgi_id "MGI:gene_id";'
#insert 'havana_gene "gene_id";' 
#remove all other attributes

awk 'BEGIN{FS=OFS="\t"} $3=="gene"{match($9,/gene_id "[^"]+"/); id=substr($9,RSTART+9,RLENGTH-10); n=id; sub(/\.[0-9]+$/,"",n); $9="gene_id \""id"\"; gene_type \"lncRNA\"; gene_name \""n"\"; level 0; mgi_id \"MGI:"id"\"; havana_gene \""id"\";"} 1' noncode_curated_lincRNA.gtf > noncode_curated_lincRNA_1.gtf

#TRANSCRIPT

# add 'gene_type "lncRNA";' after transcript_id 
# add 'gene_name "gene_id";' after 'gene_type "lncRNA";'
# add 'transcript_type "lncRNA";' after 'gene_name "gene_id";
# add 'transcript_name "transcript_id";' after 'transcript_type "lncRNA";' 
# add 'level 0;' after 'transcript_name "transcript_id";'
# insert 'transcript_support_level "0"'' after level 0 in all transcript entries
# insert 'mgi_id = "MGI:gene_id"'' after transcript_support_level in all transcript entries 
# insert 'tag "NONCODE_v6";' after 'mgi_id "gene_id";' in all transcript entries
#insert 'havana_gene "gene_id"'' after tag in all transcript entries 
# insert 'havana_transcript "transcript_id"'' following havana_gene in all transcript entries 

awk 'BEGIN{FS=OFS="\t"} $3=="transcript"{match($9,/gene_id "[^"]+"/); g=substr($9,RSTART+9,RLENGTH-10); match($9,/transcript_id "[^"]+"/); t=substr($9,RSTART+15,RLENGTH-16); gn=g; sub(/\.[0-9]+$/,"",gn); tn=t; sub(/\.[0-9]+$/,"",tn); $9="gene_id \""g"\"; transcript_id \""t"\"; gene_type \"lncRNA\"; gene_name \""gn"\"; transcript_type \"lncRNA\"; transcript_name \""tn"\"; level 0; transcript_support_level \"0\"; mgi_id \"MGI:"g"\"; tag \"NONCODE_v6\"; havana_gene \""g"\"; havana_transcript \""t"\";"} 1' noncode_curated_lincRNA_1.gtf > noncode_curated_lincRNA_2.gtf



#EXON 

# add 'gene_type "lncRNA"; 'gene_name "gene_id"; 'transcript_type "lncRNA"; 'transcript_name "transcript_id";'  after transcript_id
#place 'exon_number "number";' after 'transcript_name "transcript_id";' and remove "" in 'exon_number "number"
#After 'exon_number "number";' , add 'exon_id "exon_id"; level 0; hgnc_id "gene_id"; tag "NONCODE_v6"; havane_gene "gene";' wherein exon_id is derived by replacing NONHSAT with NONHSAE in transcript_id

awk 'BEGIN{FS=OFS="\t"} $3=="exon"{match($9,/gene_id "[^"]+"/); g=substr($9,RSTART+9,RLENGTH-10); match($9,/transcript_id "[^"]+"/); t=substr($9,RSTART+15,RLENGTH-16); match($9,/exon_number "[^"]+"/); en=substr($9,RSTART+13,RLENGTH-14); gn=g; sub(/\.[0-9]+$/,"",gn); tn=t; sub(/\.[0-9]+$/,"",tn); e=t; sub(/^NONMMUT/,"NONMMUE",e); $9="gene_id \""g"\"; transcript_id \""t"\"; gene_type \"lncRNA\"; gene_name \""gn"\"; transcript_type \"lncRNA\"; transcript_name \""tn"\"; exon_number "en"; exon_id \""e"\"; level 0; mgi_id \"MGI:"g"\"; tag \"NONCODE_v6\"; havana_gene \""g"\";"} 1' noncode_curated_lincRNA_2.gtf > noncode_curated_lincRNA_3.gtf


#renaming orphan gene and transcripts


sed -E 's/NONMMUT([0-9.]+)_gene/NONMMUOG\1/g; s/NONMMUG([0-9.]+)_tx/NONMMUOT\1/g' noncode_curated_lincRNA_3.gtf > noncode_curated_lincRNA_4.gtf

#Verify 
F=noncode_curated_lincRNA_4.gtf
# an orphan-transcript gene (NONMMUOG): show all its lines
id=$(awk -F'\t' '$3=="gene" && $9~/NONMMUOG/{match($9,/gene_id "[^"]+"/);print substr($9,RSTART+9,RLENGTH-10);exit}' $F); grep -F "\"$id\"" $F | cut -f3,9
# an orphan-gene transcript (NONMMUOT): show all its lines
id=$(awk -F'\t' '$3=="transcript" && $9~/transcript_id "NONMMUOT/{match($9,/transcript_id "[^"]+"/);print substr($9,RSTART+15,RLENGTH-16);exit}' $F); grep -F "\"$id\"" $F | cut -f3,9





