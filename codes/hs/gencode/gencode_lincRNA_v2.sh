wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/gencode.v49.primary_assembly.annotation.gtf.gz
gunzip gencode.v49.primary_assembly.annotation.gtf.gz

#Check the gene types present
awk '$3=="gene"' gencode.v49.primary_assembly.annotation.gtf \
| grep -o 'gene_type "[^"]*"' \
| sort | uniq -c | sort -nr


#Step 1. Extract all "lncRNA" - gene name in col4

awk 'BEGIN{FS=OFS="\t"} $3=="gene" && $1~/^chr([0-9]+|X|Y)$/ && /gene_type "lncRNA"/ {n="NA"; if(match($0,/gene_name "[^"]+"/)) n=substr($0,RSTART+11,RLENGTH-12); print $1,$4-1,$5,n}' gencode.v49.primary_assembly.annotation.gtf | sort -k1,1 -k2,2n > all_lncRNA.bed

#Step 2. Extract all NON-lncRNA genes - nonlncRNA gene name in col4

awk 'BEGIN{FS=OFS="\t"} $3=="gene" && $1~/^chr([0-9]+|X|Y)$/ && !/gene_type "lncRNA"/ {n="NA"; if(match($0,/gene_name "[^"]+"/)) n=substr($0,RSTART+11,RLENGTH-12); print $1,$4-1,$5,n}' gencode.v49.primary_assembly.annotation.gtf | sort -k1,1 -k2,2n > gencode_non_lncRNA_genes.bed

#Step 3. Remove overlaps #bedtools must be installed
conda activate myenv

bedtools intersect -v \
-a all_lncRNA.bed \
-b gencode_non_lncRNA_genes.bed \
> nonoverlapping_lncRNA.bed

#Step 4. Remove unwanted lncRNA subclasses

awk 'BEGIN{FS=OFS="\t"} $4~/-AS[0-9]*$/ || toupper($4)~/^(MIR|SNHG|SNORD|SNORA|RNU)/ {next} 1' nonoverlapping_lncRNA.bed | sort -u > gencode_curated_lincRNA.bed

#segregate unique gene ids 
awk -F'\t' 'NR==FNR{k[$1"\t"$2"\t"$3"\t"$4];next} $3=="gene"{match($9,/gene_id "([^"]+)"/,i); match($9,/gene_name "([^"]+)"/,n); if(($1"\t"$4-1"\t"$5"\t"n[1]) in k) print i[1]}' gencode_curated_lincRNA.bed gencode.v49.primary_assembly.annotation.gtf | sort -u > gencode_unique_lincRNA_gene_ids.txt

#create gtf
awk '
NR==FNR{
    keep[$1]
    next
}

{
    if(match($0,/gene_id "([^"]+)"/,a))
        if(a[1] in keep)
            print
}
' gencode_unique_lincRNA_gene_ids.txt  gencode.v49.primary_assembly.annotation.gtf > gencode_curated_lincRNA.gtf 

#Remove:
#1. Single Exon Transcripts
#2. Remove transcripts with length less than 200 nt (including genes less than 200bp)


  
  
#check feature types
 awk '{print $3}' gencode_curated_lincRNA.gtf | sort | uniq -c
#177287 exon
#  16003 gene
#  56497 transcript

grep -o 'gene_id "[^"]*"' gencode_curated_lincRNA.gtf | sort -u | wc -l
#16003


