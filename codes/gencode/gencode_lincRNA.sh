#Download gtf file
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/gencode.v49.primary_assembly.annotation.gtf.gz

gunzip gencode.v49.primary_assembly.annotation.gtf.gz

#Check the gene types present
awk '$3=="gene"' gencode.v49.primary_assembly.annotation.gtf \
| grep -o 'gene_type "[^"]*"' \
| sort | uniq -c | sort -nr


#Step 1. Extract all "lncRNA" - gene_id in col4

awk '
BEGIN{OFS="\t"}
$3=="gene" && $0~/gene_type "lncRNA"/ {
    gene="NA"
    if(match($0,/gene_id "[^"]+"/))
        gene=substr($0,RSTART+9,RLENGTH-10)
    print $1,$4-1,$5,gene
}
' gencode.v49.primary_assembly.annotation.gtf \
| grep -E "^chr([0-9]+|X|Y)[[:space:]]" \
| sort -k1,1 -k2,2n \
> all_lncRNA.bed

#Step 2. Extract all NON-lncRNA genes - nonlncRNA gene ids in col4

awk '
BEGIN{OFS="\t"}
$3=="gene" && $0!~/gene_type "lncRNA"/ {
    gene="NA"
    if(match($0,/gene_id "[^"]+"/))
        gene=substr($0,RSTART+9,RLENGTH-10)
    print $1,$4-1,$5,gene
}
' gencode.v49.primary_assembly.annotation.gtf \
| grep -E "^chr([0-9]+|X|Y)[[:space:]]" \
| sort -k1,1 -k2,2n \
> gencode_non_lncRNA_genes.bed

#Step 3. Remove overlaps #bedtools must be installed
conda activate myenv

bedtools intersect -v \
-a all_lncRNA.bed \
-b gencode_non_lncRNA_genes.bed \
> nonoverlapping_lncRNA.bed

#Step 4. Remove unwanted lncRNA subclasses

awk '
BEGIN{FS=OFS="\t"}

$4 ~ /-AS[0-9]*$/ {next}

toupper($4) ~ /^MIR/ {next}

toupper($4) ~ /^SNHG/ {next}

toupper($4) ~ /^SNORD/ {next}
toupper($4) ~ /^SNORA/ {next}

toupper($4) ~ /^RNU/ {next}

{print}

' nonoverlapping_lncRNA.bed \
| sort -u \
> gencode_curated_lincRNA.bed

#segregate unique gene ids 
awk '{print $4}' gencode_curated_lincRNA.bed | sort -u > gencode_unique_lincRNA_gene_ids.txt

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

#check feature types
 awk '{print $3}' gencode_curated_lincRNA.gtf | sort | uniq -c
# 178089 exon
#  16029 gene
#  56707 transcript
#check unique gene ids
grep -o 'gene_id "[^"]*"' gencode_curated_lincRNA.gtf | sort -u | wc -l
#16029
