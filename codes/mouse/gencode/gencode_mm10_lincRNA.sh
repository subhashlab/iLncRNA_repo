#Download gtf file
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M23/gencode.vM23.primary_assembly.annotation.gtf.gz

gunzip gencode.vM23.primary_assembly.annotation.gtf.gz


#Check the gene types present
awk '$3=="gene"' gencode.vM23.primary_assembly.annotation.gtf \
| grep -o 'gene_type "[^"]*"' \
| sort | uniq -c | sort -nr


#Step 1. Extract all "lncRNA" - gene_id in col4, gene_name in col5

awk '
BEGIN{OFS="\t"}
$3=="gene" && $0~/gene_type "lncRNA"/ {
    gene="NA"
    if(match($0,/gene_id "[^"]+"/))
        gene=substr($0,RSTART+9,RLENGTH-10)
    name="NA"
    if(match($0,/gene_name "[^"]+"/))
        name=substr($0,RSTART+11,RLENGTH-12)
    print $1,$4-1,$5,gene,name
}
' gencode.vM23.primary_assembly.annotation.gtf \
| grep -E "^chr([0-9]+|X|Y)[[:space:]]" \
| sort -k1,1 -k2,2n \
> all_lncRNA.bed

#Step 2. Extract all NON-lncRNA genes - nonlncRNA gene ids in col4, gene_name in col5

awk '
BEGIN{OFS="\t"}
$3=="gene" && $0!~/gene_type "lncRNA"/ {
    gene="NA"
    if(match($0,/gene_id "[^"]+"/))
        gene=substr($0,RSTART+9,RLENGTH-10)
    print $1,$4-1,$5,gene
}
' gencode.vM23.primary_assembly.annotation.gtf \
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

$5 ~ /-os[0-9]*$/ {next}


$5 ~ /as[0-9]*$/ {next}

toupper($5) ~ /^MIR/ {next}

toupper($5) ~ /^SNHG/ {next}

toupper($5) ~ /^SNORD/ {next}
toupper($5) ~ /^SNORA/ {next}

toupper($5) ~ /^RNU/ {next}

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
' gencode_unique_lincRNA_gene_ids.txt  gencode.vM23.primary_assembly.annotation.gtf > gencode_curated_lincRNA.gtf 



#verified on 12/9/26
#check feature types
 awk '{print $3}' gencode_curated_lincRNA.gtf | sort | uniq -c
#23309 exon
#5219 gene
#7598 transcript



#check unique gene ids
grep -o 'gene_id "[^"]*"' gencode_curated_lincRNA.gtf | sort -u | wc -l
#5219

