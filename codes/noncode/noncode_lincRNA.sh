#download NONCODEv6_human_hg38_lncRNA.gtf 
wget -c https://v7.noncode.org/datadownload/NONCODEv6_human_hg38_lncRNA.gtf.gz

#NONCODEv6_human_hg38_lncRNA.gtf doesn't consist of gene entries, in order to include gene entries, use gffutils tool. Output file: NONCODE_hg38_with_genes.gtf
# run noncode_gffutils_5.sh

#NONCODE_hg38_with_genes.gtf only consists of lncRNA 

#Step 1. Extract all "lncRNA" - gene_id in col4

awk -F'\t' 'BEGIN{OFS="\t"} $3=="gene" {
    match($9, /gene_id "([^"]+)"/, a)
    print $1, $4-1, $5, a[1]
}' NONCODE_hg38_with_genes.gtf > NONCODE_gene_level.bed

#Step 2. Use non lncRNA bed file from gencode and rna central since NONCODE doesn't have any non lncRNA - keep them in the same directory where you are working with noncode gtf

#Step 3. Remove overlaps #bedtools must be installed , IMP: create non_lncRNA bed file for rna central first (refer script rnacentral_lincRNA.sh)
conda activate myenv

bedtools intersect -v -a NONCODE_gene_level.bed -b gencode_non_lncRNA_genes.bed rnacentral_non_lncRNA_genes.bed > nonoverlapping_genes.bed


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

' nonoverlapping_genes.bed \
| sort -u \
> noncode_curated_lincRNA.bed
#remove gencode_curated_lincRNA.bed gene entries
bedtools intersect -v -a noncode_curated_lincRNA.bed -b gencode_curated_lincRNA.bed > nc_minus_gc_curated_lincRNA.bed

bedtools intersect -a nc_minus_gc_curated_lincRNA.bed -b gencode_curated_lincRNA.bed -wa -wb | head
bedtools intersect -a nc_minus_gc_curated_lincRNA.bed -b gencode_curated_lincRNA.bed -u | wc -l
#segregate unique gene ids 
awk '{print $4}' nc_minus_gc_curated_lincRNA.bed | sort -u > noncode_unique_lincRNA_gene_ids.txt

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
' noncode_unique_lincRNA_gene_ids.txt NONCODE_hg38_with_genes.gtf > noncode_curated_lincRNA.gtf

#check feature types
awk '{print $3}' noncode_curated_lincRNA.gtf | sort | uniq -c
# 22709 exon
#  12530 gene
#  14545 transcript


#check unique gene ids
grep -o 'gene_id "[^"]*"' noncode_curated_lincRNA.gtf | sort -u | wc -l
#12530
