#download rna central gff3 
wget -c https://ftp.ebi.ac.uk/pub/databases/RNAcentral/releases/26.0/genome_coordinates/gff3/homo_sapiens.GRCh38.gff3.gz

gunzip homo_sapiens.GRCh38.gff3.gz

#replace noncoding_exon and predicted_gene with exon and gene. 
 # Step 1: noncoding_exon -> exon
awk -F'\t' 'BEGIN{OFS="\t"} $3=="noncoding_exon"{$3="exon"} {print}' \
  homo_sapiens.GRCh38.gff3 > homo_sapiens.GRCh38_fixed.gff3

 # Step 2: predicted_gene -> gene
awk -F'\t' 'BEGIN{OFS="\t"} $3=="predicted_gene"{$3="gene"} {print}' \
  homo_sapiens.GRCh38_fixed.gff3 > homo_sapiens.GRCh38_fixed_v2.gff3

#Create new env for AGAT. The AGAT tool converts rna central gff3 into gtf use following sh files, Final Output file to use: homo_sapiens.GRCh38_fixed_v2.gtf
bash agat_convert3.sh 

 
awk -F'\t' '
$3=="gene"{
    if(match($9,/type "([^"]+)"/,a))
        print a[1]
}
' homo_sapiens.GRCh38_fixed_v2.gtf | sort | uniq -c | sort -nr > gene_type_counts.txt

#should have chromToUCSC installed (homo_sapiens.GRCh38_fixed_v2.gtf has unconventional chr naming: 1,2,3 and not chr1, chr 2 etc) 
conda activate myenv
~/chromToUcsc --get hg38
 ~/chromToUcsc -a ./hg38.chromAlias.tsv -i homo_sapiens.GRCh38_fixed_v2.gtf -o homo_sapiens.GRCh38_fixed_v2_chr.gtf -s


#Step 1. Extract all "lncRNA" - gene_id in col4

awk -F'\t' 'BEGIN{OFS="\t"} $3=="gene" && (/type "lncRNA"/ || /type "SO:0001877"/ || /type "SO:0001463"/) {
    match($9, /gene_id "([^"]+)"/, a)
    print $1, $4-1, $5, a[1]
}' homo_sapiens.GRCh38_fixed_v2_chr.gtf > rnacentral_lncRNA_genes.bed


#Step 2. Extract all NON-lncRNA genes - nonlncRNA gene ids in col4
awk -F'\t' 'BEGIN{OFS="\t"} $3=="gene" && !(/type "lncRNA"/ || /type "SO:0001877"/ || /type "SO:0001463"/) {
    match($9, /gene_id "([^"]+)"/, a)
    print $1, $4-1, $5, a[1]
}' homo_sapiens.GRCh38_fixed_v2_chr.gtf > rnacentral_non_lncRNA_genes.bed

#Step 3. Remove overlaps #bedtools must be installed
conda activate myenv

bedtools intersect -v -a rnacentral_lncRNA_genes.bed -b rnacentral_non_lncRNA_genes.bed gencode_non_lncRNA_genes.bed > nonoverlapping_genes.bed


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
> rnacentral_curated_lincRNA.bed


#remove gencode_curated_lincRNA.bed and nc_minus_gc_curated_lincRNA.bed gene entries
bedtools intersect -v -a rnacentral_curated_lincRNA_1.bed -b nc_minus_gc_curated_lincRNA.bed gencode_curated_lincRNA.bed > rc_minus_nc_gc_curated_lincRNA.bed
bedtools intersect -a rc_minus_nc_gc_curated_lincRNA.bed -b nc_minus_gc_curated_lincRNA.bed gencode_curated_lincRNA.bed -wa -wb | head
bedtools intersect -a rc_minus_nc_gc_curated_lincRNA.bed -b nc_minus_gc_curated_lincRNA.bed gencode_curated_lincRNA.bed -u | wc -l


#segregate unique gene ids
awk '{print $4}' rc_minus_nc_gc_curated_lincRNA.bed | sort -u > rnacentral_unique_lincRNA_gene_ids.txt
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
' rnacentral_unique_lincRNA_gene_ids.txt homo_sapiens.GRCh38_fixed_v2_chr.gtf > rnacentral_curated_lincRNA_v2.gtf


#check feature types 

awk '{print $3}' rnacentral_curated_lincRNA_v2.gtf | sort | uniq -c
  # 6632 exon
  # 1124 gene
  # 2730 transcript


#check unique gene ids 

grep -o 'gene_id "[^"]*"' rnacentral_curated_lincRNA_v2.gtf | sort -u | wc -l
#1124

