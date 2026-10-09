#download rna central gff3 
wget -c https://ftp.ebi.ac.uk/pub/databases/RNAcentral/releases/17.0/genome_coordinates/gff3/mus_musculus.GRCm38.gff3.gz

gunzip mus_musculus.GRCm38.gff3.gz

#Find out feature type  
awk -F'\t' '!/^#/ {print $3}' mus_musculus.GRCm38.gff3 | sort | uniq -c | sort -rn

# Above command gives noncoding_exon for exon and no gene
#replace noncoding_exon with exon
awk -F'\t' 'BEGIN{OFS="\t"} /^#/{print; next} $3=="noncoding_exon"{$3="exon"} {print}' \
  mus_musculus.GRCm38.gff3 > mus_musculus.GRCm38_fixed_exon.gff3
  
#verify
awk -F'\t' '!/^#/ {print $3}' mus_musculus.GRCm38_fixed_exon.gff3 | sort | uniq -c | sort -rn

#convert to gtf using AGAT
#Create new env for AGAT. The AGAT tool converts rna central gff3 into gtf use following sh files, Final Output file to use: mus_musculus.GRCm38_fixed_exon.gtf


awk -F'\t' '
$3=="gene"{
    if(match($9,/type "([^"]+)"/,a))
        print a[1]
}
' mus_musculus.GRCm38_fixed_exon.gtf | sort | uniq -c | sort -nr > gene_type_counts.txt

#should have chromToUCSC installed (homo_sapiens.GRCh38_fixed_v2.gtf has unconventional chr naming: 1,2,3 and not chr1, chr 2 etc) 
wget -c https://hgdownload.soe.ucsc.edu/goldenPath/mm10/database/chromAlias.txt.gz
gunzip chromAlias.txt.gz


conda activate myenv

~/chromToUcsc -a chromAlias.txt -i mus_musculus.GRCm38_fixed_exon.gtf -o mus_musculus.GRCm38_fixed_exon_chr.gtf -s

#remove chr patches
awk -F'\t' '$1 ~ /^chr([1-9]|1[0-9]|X|Y|M)$/' mus_musculus.GRCm38_fixed_exon_chr.gtf > mus_musculus.GRCm38_fixed_exon_chr.main.gtf


#Step 1. Extract all "lncRNA" - gene_id in col4

awk -F'\t' 'BEGIN{OFS="\t"} $3=="gene" && (/type "lncRNA"/ || /type "SO:0001877"/ || /type "SO:0001463"/) {
    match($9, /gene_id "([^"]+)"/, a)
    print $1, $4-1, $5, a[1]
}' mus_musculus.GRCm38_fixed_exon_chr.main.gtf > rnacentral_lncRNA_genes.bed


#Step 2. Extract all NON-lncRNA genes - nonlncRNA gene ids in col4
awk -F'\t' 'BEGIN{OFS="\t"} $3=="gene" && !(/type "lncRNA"/ || /type "SO:0001877"/ || /type "SO:0001463"/) {
    match($9, /gene_id "([^"]+)"/, a)
    print $1, $4-1, $5, a[1]
}' mus_musculus.GRCm38_fixed_exon_chr.main.gtf > rnacentral_non_lncRNA_genes.bed

#Step 3. Remove overlaps #bedtools must be installed
conda activate myenv

bedtools intersect -v -a rnacentral_lncRNA_genes.bed -b rnacentral_non_lncRNA_genes.bed gencode_non_lncRNA_genes.bed > nonoverlapping_genes.bed

#remove lncRNA subclasses (not applicable here but added to script to maintian consistency in pipeliness)
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
bedtools intersect -v -a rnacentral_curated_lincRNA.bed -b nc_minus_gc_curated_lincRNA.bed gencode_curated_lincRNA.bed > rc_minus_nc_gc_curated_lincRNA.bed

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
' rnacentral_unique_lincRNA_gene_ids.txt mus_musculus.GRCm38_fixed_exon_chr.main.gtf > rnacentral_curated_lincRNA.gtf

#21/9/26
#check feature types 
awk '{print $3}' rnacentral_curated_lincRNA.gtf | sort | uniq -c
# 18293 exon
#   8440 gene
#   8440 transcript

#7/10/26 
#check feature types

#  21572 exon
#   9932 gene
#   9932 transcript