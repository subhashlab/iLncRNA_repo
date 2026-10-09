#download rna central gff3 
wget -c https://ftp.ebi.ac.uk/pub/databases/RNAcentral/releases/26.0/genome_coordinates/gff3/homo_sapiens.GRCh38.gff3.gz

gunzip homo_sapiens.GRCh38.gff3.gz

#replace noncoding_exon and predicted_gene with exon and gene. 


#Create new env for AGAT. The AGAT tool converts rna central gff3 into gtf use following sh files, Final Output file to use: homo_sapiens.GRCh38_fixed_v2.gtf

#verify type of genes
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

#prepare description.txt for lncRNA type 
awk -F'\t' 'BEGIN{OFS="\t"} $3=="gene" && (/type "lncRNA"/ || /type "SO:0001877"/ || /type "SO:0001463"/) {
    match($9, /gene_id "([^"]+)"/, a)
    match($9, /description "([^"]+)"/, b)
    print b[1]
}' homo_sapiens.GRCh38_fixed_v2_chr.gtf > description_lncRNA.txt


#Step 1. Extract all "lncRNA" - gene_id in col4, description in col5

awk -F'\t' 'BEGIN{OFS="\t"} $3=="gene" && (/type "lncRNA"/ || /type "SO:0001877"/ || /type "SO:0001463"/) {
    match($9, /gene_id "([^"]+)"/, a)
    match($9, /description "([^"]+)"/, b)
    print $1, $4-1, $5, a[1], b[1]
}' homo_sapiens.GRCh38_fixed_v2_chr.gtf > rnacentral_lncRNA_genes.bed

#Step 2. Extract all NON-lncRNA genes - nonlncRNA gene ids in col4, description in col5
awk -F'\t' 'BEGIN{OFS="\t"} $3=="gene" && !(/type "lncRNA"/ || /type "SO:0001877"/ || /type "SO:0001463"/) {
    match($9, /gene_id "([^"]+)"/, a)
    match($9, /description "([^"]+)"/, b)
    print $1, $4-1, $5, a[1], b[1]
}' homo_sapiens.GRCh38_fixed_v2_chr.gtf > rnacentral_non_lncRNA_genes.bed

#Step 3. Remove overlaps #bedtools must be installed
conda activate myenv

bedtools intersect -v -a rnacentral_lncRNA_genes.bed -b rnacentral_non_lncRNA_genes.bed gencode_non_lncRNA_genes.bed > nonoverlapping_genes.bed

#prepare description_removed.txt - list of lncRNA subtypes which are not purely intergenic

grep -iE -e 'antisense' -e 'opposite strand' -e '-AS[0-9]+([: ]|$)' -e '\bMIR[0-9][A-Z0-9]*\b' -e 'microRNA' -e '\bSNHG[0-9]*\b' -e 'small nucleolar RNA host gene' -e '\bSNORD[0-9]*\b' -e '\bSNORA[0-9]*\b' -e '\bRNU[0-9]+[A-Z0-9]*\b' -e 'U[0-9]+ small nuclear RNA' description_lncRNA.txt > descriptions_removed.txt


#remove lncRNA subgroups
awk -F'\t' 'BEGIN{OFS="\t"; while((getline line < "descriptions_removed.txt") > 0) rm[line]=1}
    !($5 in rm)' nonoverlapping_genes.bed > rnacentral_curated_lincRNA.bed


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
' rnacentral_unique_lincRNA_gene_ids.txt homo_sapiens.GRCh38_fixed_v2_chr.gtf > rnacentral_curated_lincRNA.gtf

#check feature types 
awk '{print $3}' rnacentral_curated_lincRNA.gtf | sort | uniq -c


