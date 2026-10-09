#sort the bed12
LC_ALL=C sort -t$'\t' -k1,1 -k2,2n -k3,3nr -k4,4 NONCODEv6_mm10.lncAndGene.bed > NONCODEv6_mm10.resorted.bed

#Step 1: split genes and transcripts

BED=NONCODEv6_mm10.resorted.bed     # any order works now
awk '$4~/^NONMMUG/' $BED | cut -f1-6 > genes6.bed
awk '$4~/^NONMMUT/' $BED | cut -f1-6 > tx6.bed

#Step 2: find each transcript's gene (the parent)
bedtools intersect -s -f 1.0 -wa -wb -a tx6.bed -b genes6.bed \
 | awk 'BEGIN{FS=OFS="\t"}{l=$9-$8; if(!($4 in b)||l<bl[$4]||(l==bl[$4]&&$10<b[$4])){b[$4]=$10;bl[$4]=l}} END{for(t in b)print b[t],t}' > map_contained.tsv

#Step 3: orphan transcripts get their own gene

awk 'BEGIN{FS=OFS="\t"} NR==FNR{h[$2];next} !($4 in h){print $4"_gene",$4}' map_contained.tsv tx6.bed > map_orphan_tx.tsv

#Step 4: orphan genes become one-transcript genes

awk 'NR==FNR{h[$1];next} !($4 in h){print $4}' map_contained.tsv genes6.bed > orphan_genes.txt
awk '{print $1"\t"$1"_tx"}' orphan_genes.txt > map_orphan_genes.tsv

#Step 5: combine into one map and flip it for bed2gtf

cat map_contained.tsv map_orphan_tx.tsv map_orphan_genes.tsv > isoforms_final.tsv
awk 'BEGIN{FS=OFS="\t"}{print $2,$1}' isoforms_final.tsv > isoforms_final.flip.txt
cut -f1 isoforms_final.flip.txt | sort | uniq -d | wc -l      # expect 0

#Step 6: build the BED12 that bed2gtf will read
awk 'BEGIN{FS=OFS="\t"} NR==FNR{o[$1];next} $4~/^NONMMUT/ || ($4 in o){if($4 in o)$4=$4"_tx"; $7=$2;$8=$2;$9="0,0,0";print}' orphan_genes.txt $BED > T.nc.rgb.bed

#convert 
bed2gtf -i T.nc.rgb.bed -I isoforms_final.flip.txt -o NONCODEv6_mm10.approachA.gtf -T 8

#Step 8: merge the duplicate gene lines
F=NONCODEv6_mm10.approachA.gtf
awk -F'\t' -v OFS='\t' 'NR==FNR{if($3=="gene"){k=$9; if(!(k in s)||$4+0<s[k])s[k]=$4+0; if(!(k in e)||$5+0>e[k])e[k]=$5+0} next} $3=="gene"{if(d[$9]++)next; $4=s[$9]; $5=e[$9]} 1' $F $F > NONCODEv6_mm10.final.gtf

#Step 9: check the result

F=NONCODEv6_mm10.final.gtf
cut -f3 $F | sort | uniq -c
cut -f1 isoforms_final.tsv | sort -u | wc -l     # expected genes
wc -l < T.nc.rgb.bed                             # expected transcripts
awk '{s+=$10}END{print s}' T.nc.rgb.bed          # expected exons
wc -l < map_orphan_tx.tsv; wc -l < orphan_genes.txt   # orphan transcripts, orphan genes