#create bed of gtf using for loop 

for gtf in *.gtf; do
    awk -F'\t' 'BEGIN{OFS="\t"} $3=="gene" {
        match($9, /gene_id "([^"]+)"/, a)
        print $1, $4-1, $5, a[1]
    }' "$gtf" > "${gtf%.gtf}_gene_level.bed"
done