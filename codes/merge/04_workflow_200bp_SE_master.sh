# Step 1: Count exons per transcript_id
# Input : master_lincRNA_nonversioned.gtf
# Output: exon_counts.tsv   (columns: transcript_id, exon_count)


GTF="master_lincRNA_nonversioned.gtf"
OUT="exon_counts.tsv"

gawk -F'\t' '
    $3=="exon" {
        match($9, /transcript_id "([^"]+)"/, a)
        cnt[a[1]]++
    }
    END { for (t in cnt) print t"\t"cnt[t] }
' "$GTF" > "$OUT"

echo "Wrote $OUT ($(wc -l < "$OUT") transcripts with exons)"


# Step 2: Compute MATURE SPLICED transcript length
#   spliced_length = sum of (exon_end - exon_start + 1) over all exons
#                    belonging to that transcript_id
#   (NOT the genomic span of the "transcript" feature line -- span
#    includes introns and can look >=200bp even when the actual mature
#    transcript that reads get counted against is under 200bp. This is
#    why featureCounts matrices were showing genes with <200bp effective
#    length even after filtering on genomic span.)
# Input : master_lincRNA_nonversioned.gtf
# Output: transcript_lengths.tsv   (columns: transcript_id, spliced_len_bp)


GTF="master_lincRNA_nonversioned.gtf"
OUT="transcript_lengths.tsv"

gawk -F'\t' '
    $3=="exon" {
        match($9, /transcript_id "([^"]+)"/, a)
        len[a[1]] += ($5 - $4 + 1)
    }
    END { for (t in len) print t"\t"len[t] }
' "$GTF" > "$OUT"

echo "Wrote $OUT ($(wc -l < "$OUT") transcripts, spliced length)"



# Step 3: Build failing-transcript list
#   Fails if: single-exon (exon_count == 1) OR spliced_length < MINLEN
# Inputs : exon_counts.tsv, transcript_lengths.tsv (now spliced length)
# Output : failing_transcripts.txt
#
# NOTE: transcript_lengths.tsv is now built only from exon lines
# (step2), so a transcript with zero exon lines (shouldn't exist in a
# valid GTF, but worth a sanity check) would not appear in this file
# and therefore would NOT be flagged here. Cross-check against the
# "transcript" feature count if that's a concern:
#   comm -23 <(awk -F'\t' '$3=="transcript"{match($9,/transcript_id "([^"]+)"/,a);print a[1]}' master_lincRNA_nonversioned.gtf | sort -u) \
#            <(cut -f1 transcript_lengths.tsv | sort -u)


EXON_COUNTS="exon_counts.tsv"
TRANSCRIPT_LENGTHS="transcript_lengths.tsv"
OUT="failing_transcripts.txt"
MINLEN=200

gawk -v minlen="$MINLEN" -F'\t' '
    NR==FNR { exoncount[$1]=$2; next }
    {
        tid=$1; len=$2
        ec = (tid in exoncount) ? exoncount[tid] : 0
        if (ec == 1 || len < minlen) print tid
    }
' "$EXON_COUNTS" "$TRANSCRIPT_LENGTHS" | sort -u > "$OUT"

N_FAIL=$(wc -l < "$OUT")
N_TOTAL=$(wc -l < "$TRANSCRIPT_LENGTHS")
echo "Wrote $OUT: ${N_FAIL} / ${N_TOTAL} transcripts flagged for removal"


# Step 4: Identify genes that would become empty
#   A gene fails only if ALL of its transcripts are in failing_transcripts.txt
# Inputs : failing_transcripts.txt, master_lincRNA_nonversioned.gtf
# Output : failing_genes.txt


GTF="master_lincRNA_nonversioned.gtf"
FAILING_T="failing_transcripts.txt"
OUT="failing_genes.txt"

gawk -F'\t' '
    FNR==NR { fail[$1]=1; next }
    $3=="transcript" {
        match($9, /gene_id "([^"]+)"/, g)
        match($9, /transcript_id "([^"]+)"/, t)
        all_genes[g[1]]=1
        if (!(t[1] in fail)) survive_genes[g[1]]=1
    }
    END {
        for (gid in all_genes) if (!(gid in survive_genes)) print gid
    }
' "$FAILING_T" "$GTF" | sort -u > "$OUT"

echo "Wrote $OUT ($(wc -l < "$OUT") genes will be dropped)"


# Step 5: Write master_lincRNA_hg38.gtf
#   - Drop gene lines whose gene_id is in failing_genes.txt
#   - Drop any non-gene line (transcript/exon/CDS/UTR/start_codon/
#     stop_codon/...) whose transcript_id is in failing_transcripts.txt
#   - Single pass over the original file, so line order is preserved
#     and shared exon lines (distinct lines per transcript_id in GTF)
#     are kept whenever their own transcript survives.
# Inputs : master_lincRNA_nonversioned.gtf, failing_transcripts.txt, failing_genes.txt
# Output : master_lincRNA_hg38.gtf


GTF="master_lincRNA_nonversioned.gtf"
FAILING_T="failing_transcripts.txt"
FAILING_G="failing_genes.txt"
OUT="master_lincRNA_hg38.gtf"

gawk -F'\t' -v OFS='\t' -v failT="$FAILING_T" -v failG="$FAILING_G" '
BEGIN {
    while ((getline line < failT) > 0) fail_t[line]=1
    close(failT)
    while ((getline line < failG) > 0) fail_g[line]=1
    close(failG)
}
{
    gid=""; tid=""
    if (match($9, /gene_id "([^"]+)"/, g)) gid=g[1]
    if (match($9, /transcript_id "([^"]+)"/, t)) tid=t[1]

    if ($3 == "gene") {
        if (!(gid in fail_g)) print
        next
    }
    if (tid=="" || !(tid in fail_t)) print
}
' "$GTF" > "$OUT"

echo "Wrote $OUT"


# Step 6: Validate master_lincRNA_hg38.gtf
#   - before/after feature counts
#   - orphan child check (child feature whose transcript_id doesn't
#     exist as a transcript line in the output)
#   - empty gene check (gene line with no surviving transcript)
#   - spliced-length check: confirm NO surviving transcript has
#     mature spliced length < 200bp (this is the check that matters
#     for the featureCounts <200bp issue you were seeing)
# All three checks should print nothing except the final PASS line.


GTF="master_lincRNA_nonversioned.gtf"
OUT="master_lincRNA_hg38.gtf"
MINLEN=200

echo "--- feature counts ---"
for feat in gene transcript exon CDS; do
    before=$(awk -F'\t' -v f="$feat" '$3==f' "$GTF" | wc -l)
    after=$(awk -F'\t' -v f="$feat" '$3==f' "$OUT" | wc -l)
    printf "%-10s before=%-8s after=%-8s removed=%s\n" "$feat" "$before" "$after" "$((before-after))"
done

echo "--- orphan child check (should print nothing) ---"
gawk -F'\t' '
    NR==FNR {
        if ($3=="transcript") { match($9,/transcript_id "([^"]+)"/,t); tset[t[1]]=1 }
        next
    }
    $3!="gene" {
        match($9,/transcript_id "([^"]+)"/,t)
        if (!(t[1] in tset)) print "ORPHAN CHILD:", $0
    }
' "$OUT" "$OUT"

echo "--- empty gene check (should print nothing) ---"
gawk -F'\t' '
    NR==FNR {
        if ($3=="transcript") { match($9,/gene_id "([^"]+)"/,g); gset[g[1]]=1 }
        next
    }
    $3=="gene" {
        match($9,/gene_id "([^"]+)"/,g)
        if (!(g[1] in gset)) print "EMPTY GENE:", $0
    }
' "$OUT" "$OUT"

echo "--- spliced-length check (should print nothing) ---"
gawk -v minlen="$MINLEN" -F'\t' '
    $3=="exon" {
        match($9, /transcript_id "([^"]+)"/, a)
        len[a[1]] += ($5 - $4 + 1)
    }
    END {
        for (t in len) if (len[t] < minlen) print "SHORT SURVIVOR:", t, len[t]"bp"
    }
' "$OUT"

echo "Validation complete."

#gene       before=26202    after=17654    removed=8548
#transcript before=68671    after=58925    removed=9746
#exon       before=198134   after=188377   removed=9757
#CDS        before=0        after=0        removed=0


#gene       before=26158    after=17610    removed=8548
#transcript before=68579    after=58832    removed=9747
#exon       before=197946   after=188187   removed=9759

#gene       before=37727    after=27444    removed=10283
#transcript before=143789   after=101990   removed=41799
#exon       before=360654   after=318018   removed=42636

#12/9/26
#gene       before=26111    after=17561    removed=8550
#transcript before=68146    after=58452    removed=9694
#exon       before=196577   after=186873   removed=9704

