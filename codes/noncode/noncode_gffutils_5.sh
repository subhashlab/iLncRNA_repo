#!/bin/bash
#PBS -N noncode_gffutils
#PBS -q cpu-test
#PBS -l nodes=1:ppn=4
#PBS -l mem=16gb

exec > /tmp/noncode_gffutils_${PBS_JOBID}.out 2>&1
set -x

cd /iitjfs/home/2025pmd0062/ref/hg38/master_gtf/gffutils || { echo "cd failed"; exit 1; }  #replace file path with your own file path

/iitjfs/home/2025pmd0062/miniconda3/envs/myenv/bin/python - << 'PYEOF'
import gffutils
db = gffutils.create_db(
    "/iitjfs/home/2025pmd0062/ref/hg38/RNAcentragtf/all_gtf/NONCODEv6_human_hg38_lncRNA.gtf", #replace file path with your own file path
    dbfn="/iitjfs/home/2025pmd0062/ref/hg38/master_gtf/gffutils/noncode.db",  #replace file path with your own file path
    force=True,
    keep_order=True,
    merge_strategy="merge",
    disable_infer_genes=False,
    disable_infer_transcripts=False
)
with open("/iitjfs/home/2025pmd0062/ref/hg38/master_gtf/gffutils/NONCODE_hg38_with_genes.gtf", "w") as out:  #replace file path with your own file path
    for gene in db.features_of_type("gene", order_by=("seqid", "start")):
        out.write(str(gene) + "\n")
        for tx in db.children(gene, featuretype="transcript", order_by="start"):
            out.write(str(tx) + "\n")
            for exon in db.children(tx, featuretype="exon", order_by="start"):
                out.write(str(exon) + "\n")
print("Done.")
PYEOF

cp /tmp/noncode_gffutils_${PBS_JOBID}.out /iitjfs/home/2025pmd0062/ref/hg38/master_gtf/gffutils/  #replace file path with your own file path
