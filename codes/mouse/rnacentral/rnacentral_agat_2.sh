#!/bin/bash
#PBS -N agat_gff2gtf_rnacentral_mm10
#PBS -q cpu-test
#PBS -l select=1:ncpus=30:mem=64gb
#PBS -j oe
#PBS -o /iitjfs/home/2025pmd0062/master_gtf_2/mouse/rnacentral_mm10/agat_gff2gtf.log

# Activate conda env with AGAT
source ~/miniconda3/etc/profile.d/conda.sh
conda activate agat_env
hash -r

# Verify AGAT is available before running
which agat_convert_sp_gff2gtf.pl

# Run conversion
agat_convert_sp_gff2gtf.pl \
  --gff /iitjfs/home/2025pmd0062/master_gtf_2/mouse/rnacentral_mm10/mus_musculus.GRCm38_fixed_exon.gff3 \
  -o /iitjfs/home/2025pmd0062/master_gtf_2/mouse/rnacentral_mm10/mus_musculus.GRCm38_fixed_exon.gtf

echo "Done. Feature type summary:"
awk -F'\t' '!/^#/ {print $3}' /iitjfs/home/2025pmd0062/master_gtf_2/mouse/rnacentral_mm10/mus_musculus.GRCm38_fixed_exon.gtf | sort | uniq -c