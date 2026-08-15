#!/bin/bash
#PBS -N agat_convert2
#PBS -q cpu-test
#PBS -l nodes=1:ppn=4
#PBS -l mem=64gb

exec > /tmp/agat_convert2_${PBS_JOBID}.out 2>&1
set -x

rm -rf /iitjfs/home/2025pmd0062/ref/hg38/RNAcentragtf/all_gff3/agat_convert2/agat_tmp_homo_sapiens.GRCh38_fixed/
rm -f /iitjfs/home/2025pmd0062/ref/hg38/RNAcentragtf/all_gff3/agat_convert2/homo_sapiens.GRCh38_fixed.gtf

env -i \
    PERL5LIB=/iitjfs/home/2025pmd0062/miniconda3/envs/agat_env/lib/perl5/site_perl/x86_64-linux-thread-multi:/iitjfs/home/2025pmd0062/miniconda3/envs/agat_env/lib/perl5/site_perl \
    PATH=/iitjfs/home/2025pmd0062/miniconda3/envs/agat_env/bin:/usr/bin:/bin \
    /iitjfs/home/2025pmd0062/miniconda3/envs/agat_env/bin/agat_convert_sp_gff2gtf.pl \
    --gff /iitjfs/home/2025pmd0062/ref/hg38/RNAcentragtf/all_gff3/agat_convert2/homo_sapiens.GRCh38_fixed_v2.gff3 \
    --out /iitjfs/home/2025pmd0062/ref/hg38/RNAcentragtf/all_gff3/agat_convert2/homo_sapiens.GRCh38_fixed_v2.gtf

cp /tmp/agat_convert2_${PBS_JOBID}.out /iitjfs/home/2025pmd0062/ref/hg38/RNAcentragtf/all_gff3/agat_convert2/