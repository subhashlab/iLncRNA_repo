#!/bin/bash
#PBS -N agat_convert2
#PBS -q cpu-test
#PBS -o agat_convert2.log
#PBS -e agat_error2.err
#PBS -l select=1:ncpus=30:mem=64gb
#PBS -j oe

source ~/miniconda3/etc/profile.d/conda.sh
conda activate agat_env

agat_convert_sp_gff2gtf.pl \
  --gff /iitjfs/home/2025pmd0062/master_gtf_2/mouse/rnacentral_mm10/mus_musculus.GRCm38.gff3 \      #path to input file
  --out /iitjfs/home/2025pmd0062/master_gtf_2/mouse/rnacentral_mm10/mus_musculus.GRCm38_fixed_exon.gff3        #path to output file
