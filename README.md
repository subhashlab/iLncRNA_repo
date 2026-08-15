# Master GTF — Comprehensive intergenic lncRNA (lincRNA) Annotation Pipeline (hg38)

A multi-database lincRNA annotation pipeline that integrates reference annotations from GENCODE, NONCODE and RNAcentral into unified master GTF files for the human genome assembly hg38.

---
Sources:   
GENCODE v49 :   https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/gencode.v49.primary_assembly.annotation.gtf.gz  
NONCODE v6  :   https://v7.noncode.org/datadownload/NONCODEv6_human_hg38_lncRNA.gtf.gz  
RNA Central v26:   https://ftp.ebi.ac.uk/pub/databases/RNAcentral/releases/26.0/genome_coordinates/gff3/homo_sapiens.GRCh38.gff3.gz  

---
Tools applied to complete gene, transcript, exon hierarchy in NONCODEv6_human_hg38_lncRNA.gtf : gffutils (v0.14)  
Tools applied for the conversion of homo_sapiens.GRCh38.gff3 format into GTF format using batch job: AGAT tool (v1.7.0)  
Tools applied to remove overlapping gene entries between lncRNA and non-lncRNA: bedtools (v2.31.1)  
## Overview

This pipeline consolidates lincRNA annotations from three major databases into one master GTF file:

| Output file | Description |
|---|---|
| `master_lincRNA_1.gtf` | Merged GTF (GENCODE + NONCODE + RNAcentral) |

Priority order: GENCODE_v49 >NONCODE_v6 >RNACentral_v26
---

## Database Sources

| Database | Version | Format | Notes |
|---|---|---|---|
| GENCODE | v49 | GTF | Primary reference; lncRNA biotype used |
| NONCODE | v6 | GTF | lncRNA gtf with transcript and exons entries |
| RNAcentral | v26 | GFF3 → GTF | Converted via AGAT|

---

## Repository Structure

```
iLncRNA_repo/
├── GTF/
│   └── master_lincRNA_1.gtf
├── codes/
│   ├── gencode/
│   │   ├── README.md
│   │   └── gencode_lincRNA.sh
│   ├── merge/
│   │   ├── merge.sh
│   │   └── remove_chr_patch.sh
│   ├── noncode/
│   │   └── modif/
│   │       ├── noncode_modif.sh
│   │       ├── noncode_gffutils_5.sh
│   │       └── noncode_lincRNA.sh
│   ├── rnacentral/
│   │   └── modif/
│   │       ├── changes_rc_v2.sh
│   │       ├── agat_convert3.sh
│   │       └── rnacentral_lincRNA.sh
│   └── verify_overlap/
│       ├── 00_create_bed.sh
│       ├── 01_unique_entries.sh
│       └── 02_verify_unique_enties_bed.sh
└── README.md
```
---

## Pipeline Steps

### 1. Format Conversion (RNAcentral)

RNAcentral is distributed as GFF3 and must be converted to GTF using [AGAT](https://github.com/NBISweden/AGAT):

```bash
agat_convert_sp_gff2gtf.pl \
  --gff homo_sapiens.GRCh38.gff3 \      #path to input file
  --out homo_sapiens.GRCh38.gtf         #path to output file
```

> HPC note: Run on the Agastya cluster (`cpu-test` queue via PBS/qsub). Use full binary paths rather than `conda activate` on compute nodes. AGAT progress is written to stderr; PBS logs are written only after job completion.

### 2. lncRNA Extraction

Extract lncRNA features from each source GTF

### 3. nonlncRNA extraction
Extract lncRNA features from each source GTF

### 4. Remove overlaps using bedtools (Separate environment for bedtools and agat)

### 5. Generate curated_lincRNA.bed (must be bed4) 

bed4 has chr,start and end which are used to pull exactly required entries from respective gtf 

### 6. Generate GTF from bed4 
Pull gene, transcript and exon entries from GTF using bed as reference. 

### 7. Harmonize curated GTF using GENCODE v49 as reference. Modify col9 with attributes:
For gene level:
gene_id, gene_type, gene_name, level, tag, havana_gene

For transcript level: 
gene_id, transcript_id, gene_type, gene_name, transcript_type, transcript_name,
level, transcript_support_level (tsl), tag, havana_transcript, havana_gene

For exon level: 
gene_id, transcript_id, gene_type, gene_name, transcript_type, transcript_name,
exon_number, exon_id, level, tag, havana_transcript, havana_gene

# Merge harmonized GTFs into master_lincRNA GTF 
Outputs:
- `master_lincRNA_1.gtf`
---

## Dependencies

| Tool | Version | Purpose |
|---|---|---|
| AGAT | ≥ 1.0 | GFF3 → GTF conversion |
| bedtools | ≥ 2.30 | Interval intersection |
| awk | system | Coordinate manipulation |

---

## Known Issues / Caveats

- **Coordinate system:** BED files use 0-based half-open coordinates; GTF files use 1-based closed coordinates. All conversions add +1 to BED start positions before GTF output.

---

## Reference Genome

All annotations are based on **hg38 (GRCh38)**.

---
