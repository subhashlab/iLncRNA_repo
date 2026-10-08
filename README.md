# Master GTF — Comprehensive intergenic lncRNA (lincRNA) Annotation Pipeline (hg38)

A multi-database lincRNA annotation pipeline that integrates reference annotations from **GENCODE**, **NONCODE**, and **RNAcentral** into unified master GTF file for the human genome assembly **hg38 (GRCh38)**.

## Human annotation

| iLncRNA version	| iLncRNA release date | Organism	| Genome assembly version	| Corresponding GENCODE release	| Corresponding Ensembl release |
|-------|-----------------|---------------|---------|----------------|----------------|
|[v1.1](GTF/human/hs_master_iLncRNAs_annotation_v1.1.gtf)| 08/11/2026 | *Homo sapiens*	|	GRCh38.p14 |	49	|115|

## Mouse annotation

| iLncRNA version	| iLncRNA release date | Organism	| Genome assembly version	| Corresponding GENCODE release	| Corresponding Ensembl release |
|-------|-----------------|---------------|---------|----------------|----------------|
|[v1.0](GTF/human/mouse/mm_master_iLncRNA_annotation_v1.0.gtf)| 09/10/2026 | *Mus musculus*	|	GRCm38.p6 |	23	|98|
---

## Table of Contents

1. [Project Objectives](#1-project-objectives)
2. [Current Project Status](#2-current-project-status)
3. [Dataset and Source Information](#3-dataset-and-source-information)
4. [Repository Structure](#4-repository-structure)
5. [Description of Major Files and Subdirectories](#5-description-of-major-files-and-subdirectories)
6. [Analysis Workflow](#6-analysis-workflow)
7. [Software and Tool Versions](#7-software-and-tool-versions)
8. [Important Parameters and Conventions](#8-important-parameters-and-conventions)
9. [Relevant GitHub Repositories](#9-relevant-github-repositories)
10. [Reproducing the Analysis](#10-reproducing-the-analysis)
11. [Known Issues and Caveats](#11-known-issues-and-caveats)

---

## 1. Project Objectives

The goal of this project is to produce a **comprehensive, non-redundant set of intergenic lncRNA (lincRNA) annotations** for the human genome (hg38) by consolidating entries from three major public databases:

-  **GENCODE v49**, **NONCODE v6**, and **RNAcentral v26**
- Resolve format inconsistencies (GFF3 → GTF conversion, incomplete gene/transcript/exon hierarchies)
- Remove overlaps between lncRNA and non-lncRNA loci using interval arithmetic
- Harmonize GTF attribute fields (column 9) to match GENCODE v49 conventions
- Merge harmonized annotations into a single **master lincRNA GTF**

The resulting resource is intended for use as a reference annotation in RNA-seq quantification, differential expression analysis, and lincRNA functional studies.

---

## 2. Current Project Status

| Component | Status |
|---|---|
| GENCODE v49 lincRNA extraction | ✅ Complete |
| NONCODE v6 hierarchy repair (gffutils) | ✅ Complete |
| RNAcentral v26 GFF3 → GTF conversion (AGAT) | ✅ Complete |
| lincRNA extraction from all three sources | ✅ Complete |
| Non-lncRNA overlap removal (bedtools) | ✅ Complete |
| BED generation and GTF reconstruction | ✅ Complete |
| GTF attribute harmonization (GENCODE v49 schema) | ✅ Complete |
| Master GTF merge (master_lincRNA_1.1.gtf) | ✅ Complete |

---

## 3. Dataset and Source Information

### Reference Genome

All annotations are mapped to **hg38 / GRCh38** (primary assembly).

### Source Databases

| Database | Version | Format | Download URL |
|---|---|---|---|
| GENCODE | v49 | GTF | https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/gencode.v49.primary_assembly.annotation.gtf.gz |
| NONCODE | v6 | GTF | https://v7.noncode.org/datadownload/NONCODEv6_human_hg38_lncRNA.gtf.gz |
| RNAcentral | v26 | GFF3 | https://ftp.ebi.ac.uk/pub/databases/RNAcentral/releases/26.0/genome_coordinates/gff3/homo_sapiens.GRCh38.gff3.gz |

### Priority Order for Merging

When gene loci overlap across databases, annotations are retained in the following priority order:

```
GENCODE v49  >  NONCODE v6  >  RNAcentral v26
```

### Output Files

| Output File | Description |
|---|---|
| master_lincRNA_1.1.gtf | Merged, harmonized lincRNA GTF (GENCODE + NONCODE + RNAcentral) |

---

## 4. Repository Structure

```
iLncRNA_repo/
├── GTF/
│   └── master_lincRNA_1.1.gtf          # Final merged output
├── codes/
│   ├── gencode/
│   │   ├── README.md                  # GENCODE-specific notes
│   │   └── gencode_lincRNA.sh         # lincRNA extraction from GENCODE v49
│   ├── merge/
│   │   ├── merge.sh                   # Merge harmonized GTFs into master
│   │   └── remove_chr_patch.sh        # Remove alt/patch chromosome entries
│   ├── noncode/
│   │   ├── noncode_modif.sh # General NONCODE modifications
│   │   └── modif/
│   │       ├── noncode_gffutils_5.sh  # Repair gene/transcript/exon hierarchy
│   │       └── noncode_lincRNA.sh     # lincRNA extraction from NONCODE v6
│   ├── rnacentral/
│   │    ├── agat_convert3.sh       # GFF3 → GTF conversion via AGAT (PBS job)
│   │    ├── rnacentral_lincRNA.sh  # lincRNA extraction from RNAcentral v26
│   │    └── modif/
│   │       ├── changes_rc_v2.sh       # General RNA Central modifications
│   └── verify_overlap/
│       ├── 00_create_bed.sh           # Convert curated GTF from each database to BED4 format
│       ├── 01_unique_entries.sh       # Identify non-redundant loci
│       └── 02_verify_unique_enties_bed.sh  # Validate uniqueness of BED entries
└── README.md
```

---

## 5. Description of Major Files and Subdirectories

### GTF/
Contains the final output GTF file(s). master_lincRNA_1.1.gtf is the primary deliverable — a merged, harmonized annotation of intergenic lncRNAs from all three sources.

### codes/gencode/
Scripts for extracting lncRNA biotype entries from the GENCODE v49 primary assembly annotation. GENCODE serves as the **reference schema** for attribute harmonization.

### codes/noncode/modif/
Script for maintaining the NONCODE v6 curated GTF as per GENCODE, The raw NONCODE GTF contains entries but lacks complete col9 attributes. noncode_gffutils_5.sh uses gffutils to infer and insert missing gene entries, producing a complete gene → transcript → exon hierarchy.

### codes/rnacentral/modif/
Scripts for converting RNAcentral v26 GFF3 to GTF format (via AGAT) and subsequent cleanup. agat_convert3.sh is a PBS batch job submitted to the Agastya HPC cluster.

### codes/merge/
Scripts for the final merge step. remove_chr_patch.sh strips non-primary assembly contigs (alternate loci, patches) before merging. merge.sh concatenates the harmonized per-source GTFs into`master_lincRNA_1.1.gtf.

### codes/verify_overlap/
A three-step QC workflow to confirm that gene loci in the final GTF are non-redundant. Converts GTF entries to BED4, computes unique intervals, and verifies them against the source GTFs.

---

## 6. Analysis Workflow

The pipeline proceeds in eight ordered steps:

### Step 1 — Format Conversion (RNAcentral GFF3 → GTF)

RNAcentral is distributed as GFF3 and must be converted to GTF using [AGAT](https://github.com/NBISweden/AGAT):

```bash
agat_convert_sp_gff2gtf.pl \
  --gff homo_sapiens.GRCh38.gff3 \
  --out homo_sapiens.GRCh38.gtf
```

> **HPC note:** Run on the Agastya cluster (`cpu-test` queue via PBS/qsub). Use full binary paths rather than `conda activate` on compute nodes. AGAT writes progress to stderr; PBS logs are written only after job completion.

### Step 2 — Hierarchy Repair (NONCODE)

The NONCODE v6 GTF lacks gene-level records. Use gffutils to reconstruct the complete gene → transcript → exon hierarchy:

```bash
bash codes/noncode/modif/noncode_gffutils_5.sh
```

### Step 3 — lincRNA Extraction

Extract entries with lincRNA (intergenic lncRNA) biotype from each source GTF:

```bash
bash codes/gencode/gencode_lincRNA.sh
bash codes/noncode/modif/noncode_lincRNA.sh
bash codes/rnacentral/modif/rnacentral_lincRNA.sh
```

### Step 4 — Non-lncRNA Overlap Removal (bedtools)

Remove lncRNA loci that overlap non-lncRNA gene regions using bedtools intersect. This produces a 4-column BED file (chr, start, end, name) AGAT and bedtools require separate conda environments:

```bash
# Run in bedtools environment
bedtools intersect -v -a lncRNA_candidates.bed -b nonlncRNA.bed > nonoverlap.bed
```

> **Coordinate note:** BED uses 0-based half-open coordinates. GTF uses 1-based closed coordinates. All BED → GTF conversions add **+1 to BED start positions**.

### Step 5 — GTF Reconstruction from BED4

Pull gene, transcript, and exon entries from each source GTF using the curated BED4 as a coordinate filter.

### Step 6 — GTF Attribute Harmonization

Standardize column 9 attributes across all three sources to match the **GENCODE v49 schema**:

| Feature Level | Required Attributes |
|---|---|
| gene | gene_id, gene_type, gene_name, level, tag, havana_gene |
| transcript | gene_id, transcript_id, gene_type, gene_name, transcript_type, transcript_name, level, transcript_support_level, tag, havana_transcript, havana_gene |
| exon | gene_id, transcript_id, gene_type, gene_name, transcript_type, transcript_name, exon_number, exon_id, level, tag, havana_transcript, havana_gene |

### Step 7 — Merge into Master GTF

Concatenate harmonized per-source GTFs and remove alt/patch contigs to produce the final output:

```bash
bash codes/merge/merge.sh
bash codes/merge/remove_chr_patch.sh
# Output: GTF/master_lincRNA_1.gtf
```

---

## 7. Software and Tool Versions

| Tool | Version | Purpose | Environment |
|---|---|---|---|
| AGAT | v1.7.0 | GFF3 → GTF conversion | Separate conda env |
| bedtools | v2.31.1 | Interval intersection and overlap removal | Separate conda env |
| gffutils | v0.14 | Gene/transcript/exon hierarchy repair (NONCODE) | Python |
| awk | system | Coordinate manipulation and attribute parsing | Shell |
| Python | ≥ 3.8 | gffutils scripting | — |

> **Important:** AGAT and bedtools must be run in **separate conda environments** due to dependency conflicts.

---

## 8. Important Parameters and Conventions

| Parameter / Convention | Value / Rule |
|---|---|
| Reference genome assembly | hg38 / GRCh38 (primary assembly only) |
| Coordinate system — BED | 0-based, half-open [start, end) |
| Coordinate system — GTF | 1-based, closed [start, end] |
| BED → GTF start correction | GTF_start = BED_start + 1 |
| GTF field delimiter | **Tab only** (spaces in column 9 will break parsers) |
| Merge priority | GENCODE v49 > NONCODE v6 > RNAcentral v26 |
| Biotype filter | lncRNA (lncRNA only) |
| Chromosomes retained | Primary assembly chromosomes only (chr1–22, chrX, chrY, chrM); alt/patch contigs removed |

---

## 9. Relevant GitHub Repositories

| Repository | Description |
|---|---|
| [NBISweden/AGAT](https://github.com/NBISweden/AGAT) | GFF3/GTF conversion and manipulation toolkit |
| [arq5x/bedtools2](https://github.com/arq5x/bedtools2) | Genome arithmetic (overlap, intersect, subtract) |
| [daler/gffutils](https://github.com/daler/gffutils) | GFF/GTF parsing and database construction in Python |

---

## 10. Reproducing the Analysis

### Prerequisites

1. Clone this repository:
   ```bash
   git clone https://github.com/teenu2207/iLncRNA_repo.git
   cd iLncRNA_repo
   ```

2. Download source annotation files (see [Section 3](#3-dataset-and-source-information)) into your working reference directory.

3. Set up conda environments:
   ```bash
   # Environment 1: AGAT
   conda create -n agat_env -c bioconda agat=1.7.0
   
   # Environment 2: bedtools
   conda create -n bedtools_env -c bioconda bedtools=2.31.1
   
   # Environment 3: gffutils (Python)
   pip install gffutils==0.14
   ```

4. If running on an HPC cluster (PBS/qsub), use full binary paths inside job scripts rather than conda activate.

### Execution Order

Run scripts in the following order:

```
Step 1:  codes/rnacentral/modif/agat_convert3.sh       # GFF3 → GTF (AGAT env)
Step 2:  codes/noncode/modif/noncode_gffutils_5.sh     # Hierarchy repair
Step 3:  (i)codes/gencode/gencode_lincRNA.sh           # Extract GENCODE lincRNA
         (ii.a)codes/noncode/noncode_lincRNA.sh        # Extract NONCODE lincRNA
         (ii.b)codes/noncode/modif/noncode_modif.sh    # Modify col9 attributes
         (iii.a)codes/rnacentral/rnacentral_lincRNA.sh # Extract RNAcentral lincRNA
         (iii.b) codes/rnacentral/modif/rnacentral_modif.sh #Modify col9 attributes
Step 4:  (i)codes/verify_overlap/00_create_bed.sh         # Generate BED4
         (ii)codes/verify_overlap/01_unique_entries.sh     
         (iii)codes/verify_overlap/02_unique_enties_bed.sh
Step 5:  codes/merge/merge.sh                          # Final merge → master_lincRNA.gtf
Step 6:  codes/merge/remove_chr_patch.sh               # Remove alt contigs → master_lincRNA_1.1.gtf
         
```

### Validation

After merge, verify non-redundancy of the output:

```bash
bash codes/verify_overlap/02_verify_unique_enties_bed.sh
```

---

## 11. Known Issues and Caveats

| Issue | Details |
|---|---|
| **Coordinate system mismatch** | BED files use 0-based coordinates; GTF uses 1-based. All BED → GTF conversions must add +1 to start positions. |
| **NONCODE incomplete hierarchy** | Raw NONCODE v6 GTF lacks gene-level records; gffutils repair is required before any downstream processing. |
| **AGAT/bedtools environment conflict** | AGAT and bedtools cannot share a conda environment; always activate the correct environment before running each step. |
| **GTF tab delimiter strictness** | Column 9 attributes must be tab-delimited at the field level. Space delimiters will cause silent parsing failures in downstream tools. |
| **Alt/patch contigs** | RNAcentral and NONCODE may include entries on alternate loci or patch chromosomes; remove_chr_patch.sh filters these before the final merge. |

---

## Citation / Acknowledgements

- **GENCODE:** Mudge, J. M., Carbonell-Sala, S., Diekhans, M., Martinez, J. G., Hunt, T., Jungreis, I., Loveland, J. E., Arnan, C., Barnes, I., Bennett, R., Berry, A., Bignell, A., Cerdán-Vélez, D., Cochran, K., Cortés, L. T., Davidson, C., Donaldson, S., Dursun, C., Fatima, R., Hardy, M., … Frankish, A. (2025). GENCODE 2025: reference gene annotation for human and mouse. Nucleic acids research, 53(D1), D966–D975. https://doi.org/10.1093/nar/gkae1078
- **NONCODE:** Zhao, L., Wang, J., Li, Y., Song, T., Wu, Y., Fang, S., Bu, D., Li, H., Sun, L., Pei, D., Zheng, Y., Huang, J., Xu, M., Chen, R., Zhao, Y., & He, S. (2021). NONCODEV6: an updated database dedicated to long non-coding RNA annotation in both animals and plants. Nucleic acids research, 49(D1), D165–D171. https://doi.org/10.1093/nar/gkaa1046  
- **RNAcentral:** The RNAcentral Consortium , RNAcentral in 2026: genes and literature integration, Nucleic Acids Research, Volume 54, Issue D1, 6 January 2026, Pages D303–D313, https://doi.org/10.1093/nar/gkaf1329
- **AGAT:** Jacques Dainat, Robrecht Cannoodt, André Soares, Daniel García Ruano, Darío Hereñú, Dr. K. D. Murray, Ed Davis, Ivan Ugrin, Kathryn Crouch, Lucile Soler, pascal-git, Zachary Zollman& tayyrov. (2026). NBISweden/AGAT: AGAT v1.7.0 (Version v1.7.0) [Computer software]. Zenodo. https://doi.org/10.5281/zenodo.19499560
- **gffutils:** 	https://github.com/daler/gffutils

---

*Maintained by the Subhash lab, IIT Jammu. For questions or issues, please open a GitHub Issue.*
