# Phylogenetic and Population Structure Analysis

This directory contains reproducible workflows for SNP-based phylogenetic tree construction and population structure analysis of GDL (Global Diversity Lines) using both SNP and deletion datasets.

## Table of Contents
- [Software Requirements](#software-requirements)
- [Input Data](#input-data)
- [Analysis Workflows](#analysis-workflows)
  - [SNP Analysis Pipeline](#1-snp-analysis-pipeline)
  - [Deletion Analysis Pipeline](#2-deletion-analysis-pipeline)
- [Scripts and Notebooks](#scripts-and-notebooks)
- [References](#references)

## Software Requirements

- vcftools v0.1.16 ([Danecek et al. 2011](https://doi.org/10.1093/bioinformatics/btr330))
- bcftools v1.21 ([Danecek et al. 2021](https://doi.org/10.1093/gigascience/giab008))
- plink2 v2.0.0 ([Chang et al. 2015](https://doi.org/10.1186/s13742-015-0047-8))
- iqtree v2.4.0 ([Nguyen et al. 2015](https://doi.org/10.1093/molbev/msu300))
- Python 3 with Biopython, pandas, numpy, plotly, scipy

## Input Data

### SNP Dataset
- **Source**: [Lack et al. 2016](https://doi.org/10.1093/gbe/evz022)
- **File**: `autosomes_no4_IBD_masked_biallelic_ZW184removed_Callability_masked_biallelic_neutral_sites.vcf.gz`
- **Description**: SNPs limited to autosomes (excluding chromosome 4). Only small intronic (positions 32–65 bp) and 4-fold degenerate positions were used based on genomic annotations generated using SNPeff (167857 SNPs before filtering)
- **Availability**: See original publication.

### Deletion Dataset
- **Source**: [Lack et al. 2016](https://doi.org/10.1093/gbe/evz022)
- **File**: `GDL_Indels.vcf.gz`
- **Description**: Insertion/deletion variants from GDL lines
- **Availability**: NCBI SRA [SRP050151](https://trace.ncbi.nlm.nih.gov/Traces/index.html?view=study&acc=SRP050151)

## Analysis Workflows

### 1. SNP Analysis Pipeline

#### Step 1.1: Filter SNPs
Filter VCF to remove indels, low-quality variants, and rare alleles:

```bash
vcftools --gzvcf autosomes_no4_IBD_masked_biallelic_ZW184removed_Callability_masked_biallelic_neutral_sites.vcf.gz \
  --remove-indels \
  --mac 1 \
  --minQ 30 \
  --min-meanDP 5 \
  --max-meanDP 100 \
  --minDP 5 \
  --maxDP 100 \
  --max-missing 0.9 \
  --recode \
  --stdout | gzip > dme_filtered.vcf.gz
```

**Output**: Filtered VCF with high-quality SNPs

#### Step 1.2: Add unique variant IDs
```bash
bcftools annotate --set-id '%CHROM\_%POS\_%REF\_%ALT' dme_filtered.vcf.gz -o dme_filtered_uniqid.vcf
```

#### Step 1.3: Convert to PLINK format
```bash
plink2 --vcf dme_filtered_uniqid.vcf \
  --double-id \
  --allow-extra-chr \
  --make-bed \
  --out dme-prune
```

#### Step 1.4: LD pruning
Prune linked SNPs using a sliding window approach:

```bash
plink2 --bfile dme-prune \
  --allow-extra-chr \
  --double-id \
  --indep-pairwise 50 10 0.2 \
  --out dme-prune
```

**Parameters**: Window size=50 SNPs, step=10 SNPs, r²=0.2

#### Step 1.5: Create LD-pruned dataset and run PCA
```bash
plink2 --bfile dme-prune \
  --extract dme-prune.prune.in \
  --make-bed \
  --out dme-ldpruned \
  --allow-extra-chr \
  --pca
```

**Output**:
- `dme-ldpruned.eigenvec` - Principal component scores
- `dme-ldpruned.eigenval` - Eigenvalues (variance explained)
- Final dataset: 50,034 SNPs

#### Step 1.6: Convert to VCF for phylogenetic analysis
```bash
plink --bfile dme-ldpruned \
  --recode vcf \
  --out dme-ldpruned \
  --allow-extra-chr
```

#### Step 1.7: Convert VCF to PHYLIP format
```bash
python vcf2phylip.py -i dme-ldpruned.vcf
```

**Script**: `vcf2phylip.py` (included in this directory)
**Output**: `dme-ldpruned.min4.phy`

#### Step 1.8: Remove invariant positions
After LD pruning, some sites may become invariant. Remove these using the notebook:

```bash
jupyter notebook 2025-11-10_remove-invariant-positions-from-phylip.ipynb
```

**Script**: `2025-11-10_remove-invariant-positions-from-phylip.ipynb`
**Function**: Identifies and removes invariant sites from PHYLIP alignment
**Output**: `data/dme-ldpruned.min4.filtered-ivariant.phy`


#### Step 1.9: Build phylogenetic tree
```bash
iqtree2 -s data/dme-ldpruned.min4.filtered-ivariant.phy \
  -m MFP+ASC \
  -B 1000 \
  -T 16
```

**Parameters**:
- `-m MFP+ASC`: Model finder with ascertainment bias correction
- `-B 1000`: 1000 bootstrap replicates
- `-T 16`: 16 threads
- **Best model selected**: SYM+ASC+R6

**Output**: Maximum likelihood tree with bootstrap support

#### Step 1.10: Visualize PCA results
```bash
jupyter notebook 2025-11-10_SNP-analysis.ipynb
```

**Script**: `2025-11-10_SNP-analysis.ipynb`
**Input**: `data/dme-ldpruned.eigenvec`, `data/dme-ldpruned.eigenval`
**Output**: Interactive PCA plots, pairwise distance calculations


### 2. Deletion Analysis Pipeline


#### Step 2.1: Filter deletions
```bash
mkdir -p indel_out

vcftools --gzvcf GDL_Indels.vcf.gz \
  --minQ 30 \
  --minDP 5 \
  --maxDP 100 \
  --min-meanDP 5 \
  --max-meanDP 100 \
  --maf 0.05 \
  --max-missing 0.9 \
  --recode \
  --out indel_out/gdl_indels_filtered
```

**Output**: High-quality deletion variants

#### Step 2.3: Split multiallelic variants
```bash
bcftools norm -m -any \
  indel_out/gdl_indels_filtered.recode.vcf \
  -o indel_out/gdl_indels_filtered_split.vcf
```

#### Step 2.4: Convert to PLINK and exclude chromosome X
```bash
plink2 --vcf indel_out/gdl_indels_filtered_split.vcf \
  --allow-extra-chr \
  --make-bed \
  --out indel_out/gdl_indels_filtered_split \
  --threads 8 \
  --not-chr X
```

**Output**: 66,385 autosomal insertions/deletions

#### Step 2.5: Run PCA
```bash
plink2 --bfile indel_out/gdl_indels_filtered_split \
  --pca 10 \
  --out data/gdl_indels_filtered_pca \
  --allow-extra-chr
```

**Output**:
- `data/gdl_indels_filtered_pca.eigenvec` - PC scores
- `data/gdl_indels_filtered_pca.eigenval` - Variance explained

#### Step 2.6: Visualize results
```bash
jupyter notebook 2025-11-10_SNP-analysis.ipynb
```

**Script**: `2025-11-10_SNP-analysis.ipynb` (same notebook for both SNP and deletion PCA)
**Input**: `data/gdl_indels_filtered_pca.eigenvec`, `data/gdl_indels_filtered_pca.eigenval`
**Output**: Interactive 3D PCA plots, pairwise distances

## Scripts and Notebooks

### `vcf2phylip.py`
Converts VCF files to PHYLIP, FASTA, or NEXUS format for phylogenetic analysis.
**Author**: Edgardo M. Ortiz (v2.9)
**Usage**: `python vcf2phylip.py -i <input.vcf>`

### `2025-11-10_remove-invariant-positions-from-phylip.ipynb`
Removes invariant sites from PHYLIP alignments that may arise after filtering/pruning.
**Input**: PHYLIP alignment with potential invariant sites
**Output**: Filtered PHYLIP alignment with only variable sites

### `2025-11-10_SNP-analysis.ipynb`
Performs PCA visualization and pairwise distance calculations for both SNP and deletion datasets.
**Functions**:
- Loads PCA results from PLINK
- Creates interactive scatter plots (2D for SNPs, 3D for deletions)
- Calculates Euclidean distances between samples in PC space
- Colors samples by geographic origin



## References

- Danecek et al. (2011) The variant call format and VCFtools. *Bioinformatics* 27(15):2156-2158. https://doi.org/10.1093/bioinformatics/btr330
- Danecek et al. (2021) Twelve years of SAMtools and BCFtools. *GigaScience* 10(2):giab008. https://doi.org/10.1093/gigascience/giab008
- Chang et al. (2015) Second-generation PLINK: rising to the challenge of larger and richer datasets. *GigaScience* 4:7. https://doi.org/10.1186/s13742-015-0047-8
- Nguyen et al. (2015) IQ-TREE: a fast and effective stochastic algorithm for estimating maximum-likelihood phylogenies. *Mol Biol Evol* 32(1):268-274. https://doi.org/10.1093/molbev/msu300
- Lack et al. (2016) A thousand fly genomes: an expanded Drosophila genome nexus. *Genome Biol Evol* 8(3):875-889. https://doi.org/10.1093/gbe/evz022
- Letunic & Bork (2024) Interactive Tree of Life (iTOL) v6: recent updates to the phylogenetic tree display and annotation tool. *Nucleic Acids Res* 52(W1):W78-W82. https://doi.org/10.1093/nar/gkae268
