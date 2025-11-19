# Code accompanying Skrutl et al. 

## Installation

Clone this repository:

```bash
git clone https://github.com/MicrobiologyETHZ/skrutl_satellite_deletion.git
cd skrutl_satellite_deletion
```

## DE Analysis

To reproduce the differential expression analysis:

1. Pull the Bioconductor Docker image:
   ```bash
   docker pull bioconductor/bioconductor_docker:latest
   ```

2. Run the Docker container:
   ```bash
   docker run -d -v ./skrutl_satellite_deletion:/home/rstudio -e PASSWORD=yourpassword -p 8787:8787 bioconductor/bioconductor_docker:latest
   ```

3. Install required R packages (in a new terminal):
   ```bash
   docker exec bioconductor/bioconductor_docker:latest Rscript /home/rstudio/DE_analysis/install_packages.R
   ```

4. Run the DE analysis script:
   ```bash
   docker exec bioconductor/bioconductor_docker:latest Rscript /home/rstudio/DE_analysis/22-02-23-deleteion-deseq.R
   ```

Alternatively, you can access RStudio Server at `http://localhost:8787` (username: `rstudio`, password: `yourpassword`) and run the scripts interactively.

  
## Phylogenetics

The phylogenetic analysis constructs a maximum likelihood tree and performs population structure analysis using SNP and deletion datasets from the Drosophila GDL.

**Key analyses:**
- SNP filtering, LD pruning, and phylogenetic tree construction using IQ-TREE
- Principal component analysis (PCA) for both SNP and deletion datasets
- Population structure visualization

**Software requirements:** vcftools, bcftools, plink2, iqtree2, Python 3 (Biopython, pandas, plotly)

For detailed step-by-step instructions, see [phylogenetics/README.md](phylogenetics/README.md).
