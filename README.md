# Population Structure Analysis with Whole Genome Sequencing Data

## Overview
This script provides a global analysis of population structure using Whole Genome Sequencing (WGS) data. The main analyses performed include:

- **Principal Component Analysis (PCA)** across all chromosomes
- **F-statistics (FST) computation** for population differentiation
- **sNMF analysis** for ancestry inference
- **Site Frequency Spectrum (SFS) computation**
- **Diversity analysis (DIV)**

## Authors
- Elise Gay (EPHE)
- Romuald Laso-Jadart (EPHE)
- Stefano Mona (EPHE)

**Please inform the authors before sharing this script.**

---

## Dependencies
This script requires the following R libraries:

```r
library(LEA)
library(vcfR)
library(rlist)
library(PopGenome)
library(qvalue)
library(pegas)
library(ggplot2)
library(adegenet)
library(hierfstat)
library(withr)
library(ggrepel)
library(reshape2)
library(gridExtra)
library(pcadapt)```

# Population Structure Analysis with Whole Genome Sequencing Data

## Overview
This script provides a global analysis of population structure using Whole Genome Sequencing (WGS) data. The main analyses performed include:

- **Principal Component Analysis (PCA)** across all chromosomes
- **F-statistics (FST) computation** for population differentiation
- **sNMF analysis** for ancestry inference
- **Site Frequency Spectrum (SFS) computation**
- **Diversity analysis (DIV)**

## Authors
- Stefano Mona
- Elise Gay
- Romuald Laso-Jadart
- Pierre Lesturgie (EPHE - MNHN)
- 2023

**Please inform the authors before sharing this script.**

---

## Dependencies
This script requires the following R libraries:

```r
library(LEA)
library(vcfR)
library(rlist)
library(PopGenome)
library(qvalue)
library(pegas)
library(ggplot2)
library(adegenet)
library(hierfstat)
library(withr)
library(ggrepel)
library(reshape2)
library(gridExtra)
library(pcadapt)
```

Ensure these packages are installed before running the script.

---

## Input Data
### 1. Population Metadata
A metadata table containing sample IDs and their associated population labels is required.

**Example format (metadata/Samples_table.txt):**
```plaintext
samples  pop  
sample_1 pop1 
sample_2 pop1 
sample_3 pop1 
sample_4 pop2 
sample_5 pop2 
```

### 2. VCF Files
Filtered VCF files containing genotype data should be provided for analysis.

---

## Usage
### 1. Load Metadata
```r
table_pop=read.table("metadata/Samples_table.txt", header = TRUE, row.names = 1)
samples=row.names(table_pop)
pop=table_pop$pop
```

### 2. Compute Pairwise FST
```r
data_FST=read.vcfR("data/global_FST/super_1/WGS_21_SUPER_1_TAG_Filtered.vcf")
res_FST=calcola_fst_pairwise_bootstrap(data_FST, lista_pop, 0.05, 10)
```

### 3. Perform Global FST Analysis
```r
vcf_data=readData("data/global_FST/super_1/", format = "VCF")
res=boot_popgenome(vcf_data, lista_pop, bootstrap=10)
```

### 4. Conduct PCA
```r
respca_10 = pcadapt(input = read.pcadapt("data/PCA_FST/WGS_21_SUPER_1_TAG_Filtered.vcf", type = "vcf"), min.maf=0.05, K = 10)
scores = data.frame(respca_10$scores)
```

### 5. Run sNMF for Ancestry Analysis
```r
setwd("data/snmf/")
project = snmf("WGS_21_SUPER_1_TAG_Filtered_Plink.ped", K=1:8, entropy=TRUE, repetitions = 20, project = "new")
```

### 6. Compute Site Frequency Spectrum (SFS)
```r
sfs_result = site.spectrum(vcf_data)
```

---

## Output
- **FST Values**: Pairwise and global FST results for population differentiation.
- **PCA Plot**: Visualization of population clustering.
- **sNMF Results**: Admixture proportions for individuals.
- **SFS Plot**: Visualization of allele frequency distribution.

---

## Notes
- Ensure all input files are correctly formatted before execution.
- Adjust parameters such as `bootstrap` and `K` as needed for robust analysis.
- Results are generated in the corresponding analysis folders.

---

## License
This script is for research purposes only. Please credit the authors when using or modifying the script.

