# Phyloseq

## Quick install (BiocManager)

```
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("phyloseq")
```

## Load libraries

```
install.packages("dplyr")     # To manipulate dataframes
install.packages("readxl")    # To read Excel files into R
install.packages("ggplot2")   # for high quality graphics
source("https://bioconductor.org/biocLite.R")
biocLite("phyloseq")  
```

## Data
The aim of this project was to characterise the microbial diversity of anaerobic reactors. The data is part of the project [Effect of hydrochar from acid hydrolysis on anaerobic digestion of chicken manure](https://www.sciencedirect.com/science/article/pii/S2213343722012167), where solid residues from the production of biofuel precursor levulinic acid were generated via microwave-assisted acid hydrolysis of Miscanthus x Giganteus and investigated for the supplementation of anaerobic reactors digesting Chicken Manure. 


