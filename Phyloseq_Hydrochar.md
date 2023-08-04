# Phyloseq

## Quick install (BiocManager)

```
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("phyloseq")
```


## Data
The aim of this project was to characterise the microbial diversity of anaerobic reactors. The data is part of the project [Effect of hydrochar from acid hydrolysis on anaerobic digestion of chicken manure](https://www.sciencedirect.com/science/article/pii/S2213343722012167), where solid residues from the production of biofuel precursor levulinic acid were generated via microwave-assisted acid hydrolysis of Miscanthus x Giganteus and investigated for the supplementation of anaerobic reactors digesting Chicken Manure. 

## Load libraries

```
install.packages("dplyr")     # To manipulate dataframes
install.packages("readxl")    # To read Excel files into R
install.packages("ggplot2")   # for high quality graphics
source("https://bioconductor.org/biocLite.R")
biocLite("phyloseq")  
```

## Read the data and transform it to Phyloseq objects
Three different tables are needed (tables were obtained from QIIME2) and organised in Excel:
•	**OTU**
•	**Taxonomy**
•	**Samples**

<img width="1242" alt="Screenshot 2023-08-04 at 11 35 19" src="https://github.com/srz11d/Phyloseq_for_microbial_diversity/assets/135147161/827769cd-0f4f-44ed-9c34-83b76054756f">

_OTU Table_

<img width="1242" alt="Screenshot 2023-08-04 at 11 42 51" src="https://github.com/srz11d/Phyloseq_for_microbial_diversity/assets/135147161/92aef170-5352-433f-811d-2b39aefb1e9c">

_Taxonomy Table_

<img width="157" alt="Screenshot 2023-08-04 at 11 43 40" src="https://github.com/srz11d/Phyloseq_for_microbial_diversity/assets/135147161/5b0c5d6e-d355-4f5f-a57f-3fe404e2eb5b">

_Samples Table_

### Read data
```
setwd("/Users/user/Documents/Postdoc/Sequencing/Anaerobic_Digestion")
otu_mat<-read_excel("/Users/user/Documents/Postdoc/Sequencing/Anaerobic_Digestion/OTU_table_trd.xlsx", sheet = "OTU matrix")
tax_mat<-read_excel("/Users/user/Documents/Postdoc/Sequencing/Anaerobic_Digestion/OTU_taxonomy_trd.xlsx", sheet = "Taxonomy table")
samples_df<-read_excel("/Users/user/Documents/Postdoc/Sequencing/Anaerobic_Digestion/table_samples_trd.xlsx", sheet = "Samples")
```

## Define the row names
```
otu_mat <- otu_mat %>%
  tibble::column_to_rownames("out") 
tax_mat <- tax_mat %>% 
  tibble::column_to_rownames("out")
samples_df <- samples_df %>% 
  tibble::column_to_rownames("sample")
```

## Transform into matrixes otu to tables (sample table can be left as data frame)
```
otu_mat <- as.matrix(otu_mat)
tax_mat <- as.matrix(tax_mat)
```

## Transform to Phyloseq objects
```
OTU = otu_table(otu_mat, taxa_are_rows = TRUE)
TAX = tax_table(tax_mat)
samples = sample_data(samples_df)
carbom <- phyloseq(OTU, TAX, samples)
carbom
```
## Visualise data
```
sample_names(carbom)
rank_names(carbom)
sample_variables(carbom)
```

## Normalize the number of reads in each sample using median sequencing depth
```
total = median(sample_sums(carbom))
standf = function(x, t=total) round(t * (x / sum(x)))
carbom = transform_sample_counts(carbom, standf)
```

## Bar graph
The bar graph shown here was created at the phylum level, but depending on the data, graphs can be created at different levels

```
plot_bar(carbom, fill = "Phylum")
```
<img width="972" alt="Screenshot 2023-08-04 at 11 53 23" src="https://github.com/srz11d/Phyloseq_for_microbial_diversity/assets/135147161/4ca593c2-be10-430b-841e-d72e12e08f24">

_Relative abundance at phylum level_


## Heatmap
The heatmap also shows the (abundance) of items, however, for this specific case the map contained multiple rows and the information is not clear. 

It is possible to edit the input table, but for this specific case the heatmap was not longer considered relevant.

```
plot_heatmap(carbom, method = "NMDS", distance = "bray")
```

<img width="547" alt="Screenshot 2023-08-04 at 12 01 38" src="https://github.com/srz11d/Phyloseq_for_microbial_diversity/assets/135147161/da5bb294-ac1e-46ce-8c52-40f6b942cb5c">

## Alpha diversity
Alpha diversity measurements can be calculated individually or all together with the function **Plot_richness** 
```
plot_richness(carbom, color = "samples", measures=c("Chao1", "Shannon", "Observed", "InvSimpson"))
plot_richness(carbom, color = "samples")
estimate_richness(carbom, measures = "Observed")
```
<img width="543" alt="Screenshot 2023-08-04 at 12 07 37" src="https://github.com/srz11d/Phyloseq_for_microbial_diversity/assets/135147161/3d2417b8-230f-4d5e-bff2-6e9d47d74e84">

<img width="524" alt="Screenshot 2023-08-04 at 12 08 38" src="https://github.com/srz11d/Phyloseq_for_microbial_diversity/assets/135147161/2bc4b9c3-dca1-4bd5-bc01-ee5d5c0a9afe">

<img width="531" alt="Screenshot 2023-08-04 at 12 09 16" src="https://github.com/srz11d/Phyloseq_for_microbial_diversity/assets/135147161/024f94dc-aec5-4b13-952d-2575e434de5b">

## Ordination
Multivariate analysis is based on Bray-Curtis distance and NMDS ordination.

```
carbom.ord <- ordinate(carbom, "NMDS", "bray")
plot_ordination(carbom, carbom.ord, type = "taxa", color = "Class", shape = "Domain", title = "OTUs")
plot_ordination(carbom, carbom.ord, type = "taxa", color = "Phylum", shape = "Domain", title = "OTUs")
plot_ordination(carbom, carbom.ord, type = "taxa", color = "Class", title = "OTUs", label = "Class") + facet_wrap(~Domain, 3)
```
<img width="531" alt="Screenshot 2023-08-04 at 12 11 27" src="https://github.com/srz11d/Phyloseq_for_microbial_diversity/assets/135147161/11d4615d-8677-4212-9f81-b0bd0c718e20">



