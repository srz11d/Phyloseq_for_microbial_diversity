library("phyloseq")
library("ggplot2")      # graphics
library("readxl")       # necessary to import the data from Excel file
library("dplyr")        # filter and reformat data frames
library("tibble")       # Needed for converting column to row names
library("vegan")
library("DESeq2")
setwd("/Users/user/Documents/Postdoc/Sequencing/Anaerobic_Digestion")
otu_mat<-read_excel("/Users/user/Documents/Postdoc/Sequencing/Anaerobic_Digestion/OTU_table_trd.xlsx", sheet = "OTU matrix")
tax_mat<-read_excel("/Users/user/Documents/Postdoc/Sequencing/Anaerobic_Digestion/OTU_taxonomy_trd.xlsx", sheet = "Taxonomy table")
samples_df<-read_excel("/Users/user/Documents/Postdoc/Sequencing/Anaerobic_Digestion/table_samples_trd.xlsx", sheet = "Samples")
#Define the row names#
otu_mat <- otu_mat %>%
  tibble::column_to_rownames("out") 
tax_mat <- tax_mat %>% 
  tibble::column_to_rownames("out")
samples_df <- samples_df %>% 
  tibble::column_to_rownames("sample") 
#Transform into matrixes otu tx tables (sample table can be left as data frame)#
otu_mat <- as.matrix(otu_mat)
tax_mat <- as.matrix(tax_mat)
#Transform to phyloseq objects#
OTU = otu_table(otu_mat, taxa_are_rows = TRUE)
TAX = tax_table(tax_mat)
samples = sample_data(samples_df)
carbom <- phyloseq(OTU, TAX, samples)
carbom
#Visualise data#
sample_names(carbom)
rank_names(carbom)
sample_variables(carbom)
#Normalize number of reads in each sample using median sequencing depth#
total = median(sample_sums(carbom))
standf = function(x, t=total) round(t * (x / sum(x)))
carbom = transform_sample_counts(carbom, standf)
#Bar graph#
plot_taxa_bar(carbom, fill = "Phylum") #graph needs to be chanhed#
plot_bar(carbom, fill = "Phylum")
#Heatmap#
plot_heatmap(carbom, method = "NMDS", distance = "bray") #table needs to be changed#
#Heatmap with only the 20%#
#carbom_abun <- filter_taxa(carbom, function(x) sum(x > total*0.20) > 1, TRUE)
#carbom_abun
carbom_abun = transform_sample_counts(carbom, function(x) x/ sum(x))
carbom_abund = filter_taxa(carbom_abun, function(x) sum(x) > 0.05, TRUE)
otu_table(carbom_abun)[1:8, 1:5]
plot_heatmap(carbom_abund, method = "NMDS", distance = "bray")
plot_heatmap(carbom_abund, method = "NMDS", taxa.label = "Phylum", taxa.order = "Phylum")
plot_heatmap(carbom_abund, method = "NMDS", taxa.label = "Genus", taxa.order = "Genus")
#bar
plot_bar(carbom_abund, fill = "Genus")
plot_bar(carbom_abund, fill = "Phylum")
plot_bar(domain, fill = "Domain")
#Alpha diversity#
plot_richness(carbom, color = "samples", measures=c("Chao1", "Shannon", "Observed", "InvSimpson"))
plot_richness(carbom, color = "samples")
estimate_richness(carbom, measures = "Observed")
#Ordination
carbom.ord <- ordinate(carbom, "NMDS", "bray")
plot_ordination(carbom, carbom.ord, type = "taxa", color = "Class", shape = "Domain", title = "OTUs")
plot_ordination(carbom, carbom.ord, type = "taxa", color = "Phylum", shape = "Domain", title = "OTUs")
plot_ordination(carbom, carbom.ord, type = "taxa", color = "Class", title = "OTUs", label = "Class") + facet_wrap(~Domain, 3)

plot_ordination(carbom, ordinate(carbom, method='RDA'), type = "split", color = "sample", shape='Domain')




