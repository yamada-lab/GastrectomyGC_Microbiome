library("phyloseq")
library("ggplot2")
library("plyr")
library("vegan")
library("ggplot2")

relab <-read.table("<insert the tables for species, modules, or metabolites here: e.g ./Raw_Data/Species_mOTUAnn_5.1e-05.up.tsv>", header = T, row.names = 1, sep = "\t")
sample_id <- read.delim("<sample metadata files here: e.g ./Raw_Data/Metadata_Table.tsv",sep="\t",header=TRUE,row.names = 1,check.names = FALSE)

#converting to phyloseq format
relabP = otu_table(relab, taxa_are_rows = TRUE)
sampledata = sample_data(sample_id)
physeq1 = phyloseq(relabP,sampledata)

ordu = ordinate(physeq1, "PCoA", "bray")
ordplot=plot_ordination(physeq1, ordu,type="sample",color="<columns for viz>",axes = c(1,2,3)) 
#cols <- c("Gastrectomy" = "#DC9238", "Healthy" = "#0DA4DA")
ordplot + 
  stat_ellipse(type = "t") +
  theme_bw()#+scale_colour_manual(values = cols)
  
  
#specify the '<>' part with the path to the file or parameter of interest