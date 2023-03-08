# to install phyloseq
#if (!require("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")

#BiocManager::install("phyloseq")

library(permute)
library(lattice)
library(ape)
library(vegan)
library(dplyr)
library(phyloseq)
library(ggplot2)


#Load mapping file into R
meta <- read.table("mapping_File.txt",header=TRUE,row.names=1,sep="\t",stringsAsFactors=FALSE)
sampleData <- sample_data(meta)

#Load reformatted FUNGuild data into R
FG <- read.table("FUNGuild_example_fix_header.txt",header=T,sep="\t",row.names=1)

#Select only the column that we need
FGotus <- select(FG, -(Taxonomy:Citation.Source))
FGotumat <- as(as.matrix(FGotus), "matrix")
FGOTU <- otu_table(FGotumat, taxa_are_rows = TRUE)
FGtaxmat <- select(FG, Confidence.Ranking, Trophic.Mode, Guild, Growth.Morphology)
FGtaxmat <- as(as.matrix(FGtaxmat),"matrix")
FGTAX = tax_table(FGtaxmat)

#Creating phyloseq object
physeq = phyloseq(FGOTU,FGTAX,sampleData)
physeq.prune = prune_taxa(taxa_sums(physeq) > 1, physeq)
physeq.prune.nopossible = subset_taxa(physeq.prune, Confidence.Ranking="Highly Possible")
physeq.prune.nopossible = subset_taxa(physeq.prune.nopossible, Confidence.Ranking!="-")

#Create color palette
cbbPalette <- c("#009E73","#999999", "#E69F00", "#56B4E9", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "midnightblue", "lightgreen","saddlebrown", "brown", "aquamarine4","lavenderblush2","snow3", "darkblue", "darkgoldenrod1", "darkseagreen", "darkorchid", "darkolivegreen1", "black","lightskyblue", "darkgreen", "deeppink", "khaki2", "firebrick", "brown1", "darkorange1", "cyan1", "royalblue4", "darksalmon", "darkblue","royalblue4", "dodgerblue3", "steelblue1", "lightskyblue", "darkseagreen", "darkgoldenrod1", "darkseagreen", "darkorchid", "darkolivegreen1", "brown1", "darkorange1", "cyan1", "darkgrey")

#Create a plot
FUNGuildcom = ggplot(data = psmelt(physeq.prune.nopossible), mapping = aes_string(x = "Sample" ,y = "Abundance", fill = "Trophic.Mode" )) + geom_bar(stat="identity", position="fill") + ggtitle("Fungal Trophic Mode Composition ")+theme(axis.text.x = element_text(angle = 90, hjust = 1)) + scale_fill_manual(values = cbbPalette)

#Save a plot to PDF file
pdf("./Fungal Trophic Mode Composition.pdf", width = 8, height = 5)
FUNGuildcom
dev.off()
pdf 
