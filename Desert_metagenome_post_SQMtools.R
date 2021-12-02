##### If necessary
# install.packages("wesanderson")
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("pathview", force = TRUE)


library(wesanderson)    # Beautiful color palettes
library(zCompositions)  # cmultRepl
library(CoDaSeq)        # log-centred ratios
library(vegan)          # PCA
library(scales)         # Plotting transparency
library(eulerr)         # This is the best package I've found for generating Venn diagrams
library(tidyverse)      # Always
library(reticulate) # Bring in Python
# reticulate::use_condaenv("SqueezeMeta") # Configure which version of Python to use
library(SQMtools)       # Required to read in the SqueezeMeta database

# Phylogenetic tree packages
library(ggimage)
# library(BiocManager)            # If necessary
# BiocManager::install("ggtree")  # If necessary
library(ggtree)
library(picante)
library(phytools)
library(phylotools)
library(ape)
library(adephylo)
library(phylobase)
# library(phylosignal) # Needed for phylo4d plotting
library(treeio)
library (devtools)
source_url("https://raw.githubusercontent.com/lrjoshi/FastaTabular/master/fasta_and_tabular.R")
library(seqinr)
library(phyloseq)
library(gtools)

rm(list = ls())

load("~/Dropbox/Cyclotron Research/Carbon_Oxidation_and_Fixation_Study/Metagenome/SteveM/Curated_Rcode/SQM_tools_desert.RData")

## Read in the dataframes
## The data represent:
# 1. df0 = the full dataset when first read in using SQMtools
# 2. df1 = filtered to just archaea, bacteria, and fungi
# 3. df2 = >50% complete and <10% contaminated
# 4. df3 = >2 reads
# 5. df4 = >2 occurrences
# filter.df = occurrences of KEGG and taxa throughout the filters
# mega.sub = the megahit ids for contigs were randomly selected as there multiple ids for the same contig taxonomy.
filter.df <- read.csv("~/Dropbox/Cyclotron Research/Carbon_Oxidation_and_Fixation_Study/Metagenome/SteveM/Curated_Rcode/filter.barplot.csv", row.names = 1)
kegg.df0 <- read.csv("~/Dropbox/Cyclotron Research/Carbon_Oxidation_and_Fixation_Study/Metagenome/SteveM/Curated_Rcode/kegg.full.csv", row.names = 1)
kegg.df1 <- read.csv("~/Dropbox/Cyclotron Research/Carbon_Oxidation_and_Fixation_Study/Metagenome/SteveM/Curated_Rcode/kegg.filt0.csv", row.names = 1)
kegg.df2 <- read.csv("~/Dropbox/Cyclotron Research/Carbon_Oxidation_and_Fixation_Study/Metagenome/SteveM/Curated_Rcode/kegg.filt1.csv", row.names = 1)
kegg.df3 <- read.csv("~/Dropbox/Cyclotron Research/Carbon_Oxidation_and_Fixation_Study/Metagenome/SteveM/Curated_Rcode/kegg.filt2.csv", row.names = 1)
kegg.df4 <- read.csv("~/Dropbox/Cyclotron Research/Carbon_Oxidation_and_Fixation_Study/Metagenome/SteveM/Curated_Rcode/kegg.filt3.csv", row.names = 1)
kegg.tpm.df0 <- read.csv("~/Dropbox/Cyclotron Research/Carbon_Oxidation_and_Fixation_Study/Metagenome/SteveM/Curated_Rcode/kegg.full.tpm.csv", row.names = 1)
kegg.tpm.df1 <- read.csv("~/Dropbox/Cyclotron Research/Carbon_Oxidation_and_Fixation_Study/Metagenome/SteveM/Curated_Rcode/kegg.filt0.tpm.csv", row.names = 1)
kegg.tpm.df2 <- read.csv("~/Dropbox/Cyclotron Research/Carbon_Oxidation_and_Fixation_Study/Metagenome/SteveM/Curated_Rcode/kegg.filt1.tpm.csv", row.names = 1)
kegg.tpm.df3 <- read.csv("~/Dropbox/Cyclotron Research/Carbon_Oxidation_and_Fixation_Study/Metagenome/SteveM/Curated_Rcode/kegg.filt2.tpm.csv", row.names = 1)
kegg.tpm.df4 <- read.csv("~/Dropbox/Cyclotron Research/Carbon_Oxidation_and_Fixation_Study/Metagenome/SteveM/Curated_Rcode/kegg.filt3.tpm.csv", row.names = 1)
taxa.df0 <- read.csv("~/Dropbox/Cyclotron Research/Carbon_Oxidation_and_Fixation_Study/Metagenome/SteveM/Curated_Rcode/taxa.full.csv", row.names = 1)
taxa.df1 <- read.csv("~/Dropbox/Cyclotron Research/Carbon_Oxidation_and_Fixation_Study/Metagenome/SteveM/Curated_Rcode/taxa.filt0.csv", row.names = 1)
taxa.df2 <- read.csv("~/Dropbox/Cyclotron Research/Carbon_Oxidation_and_Fixation_Study/Metagenome/SteveM/Curated_Rcode/taxa.filt1.csv", row.names = 1)#
taxa.df3 <- read.csv("~/Dropbox/Cyclotron Research/Carbon_Oxidation_and_Fixation_Study/Metagenome/SteveM/Curated_Rcode/taxa.filt2.csv", row.names = 1)
taxa.df4 <- read.csv("~/Dropbox/Cyclotron Research/Carbon_Oxidation_and_Fixation_Study/Metagenome/SteveM/Curated_Rcode/taxa.filt3.csv", row.names = 1)
mega.sub <- read.csv("~/Dropbox/Cyclotron Research/Carbon_Oxidation_and_Fixation_Study/Metagenome/SteveM/Curated_Rcode/mega.sub1.csv", row.names = 1)

##****************************
## Curating the taxonomy

## Contig taxonomy
contig.tax <- read.csv("~/Dropbox/Cyclotron Research/Carbon_Oxidation_and_Fixation_Study/Metagenome/SteveM/Curated_Rcode/taxonomy.full.csv",row.names = 1, header = T)
names(contig.tax) <- c("kingdom","phylum","class","order","family","genus","species")
contig.tax.sub0 <- contig.tax[rownames(contig.tax) %in% mega.sub$x,]
rownames(contig.tax.sub0)
rownames(full.sub5$contigs$tax)

## Remove any taxa we're not interested in
contig.tax.sub0 <- contig.tax.sub0[contig.tax.sub0$kingdom != "Unclassified",]            # 193
contig.tax.sub0 <- contig.tax.sub0[contig.tax.sub0$kingdom != "Viruses",]                 # 193
contig.tax.sub0 <- contig.tax.sub0[contig.tax.sub0$phylum != "Streptophyta",]             # 187
contig.tax.sub0 <- contig.tax.sub0[contig.tax.sub0$phylum != "Apicomplexa",]              # 186
contig.tax.sub0 <- contig.tax.sub0[contig.tax.sub0$phylum != "Chordata",]                 # 185
contig.tax.sub0 <- contig.tax.sub0[contig.tax.sub0$phylum != "Nematoda",]                 # 185
contig.tax.sub0 <- contig.tax.sub0[contig.tax.sub0$phylum != "Unclassified Eukaryota",]   # 184

## The taxonomy is coded weird in that taxa can be unclassified all the way through, but will be codified differently at finer taxonomies
## (e.g., "no genus in NCBI" versus "unclassified bacterium"). Let's group these together and recode as necessary
contig.tax.sub0$genus[grep("no genus in NCBI", contig.tax.sub0$genus)] <- paste("Unclassified", contig.tax.sub0$genus[grep("no genus in NCBI", contig.tax.sub0$genus)])
contig.tax.sub0$genus <- gsub(" bacterium \\(no genus in NCBI\\)","", contig.tax.sub0$genus)
contig.tax.sub0$genus <- gsub(" archaeon \\(no genus in NCBI\\)","", contig.tax.sub0$genus)

##****************************
## Relative abundances

# Extract a taxonomy matrix to use for subsetting later
taxa.df4.taxonomy <- as.data.frame(contig.tax.sub0)
rownames(taxa.df4.taxonomy) <- taxa.df4.taxonomy$species
# rownames(taxa.df4.taxonomy) <- make.contig.tax.sub0$species
levels(as.factor(rownames(taxa.df4.taxonomy)))

# Make sure there are no zero-sum taxa in any of the datasets from the various filters
taxa.df0 <- taxa.df0[rowSums(taxa.df0) != 0,]
taxa.df1 <- taxa.df1[rowSums(taxa.df1) != 0,]
taxa.df2 <- taxa.df2[rowSums(taxa.df2) != 0,]
taxa.df3 <- taxa.df3[rowSums(taxa.df3) != 0,]
taxa.df4 <- taxa.df4[rowSums(taxa.df4) != 0,]

# Create a relative abundance dataset from the full abundance table
taxa.rel.abu.df0 <- decostand(taxa.df0, MARGIN = 2, "total")*100
colSums(taxa.rel.abu.df0)       # Double check the samples sum to 100%

# Subset the relative abundance table to the final filtered dataset
taxa.rel.abu.df4 <- taxa.rel.abu.df0[rownames(taxa.df4.taxonomy),]
range(colSums(taxa.rel.abu.df4))       # This will show the summed relative abundances of the remaining taxa (27 to 75%)
identical(rownames(taxa.rel.abu.df4), rownames(taxa.df4.taxonomy))

# Convert to presence absence to calculate prevalence of each taxa
taxa.pa.df4a <- decostand(taxa.rel.abu.df4, "pa")
taxa.pa.df4b <- cbind.data.frame(apply(taxa.pa.df4a[,c(8,1,4,12)], 1, FUN = max),
                                 apply(taxa.pa.df4a[,c(9,2,10,11,3,7,5,6)], 1, FUN = max))  # Reorder the columns to group like sample types
names(taxa.pa.df4b) <- c("bulk","meth")
taxa.pa.df4c <- taxa.pa.df4b[rowSums(taxa.pa.df4b) < 3,]                      # Remove singletons/doubletons (there are none in this case)
taxa.pa.df4d <- specnumber(taxa.abs.abu.df4[c(8,1,4,12,9,2,10,11,3,7,5,6)], 
                           MARGIN = 2)                                        # Calculate taxa counts (richness)

##****************************
## Diversity analyses (richness and evenness)
taxa.abs.abu.df4 <- taxa.df0[rownames(taxa.df4.taxonomy),]                    # Sort the abundance df by taxonomy so the row names match
identical(names(taxa.pa.df4d), 
          names(taxa.abs.abu.df4[c(8,1,4,12,9,2,10,11,3,7,5,6)]))             # Check that the richness and abs abundance columns match
taxa.evenness <- diversity(taxa.abs.abu.df4[c(8,1,4,12,9,2,10,11,3,7,5,6)], 
                           MARGIN = 2, 
                           index="simpson")/
                           log(taxa.pa.df4d)                        # Calculate evenness
identical(names(taxa.pa.df4d), names(taxa.evenness))

# The ch4 and ch4+co2 are similar, so we'll treat them as replicates.
# Let's put a diversity figure in the supp to justify this step to reviewers.
names(taxa.pa.df4d) <- names(taxa.evenness) <- c("D621N","D734Y","G541N","G624Y","D621NM1","D621NM2","D734YM1","D734YM2","G541NM1","G541NM2","G624YM1","G624YM2")
taxa.pa.df4d1 <- taxa.pa.df4d[c("D621N","G541N","D734Y","G624Y","D621NM1","D621NM2","G541NM1","G541NM2","D734YM1","D734YM2","G624YM1","G624YM2")]
taxa.evenness1 <- taxa.evenness[c("D621N","G541N","D734Y","G624Y","D621NM1","D621NM2","G541NM1","G541NM2","D734YM1","D734YM2","G624YM1","G624YM2")]

# Plot the richness and evenness for each sample
jpeg("Taxa_diversity_samples.jpg", width = 5.46, height = 5.17, units = "in", res = 600)
par(mar = c(5.75,3.25,0.3,3.25))
plot(taxa.pa.df4d1, pch = c(16,16,17,17,16,16,16,16,17,17,17,17), col = rep(c("#35274A", "#0B775E"), c(4,8)), axes = F, xlab = "", ylab = "", ylim = range(taxa.pa.df4d))
axis(side = 2)
par(new = T)
plot(taxa.evenness1, pch = c(1,1,2,2,1,1,1,1,2,2,2,2), col = rep(c("#35274A", "#0B775E"), c(4,8)), axes = F, xlab = "", ylab = "", ylim = c(0.14,0.17))
axis(side = 1, at = c(1:12), labels = F)
text(x = 1:length(taxa.pa.df4d), y = 0.1373, labels = names(taxa.pa.df4d1), col = rep(c("#35274A", "#0B775E"), c(4,8)), xpd = NA, adj = 0, srt = -45)
text(x = c(2.5,8.5), y = 0.1315, labels = c("Bulk", expression(paste("CH"[4]))), col = c("#35274A", "#0B775E"), xpd = NA, adj = 0.25, srt = 0)
axis(side = 4)
mtext("Richness", side = 2, line = 2.2)
mtext("Evenness", side = 4, line = 2.2)
box()
dev.off()

colSums(taxa.pa.df4b != 0)

##****************************
## Generate files for phylogenetic tree (phylogram; might be near impossible for MG sequences)
## but also make a taxonomic tree (cladogram) based on contigs taxonomy at the species level

# Build a sample dataframe for the fully filtered data set
samp.df <- cbind.data.frame(sample.id = names(taxa.df4),
                            parent = substr(names(taxa.df4),1,1),
                            frostboil = substr(names(taxa.df4),2,4),
                            diapir.presence = substr(names(taxa.df4),5,5),
                            sample.type = substr(names(taxa.df4),6,6))
samp.df$sample.type <- gsub("^$|^ $","bulk", samp.df$sample.type)
samp.df$sample.type[samp.df$sample.type == "M" | samp.df$sample.type == "B"] <- "methane"
samp.df$parent[samp.df$parent == "D"] <- "dolomite"
samp.df$parent[samp.df$parent == "G"] <- "granite"
rownames(samp.df) <- samp.df$sample.id

# I blasted the long sequences that were causing issues with the phylogenetic tree and found
# there were more precise identifications on NCBI that contained shorter snippets of the 
# sequences that may play more nicely with the multiple sequence alignments. I'm going
# to recode those in the taxonomy and remake the tree files 
taxa.df4.taxonomy[grep("Arthrobacter", taxa.df4.taxonomy[,7]),]
taxa.df4.tax <- as.data.frame(taxa.df4.taxonomy)
taxa.df4.tax$species[grep(" bacterium", taxa.df4.tax$species)] <- paste("Unclassified ", gsub(" bacterium", "", taxa.df4.tax$species[grep(" bacterium", taxa.df4.tax$species)]), sep = "")
taxa.df4.tax$species[grep(" archaeon", taxa.df4.tax$species)] <- paste("Unclassified ", gsub(" archaeon", "", taxa.df4.tax$species[grep(" archaeon", taxa.df4.tax$species)]), sep = "")
taxa.df4.tax[grep("Oxalobacteraceae",taxa.df4.tax$species),] <- cbind.data.frame("Bacteria","Proteobacteria","Betaproteobacteria","Burkholderiales",	"Oxalobacteraceae",	"Massilia",	"Massilia violaceinigra B2")
taxa.df4.tax[grep("Propionibacteriaceae",taxa.df4.tax$species),] <- cbind.data.frame("Bacteria",	"Actinobacteria",	"Actinobacteria",	"Propionibacteriales",	"Propionibacteriaceae",	"Cutibacterium",	"Cutibacterium acnes")
taxa.df4.tax[grep("Lawsonella",taxa.df4.tax$species),] <- cbind.data.frame("Bacteria",	"Actinobacteria",	"Actinobacteria",	"Corynebacteriales",	"Lawsonellaceae",	"Lawsonella",	"Lawsonella clevelandensis")
taxa.df4.tax[grep("Propionibacteriaceae",taxa.df4.tax$species),] <- cbind.data.frame("Bacteria",	"Actinobacteria",	"Actinobacteria",	"Propionibacteriales",	"Propionibacteriaceae",	"Cutibacterium",	"Cutibacterium acnes")
taxa.df4.tax[grep("Burkholderiaceae",taxa.df4.tax$species),] <- cbind.data.frame("Bacteria",	"Proteobacteria",	"Betaproteobacteria",	"Burkholderiales",	"Burkholderiaceae",	"Paraburkholderia",	"Paraburkholderia caffeinilytica")
taxa.df4.tax[grep("Ralstonia",taxa.df4.tax$species),] <- cbind.data.frame("Bacteria",	"Proteobacteria",	"Betaproteobacteria",	"Burkholderiales",	"Burkholderiaceae",	"Ralstonia",	"Ralstonia pickettii")
# taxa.df4.tax[grep("epidermidis",taxa.df4.tax$species),] <- cbind.data.frame("Bacteria",	"Firmicutes",	"Bacilli",	"Bacillales",	"Staphylococcaceae",	"Staphylococcus",	"Staphylococcus epidermidis")
taxa.df4.tax[grep("Massilia",taxa.df4.tax$species),] <- cbind.data.frame("Bacteria",	"Proteobacteria",	"Betaproteobacteria",	"Burkholderiales",	"Oxalobacteraceae",	"Massilia",	"Massilia oculi")
taxa.df4.tax[grep("Nocardioidaceae",taxa.df4.tax$species),] <- cbind.data.frame("Bacteria",	"Actinobacteria",	"Actinobacteria",	"Propionibacteriales",	"Nocardioidaceae",	"Marmoricola",	"Marmoricola scoriae")
taxa.df4.tax[grep("Pseudomonas",taxa.df4.tax$species),] <- cbind.data.frame("Bacteria",	"Proteobacteria",	"Gammaproteobacteria",	"Pseudomonadales",	"Pseudomonadaceae",	"Pseudomonas",	"Pseudomonas yamanorum")
taxa.df4.tax[grep("Unclassified Staphylococcus",taxa.df4.tax$species),] <- cbind.data.frame("Bacteria",	"Firmicutes",	"Bacilli",	"Bacillales",	"Staphylococcaceae",	"Staphylococcus",	"Staphylococcus capitis")
taxa.df4.tax[grep("Burkholderia",taxa.df4.tax$species),] <- cbind.data.frame("Bacteria",	"Proteobacteria",	"Betaproteobacteria",	"Burkholderiales",	"Burkholderiaceae",	"Burkholderia",	"Burkholderia gladioli")
taxa.df4.tax[grep("Propionibacteriales",taxa.df4.tax$species),] <- cbind.data.frame("Bacteria",	"Actinobacteria",	"Actinomycetia",	"Propionibacteriales",	"Nocardioidaceae",	"Nocardioides",	"unclassified Nocardioides")
taxa.df4.tax[grep("Proteobacteria",taxa.df4.tax$species),] <- cbind.data.frame("Bacteria",	"Proteobacteria",	"Deltaproteobacteria",	"Myxococcales",	"Myxococcaceae",	"Corallococcus", "Corallococcus coralloides")
taxa.df4.tax[grep("Blastocatellia",taxa.df4.tax$species),] <- cbind.data.frame("Bacteria",	"Acidobacteria",	"Acidobacteriia",	"Acidobacteriales",	"Acidobacteriacea",	"Acidobacterium",	"Acidobacterium capsulatum")
taxa.df4.tax[grep("Marmoricola",taxa.df4.tax$species),] <- cbind.data.frame("Bacteria",	"Actinobacteria",	"Actinobacteria",	"Propionibacteriales",	"Nocardioidaceae",	"Marmoricola",	"Marmoricola scoriae")
# taxa.df4.tax[grep("Acidobacteria bacterium",taxa.df4.tax$species),] <- cbind.data.frame("Bacteria",	"Acidobacteria",	"Acidobacteriia",	"Bryobacterales",	"Unclassified Bryobacterales",	"Unclassified Bryobacterales",	"Unclassified Bryobacterales")
taxa.df4.tax[grep("Gemmatimonadales",taxa.df4.tax$species),] <- cbind.data.frame("Bacteria", "Actinobacteria", "Actinomycetia", "Corynebacteriales", "Nocardiaceae", "Nocardia","Nocardia farcinica")
taxa.df4.tax[grep("Actinobacteria",taxa.df4.tax$species),] <- cbind.data.frame("Bacteria", "Actinobacteria", "Actinomycetia", "Micrococcales", "Microbacteriaceae", "Protaetiibacter", "Unclassified Protaetiibacter")
taxa.df4.tax$species[taxa.df4.tax$species == "Arthrobacter sp."] <- "Unclassified Arthrobacter"
taxa.df4.tax$species[taxa.df4.tax$species == "Chryseobacterium sp."] <- "Unclassified Chryseobacterium"
# taxa.df4.tax$species[taxa.df4.tax$species == "Massilia sp."] <- "Unclassified Massilia"
taxa.df4.tax$species[taxa.df4.tax$species == "Microbacterium sp."] <- "Unclassified Microbacterium"
taxa.df4.tax$species[taxa.df4.tax$species == "Nocardioides sp."] <- "Unclassified Nocardioides"
taxa.df4.tax$species[taxa.df4.tax$species == "Nostoc sp."] <- "Unclassified Nostoc"
taxa.df4.tax$species[taxa.df4.tax$species == "Rhodococcus sp."] <- "Unclassified Rhodococcus"
taxa.df4.tax$species[taxa.df4.tax$species == "Streptococcus sp."] <- "Unclassified Streptococcus"
taxa.df4.tax$species[taxa.df4.tax$species == "Noviherbaspirillum sp."] <- "Unclassified Noviherbaspirillum"
taxa.df4.tax$species[taxa.df4.tax$species == "Propionibacterium sp."] <- "Unclassified Propionibacterium"
taxa.df4.tax$species[taxa.df4.tax$species == "Micromonospora sp."] <- "Unclassified Micromonospora"
taxa.df4.tax$species[taxa.df4.tax$species == "Aeromicrobium sp."] <- "Unclassified Aeromicrobium"

# Need unique names for the taxonomy because of contig multiplicity
taxa.df4.taxonomy <- as.matrix(taxa.df4.tax)                          # Convert to matrix
rownames(taxa.df4.taxonomy) <- make.unique(taxa.df4.taxonomy[,7])     # Make names unique (this appends a number to the name)
taxa.df4.taxonomy.df <- as.data.frame(taxa.df4.taxonomy)              # Convert to data frame
cbind(rownames(taxa.rel.abu.df4), rownames(taxa.df4.taxonomy.df))     # Check the row names. taxa.df4.taxonomy.df has the finer taxonomic classification
rownames(taxa.rel.abu.df4) <- rownames(taxa.df4.taxonomy.df)          # Assign rownames(taxa.df4.taxonomy.df) to rownames(taxa.rel.abu.df4) so they match
identical(rownames(taxa.rel.abu.df4), rownames(taxa.df4.taxonomy.df)) # Check that they match

# Use phyloseq to subset the abundance, taxonomy, and sample dfs at the same time
desert.phy <- phyloseq(otu_table(taxa.rel.abu.df4, taxa_are_rows = T), 
                       sample_data(samp.df), 
                       tax_table(as.matrix(taxa.df4.taxonomy.df)))          # Full filtered: 184 x 12
desert.phy2 <- tax_glom(desert.phy, taxrank=rank_names(desert.phy)[7])      # Agglomerate 165 x 12
desert.phy3 <- filter_taxa(desert.phy2, function(x) mean(x) > 0.01, TRUE)   # Only include mean abundances > 0.01%: 81 x 12

####
# Between all these hashes, I recode species that are the same at genus (and likely finer)
# to species-level classifications (e.g., Arthrobacter/Arthrobacter sp. and Arthrobacter/Unclassified Arthrobacter are the same).
# I then agglomerate to species to sidestep the multiplicity in species in the contig data.
# The possibly sketchy part is where I randomly select one megahit designation to subset the 
# contig sequence to just choose one to represent each taxon. The end result is 76 unique species.
# I needed to extract taxonomy and sequences for phylogenetic tree construction (Zohaib). I used the filtered and 
# subset df above to subset the SqueezeMeta object. Also this step is required to match the taxa names generated earlier
tax.df <- as.data.frame(desert.phy3@tax_table)
tax.df1 <- mutate_if(tax.df, is.factor, as.character)

## Abundances, taxonomy, and sequences
# Convert contigs to a df
full.sub5.tax <- as.data.frame(full.sub5$contigs$tax)
# full.sub5.tax$genus[grep(" bacterium \\(no genus in NCBI\\)", full.sub5.tax$genus)]
full.sub5.tax$genus <- paste("Unclassified ", gsub(" bacterium \\(no genus in NCBI\\)","", full.sub5.tax$genus), sep = "")
full.sub5.tax$genus <- paste("Unclassified ", gsub(" archaeon \\(no genus in NCBI\\)","", full.sub5.tax$genus), sep = "")
full.sub5.tax$genus <- gsub("Unclassified Unclassified ","Unclassified ", full.sub5.tax$genus)
full.sub5.tax$genus <- gsub("Unclassified Unclassified ","Unclassified ", full.sub5.tax$genus)
full.sub5.tax$family <- paste("Unclassified ", gsub(" bacterium \\(no family in NCBI\\)","", full.sub5.tax$family), sep = "")
full.sub5.tax$order <- paste("Unclassified ", gsub(" bacterium \\(no order in NCBI\\)","", full.sub5.tax$order), sep = "")
full.sub5.tax$class <- paste("Unclassified ", gsub(" bacterium \\(no class in NCBI\\)","", full.sub5.tax$class), sep = "")
full.sub5.tax$genus <- gsub("Unclassified Unclassified ","Unclassified ", full.sub5.tax$genus)
full.sub5.tax$genus <- gsub("Unclassified Unclassified ","Unclassified ", full.sub5.tax$genus)
full.sub5.tax$family <- gsub("Unclassified Unclassified ","Unclassified ", full.sub5.tax$family)
full.sub5.tax$family <- gsub("Unclassified Unclassified ","Unclassified ", full.sub5.tax$family)
full.sub5.tax$order <- gsub("Unclassified Unclassified ","Unclassified ", full.sub5.tax$order)
full.sub5.tax$order <- gsub("Unclassified Unclassified ","Unclassified ", full.sub5.tax$order)
full.sub5.tax$class <- gsub("Unclassified Unclassified ","Unclassified ", full.sub5.tax$class)
full.sub5.tax$class <- gsub("Unclassified Unclassified ","Unclassified ", full.sub5.tax$class)
full.sub5.tax$species[grep(" bacterium", full.sub5.tax$species)] <- paste("Unclassified ", gsub(" bacterium", "", full.sub5.tax$species[grep(" bacterium", full.sub5.tax$species)]), sep = "")
full.sub5.tax$species[grep(" archaeon", full.sub5.tax$species)] <- paste("Unclassified ", gsub(" archaeon", "", full.sub5.tax$species[grep(" archaeon", full.sub5.tax$species)]), sep = "")
full.sub5.tax[grep("Oxalobacteraceae",full.sub5.tax$species),] <- cbind.data.frame("Bacteria","Proteobacteria","Betaproteobacteria","Burkholderiales",	"Oxalobacteraceae",	"Massilia",	"Massilia violaceinigra B2")
full.sub5.tax[grep("Propionibacteriaceae",full.sub5.tax$species),] <- cbind.data.frame("Bacteria",	"Actinobacteria",	"Actinobacteria",	"Propionibacteriales",	"Propionibacteriaceae",	"Cutibacterium",	"Cutibacterium acnes")
full.sub5.tax[grep("Lawsonella",full.sub5.tax$species),] <- cbind.data.frame("Bacteria",	"Actinobacteria",	"Actinobacteria",	"Corynebacteriales",	"Lawsonellaceae",	"Lawsonella",	"Lawsonella clevelandensis")
full.sub5.tax[grep("Propionibacteriaceae",full.sub5.tax$species),] <- cbind.data.frame("Bacteria",	"Actinobacteria",	"Actinobacteria",	"Propionibacteriales",	"Propionibacteriaceae",	"Cutibacterium",	"Cutibacterium acnes")
full.sub5.tax[grep("Burkholderiaceae",full.sub5.tax$species),] <- cbind.data.frame("Bacteria",	"Proteobacteria",	"Betaproteobacteria",	"Burkholderiales",	"Burkholderiaceae",	"Paraburkholderia",	"Paraburkholderia caffeinilytica")
full.sub5.tax[grep("Ralstonia",full.sub5.tax$species),] <- cbind.data.frame("Bacteria",	"Proteobacteria",	"Betaproteobacteria",	"Burkholderiales",	"Burkholderiaceae",	"Ralstonia",	"Ralstonia pickettii")
full.sub5.tax[grep("Massilia",full.sub5.tax$species),] <- cbind.data.frame("Bacteria",	"Proteobacteria",	"Betaproteobacteria",	"Burkholderiales",	"Oxalobacteraceae",	"Massilia",	"Massilia oculi")
full.sub5.tax[grep("Nocardioidaceae",full.sub5.tax$species),] <- cbind.data.frame("Bacteria",	"Actinobacteria",	"Actinobacteria",	"Propionibacteriales",	"Nocardioidaceae",	"Marmoricola",	"Marmoricola scoriae")
full.sub5.tax[grep("Pseudomonas",full.sub5.tax$species),] <- cbind.data.frame("Bacteria",	"Proteobacteria",	"Gammaproteobacteria",	"Pseudomonadales",	"Pseudomonadaceae",	"Pseudomonas",	"Pseudomonas yamanorum")
full.sub5.tax[grep("Unclassified Staphylococcus",full.sub5.tax$species),] <- cbind.data.frame("Bacteria",	"Firmicutes",	"Bacilli",	"Bacillales",	"Staphylococcaceae",	"Staphylococcus",	"Staphylococcus capitis")
full.sub5.tax[grep("Burkholderia",full.sub5.tax$species),] <- cbind.data.frame("Bacteria",	"Proteobacteria",	"Betaproteobacteria",	"Burkholderiales",	"Burkholderiaceae",	"Burkholderia",	"Burkholderia gladioli")
full.sub5.tax[grep("Propionibacteriales",full.sub5.tax$species),] <- cbind.data.frame("Bacteria",	"Actinobacteria",	"Actinomycetia",	"Propionibacteriales",	"Nocardioidaceae",	"Nocardioides",	"unclassified Nocardioides")
full.sub5.tax[grep("Proteobacteria",full.sub5.tax$species),] <- cbind.data.frame("Bacteria",	"Proteobacteria",	"Deltaproteobacteria",	"Myxococcales",	"Myxococcaceae",	"Corallococcus", "Corallococcus coralloides")
full.sub5.tax[grep("Blastocatellia",full.sub5.tax$species),] <- cbind.data.frame("Bacteria",	"Acidobacteria",	"Acidobacteriia",	"Acidobacteriales",	"Acidobacteriacea",	"Acidobacterium",	"Acidobacterium capsulatum")
full.sub5.tax[grep("Marmoricola",full.sub5.tax$species),] <- cbind.data.frame("Bacteria",	"Actinobacteria",	"Actinobacteria",	"Propionibacteriales",	"Nocardioidaceae",	"Marmoricola",	"Marmoricola scoriae")
full.sub5.tax[grep("Gemmatimonadales",full.sub5.tax$species),] <- cbind.data.frame("Bacteria", "Actinobacteria", "Actinomycetia", "Corynebacteriales", "Nocardiaceae", "Nocardia","Nocardia farcinica")
full.sub5.tax[grep("Actinobacteria",full.sub5.tax$species),] <- cbind.data.frame("Bacteria", "Actinobacteria", "Actinomycetia", "Micrococcales", "Microbacteriaceae", "Protaetiibacter", "Unclassified Protaetiibacter")
full.sub5.tax$species[full.sub5.tax$species == "Arthrobacter sp."] <- "Unclassified Arthrobacter"
full.sub5.tax$species[full.sub5.tax$species == "Chryseobacterium sp."] <- "Unclassified Chryseobacterium"
full.sub5.tax$species[full.sub5.tax$species == "Aeromicrobium sp."] <- "Unclassified Aeromicrobium"
full.sub5.tax$species[full.sub5.tax$species == "Microbacterium sp."] <- "Unclassified Microbacterium"
full.sub5.tax$species[full.sub5.tax$species == "Micromonospora sp."] <- "Unclassified Micromonospora"
full.sub5.tax$species[full.sub5.tax$species == "Nocardioides sp."] <- "Unclassified Nocardioides"
full.sub5.tax$species[full.sub5.tax$species == "Nostoc sp."] <- "Unclassified Nostoc"
full.sub5.tax$species[full.sub5.tax$species == "Noviherbaspirillum sp."] <- "Unclassified Noviherbaspirillum"
full.sub5.tax$species[full.sub5.tax$species == "Propionibacterium sp."] <- "Unclassified Propionibacterium"
full.sub5.tax$species[full.sub5.tax$species == "Rhodococcus sp."] <- "Unclassified Rhodococcus"
full.sub5.tax$species[full.sub5.tax$species == "Streptococcus sp." ] <- "Unclassified Streptococcus"
levels(as.factor(full.sub5.tax$species[grep(" sp.", full.sub5.tax$species)])) # Check that all of the "XXXX sp." have been recoded.

# Subset the data to the 81 species (there will be more in the resultant df because of multiplicity n = 6915)
full.sub5.tax.df <- full.sub5.tax[full.sub5.tax$species %in% tax.df1$species,]    # Subset taxonomy
full.sub5.seqs <- full.sub5$contigs$seqs[names(full.sub5$contigs$seqs) %in% 
                                           rownames(full.sub5.tax.df)]            # Subset sequences
nlevels(as.factor(full.sub5.tax.df$phylum))                                       # 16 phyla
nlevels(as.factor(full.sub5.tax.df$genus))                                        # 74 genera
nlevels(as.factor(full.sub5.tax.df$species))                                      # 76 species
full.sub5.abu <- as.data.frame(full.sub5$contigs$abund)                           # abundances to df
full.sub5.abu <- full.sub5.abu[rownames(full.sub5.abu) %in% 
                                 rownames(full.sub5.tax.df),]  # Abundances       # Subset abundances
full.sub5.tax.df2 <- full.sub5.tax.df                                             # Create a copy to work with

# Here we select one megahit id to represent each taxon
full.sub5.tax.df2$row.names <- rownames(full.sub5.tax.df2)
set.seed(42)
contig.tax1.df3 <- full.sub5.tax.df2 %>%
  group_by(species) %>%
  slice_sample(n = 1, replace = F)
contig.tax1.df3$row.names

# Now agglomerate at the species level
identical(rownames(full.sub5.abu),rownames(full.sub5.tax.df)) # Check the rownames match b/w abundances and taxonomy
full.sub5.abu.df2 <- full.sub5.abu                            # Make a new df to work with
full.sub5.abu.df2$species <- full.sub5.tax.df2$species        # Assign species names as a column
contig.abu1.df3 <- full.sub5.abu.df2 %>%
  group_by(species) %>%
  summarise(across(D734Y:G624Y, ~sum(., na.rm = T)))          # Sum non-unique species (absolute abundances)

# Now use the above df to subset the sequence file
contig.seq1.df3 <- full.sub5.seqs[names(full.sub5.seqs) %in% contig.tax1.df3$row.names]
names(contig.seq1.df3)
contig.tax1.df3 <- as.data.frame(contig.tax1.df3)             # Convert tibble to df
rownames(contig.tax1.df3) <- contig.tax1.df3$row.names        # Assign row names
contig.abu1.df3 <- as.data.frame(contig.abu1.df3)             # Convert tibble to df

# Double check for match before assigning row names
identical(contig.abu1.df3$species, contig.tax1.df3$species)
rownames(contig.abu1.df3) <- rownames(contig.tax1.df3)
nlevels(as.factor(contig.tax1.df3$species)) # Should be 76

# This step is to identify which sequences are too long and creating problems for the multiple sequence alignment that Zohaib is working on
# seq.lens <- cbind.data.frame(names(contig.seq1.df3),nchar(contig.seq1.df3))
# seq.lens[seq.lens$`nchar(contig.seq1.df3)`>9000,]

# Output the files to use for phylogenetic analyses if necessary
write.table(contig.seq1.df3, "~/Dropbox/Cyclotron Research/Carbon_Oxidation_and_Fixation_Study/Metagenome/SteveM/Curated_Rcode/seq.species.glom1.fasta")
write.csv(contig.abu1.df3, "~/Dropbox/Cyclotron Research/Carbon_Oxidation_and_Fixation_Study/Metagenome/SteveM/Curated_Rcode/abu.species.glom1.csv")
write.csv(contig.tax1.df3, "~/Dropbox/Cyclotron Research/Carbon_Oxidation_and_Fixation_Study/Metagenome/SteveM/Curated_Rcode/tax.species.glom1.csv")


##**********************
#### Here I will make the taxonomic tree based on species-level agglomeration

## Prepare the data
abu.df <- as.data.frame(desert.phy3@otu_table)                        # Abundances
sam.df <- as.data.frame(desert.phy3@sam_data)                         # Samples
tax.df <- as.data.frame(desert.phy3@tax_table)                        # Taxonomy
tax.df <- mutate_if(tax.df, is.character, as.factor)                  # as.phylo requires factors rather than strings

taxa.rel.abu.df0a <- abu.df                                           # Create a copy to work with
taxa.rel.abu.df0a$species <- gsub(".1","", 
                                  gsub(".3","", 
                                       rownames(taxa.rel.abu.df0a)))  # Nested gsub to remove 'make.unique' numbers and assign row names
taxa.rel.abu.df0a <- taxa.rel.abu.df0a %>%
  group_by(species) %>%
  summarise(across(D734Y:G624Y, ~sum(., na.rm = T)))                  # Agglomerate the relative abundances by species
taxa.rel.abu.df0a <- as.data.frame(taxa.rel.abu.df0a)                 # Convert tibble to df
rownames(taxa.rel.abu.df0a) <- taxa.rel.abu.df0a$species              # Row names
taxa.rel.abu.df0a$species <- NULL                                     # Remove the species column
contig.abu1.df3a <- taxa.rel.abu.df0a[rownames(taxa.rel.abu.df0a) 
                                      %in% contig.tax1.df3$species,]  # Subset abundance to taxonomy
range(colSums(contig.abu1.df3a))                                      # This will show the summed relative abundances of the remaining taxa
identical(rownames(contig.abu1.df3a), contig.tax1.df3$species)        # Double check matches

# Calculate mean relative abundances and fold-changes for the background and methane hotspots
abu.df1.means <- cbind.data.frame("background" = rowMeans(contig.abu1.df3a[c(1,4,8,12)]),
                                 "methane" = rowMeans(contig.abu1.df3a[c(2,6,7,11,3,5,9,10)]))
abu.df1.means$ch4.fold <- foldchange(abu.df1.means$methane, abu.df1.means$background) # Fold-changes between the background and methane
abu.df1.means$abs.ch4.fold <- abs(abu.df1.means$ch4.fold) # Absolute values
range(abs(abu.df1.means$abs.ch4.fold[!is.infinite(abu.df1.means$abs.ch4.fold)]))  # 1.024953 to 378.203447
hist(abs(abu.df1.means$abs.ch4.fold[!is.infinite(abu.df1.means$abs.ch4.fold)]), breaks = seq(0,400,5), ylim = c(0,50))
abu.df1.means[order(-abu.df1.means$ch4.fold),]  # Sort the fold changes for visual inspection of the table
abu.df1.means[abu.df1.means$ch4.fold > 100,]

contig.tax1.df3a <- contig.tax1.df3                                                       # Make a new df to work with
rownames(contig.tax1.df3a) <- contig.tax1.df3a$species                                    # Assign species as row names for the taxonomy
rownames(abu.df1.means) <- contig.tax1.df3a$species                                       # Do the same for relative abundances
merge1a <- merge(contig.tax1.df3a, abu.df1.means, by = 0)                                 # Merge the taxonomy and relative abundances
rownames(merge1a) <- merge1a$Row.names                                                    # Assign row names
merge1a$Row.names <- NULL                                                                 # Remove the row names column
contig.tax1.df3a <- mutate_if(contig.tax1.df3a, is.character, as.factor)                  # If character string, convert to factor

#### This is code to check the phylograms produced by Zohaib. More of a side bar.
# desert.tree <- read.tree("~/Dropbox/Cyclotron Research/Carbon_Oxidation_and_Fixation_Study/Metagenome/SteveM/Phylogenetic_tree/seq.species.glom-no-matrix-upgma.nwk")
# tax.glom.csv <- read.csv("~/Dropbox/Cyclotron Research/Carbon_Oxidation_and_Fixation_Study/Metagenome/SteveM/Curated_Rcode/Zohaib/tax.species.glom.csv")
# desert.tree$tip.label <- gsub(".*-","",desert.tree$tip.label)
# 
# tax.glom.csv$species <- gsub(" ","_", tax.glom.csv$species)
# rownames(tax.glom.csv) <- tax.glom.csv$species
# tax.glom.csv <- tax.glom.csv[desert.tree$tip.label,]
# identical(tax.glom.csv$species, desert.tree$tip.label)
# 
# p <- ggtree(desert.tree, layout = "circular", size = 0.2)
# p + geom_label(aes(x=branch, label=NA)) + 
#   geom_tiplab(aes(angle = angle, cex = 1, label = gsub(" ","_",label)), offset = .15, parse = T, show.legend=F)
# 
# cls <- list(Acidobacteria = c(tax.glom.csv$species[tax.glom.csv$phylum == "Acidobacteria"]),
#             Actinobacteria = c(tax.glom.csv$species[tax.glom.csv$phylum == "Actinobacteria"]),
#             Ascomycota = c(tax.glom.csv$species[tax.glom.csv$phylum == "Ascomycota"]),
#             Bacteroidetes = c(tax.glom.csv$species[tax.glom.csv$phylum == "Bacteroidetes"]),
#             Candidatus_Rokubacteria = c(tax.glom.csv$species[tax.glom.csv$phylum == "Candidatus Rokubacteria"]),
#             Candidatus_Saccharibacteria = c(tax.glom.csv$species[tax.glom.csv$phylum == "Candidatus Saccharibacteria"]),
#             Chloroflexi = c(tax.glom.csv$species[tax.glom.csv$phylum == "Chloroflexi"]),
#             Cyanobacteria = c(tax.glom.csv$species[tax.glom.csv$phylum == "Cyanobacteria"]),
#             Firmicutes = c(tax.glom.csv$species[tax.glom.csv$phylum == "Firmicutes"]),
#             Gemmatimonadetes = c(tax.glom.csv$species[tax.glom.csv$phylum == "Gemmatimonadetes"]),
#             Nitrospirae = c(tax.glom.csv$species[tax.glom.csv$phylum == "Nitrospirae"]),
#             Planctomycetes = c(tax.glom.csv$species[tax.glom.csv$phylum == "Planctomycetes"]),
#             Proteobacteria = c(tax.glom.csv$species[tax.glom.csv$phylum == "Proteobacteria"]),
#             Thaumarchaeota = c(tax.glom.csv$species[tax.glom.csv$phylum == "Thaumarchaeota"]),
#             Archaea = c(tax.glom.csv$species[tax.glom.csv$phylum == "Unclassified Archaea"]),
#             Bacteria = c(tax.glom.csv$species[tax.glom.csv$phylum == "Unclassified Bacteria"]),
#             Verrucomicrobia = c(tax.glom.csv$species[tax.glom.csv$phylum == "Verrucomicrobia"]))
# 
# 
# tree <- groupOTU(desert.tree, cls)
# col.pals <- c(wes_palette("Darjeeling1", nlevels(as.factor(tax.glom.csv$phylum)), type = "continuous"))
# ggtree(tree, aes(color=group), layout = "circular", size = 0.2) + geom_tiplab() +
#   scale_color_manual(values=col.pals) + theme(legend.position="right")


# Now can finally make the tree using ape
frm1 <- ~superkingdom/phylum/class/order/family/genus/species   # Create the form (KPCOFGS) to use for the tree resolution (i.e., to species)
tr1 <- ape::as.phylo(frm1, data = contig.tax1.df3a)             # Make the tree
tr1$tip.label <- gsub(" ","_",tr1$tip.label)
rownames(merge1a) <- gsub(" ","_",rownames(merge1a))
merge2a <- merge1a[tr1$tip.label,]                              # Order the merged file by the tip labels (i.e., species)
identical(rownames(merge2a), tr1$tip.label)                     # Always double check they match

## Make the tree and label the nodes sequentially to determine coloring, etc. below
p <- ggtree(tr1, layout = "circular", size = 0.2)
p + geom_label(aes(x=branch, label=node)) + 
  geom_tiplab(aes(angle = angle, cex = 1, label = gsub(" ","_",label)), offset = .15, parse = T, show.legend=F)

# I inspected the tree from above and manually assigned the colors based on phylum
acidobacteria <- c(1:4,79)
actinobacteria <- c(19:51,87:102)
ascomycota <- 76
bacteroidetes <- c(59:61,106)
candidatus_rokubacteria <- 62
candidatus_saccharibacteria <- 63
chloroflexi <- c(64:66,107)
cyanobacteria <- c(67:69,108,109)
firmicutes <- c(52:57,103:105)
gemmatimonadetes <- 70
# nitrospirae <- 0
planctomycetes <- 71
proteobacteria <- c(5:18,80:86)
thaumarchaeota <- c(74,75,111)
unclassified_archaea <- 73
unclassified_bacteria <- 58
verrucomicrobia <- 72
black.col <- c(77,78,110)

# This is just a check that I've got everything (the row number should match with the assigned number)
sort(c(acidobacteria,
       actinobacteria,
       ascomycota,
       bacteroidetes,
       candidatus_rokubacteria,
       candidatus_saccharibacteria,
       chloroflexi,
       cyanobacteria,
       firmicutes,
       gemmatimonadetes,
       # nitrospirae,
       planctomycetes,
       proteobacteria,
       thaumarchaeota,
       unclassified_archaea,
       unclassified_bacteria,
       verrucomicrobia,
       black.col))

# Make a color palette
col.pals <- c(wes_palette("Darjeeling1", nlevels(contig.tax1.df3a$phylum), type = "continuous"))

# Prepare the phylum color for the tree
d <- data.frame(node = 1:111, color = "black")
d$color[acidobacteria] <- col.pals[1]
d$color[actinobacteria] <- col.pals[2]
d$color[ascomycota] <- col.pals[3]
d$color[bacteroidetes] <- col.pals[4]
d$color[candidatus_rokubacteria] <- col.pals[5]
d$color[candidatus_saccharibacteria] <- col.pals[6]
d$color[chloroflexi] <- col.pals[7]
d$color[cyanobacteria] <- col.pals[8]
d$color[firmicutes] <- col.pals[9]
d$color[gemmatimonadetes] <- col.pals[10]
d$color[planctomycetes] <- col.pals[11]
d$color[proteobacteria] <- col.pals[12]
d$color[thaumarchaeota] <- col.pals[13]
d$color[unclassified_archaea] <- col.pals[14]
d$color[unclassified_bacteria] <- col.pals[15]
d$color[verrucomicrobia] <- col.pals[16]

# Now to make a new tree that includes the color information
p <- ggtree(tr1,layout = "circular", size = 0.2) %<+% d + aes(color=I(color)) + xlim(-0.5, 9.1)

# Add relative abundances
p$data$background <- NA
p$data$background[c(1:76)] <- merge2a$background
p$data$methane <- NA
p$data$methane[c(1:76)] <- merge2a$methane
p$data$phylum <- NA
p$data$phylum[c(1:76)] <- as.character(gsub(" ", "_", merge2a$phylum))
p$data$species <- NA
p$data$species[c(1:76)] <- as.character(gsub(" \\(no species in NCBI)", "", as.character(merge2a$species)))
p$data$species[c(1:76)] <- as.character(gsub(" ", "_", p$data$species[c(1:76)]))

## I import this tree to Photoshop to beautify it
jpeg("tree1.jpg", width = 8, height = 8, units = "in", res = 600)
p2 <- p + geom_tiplab(aes(angle = angle, cex = 1, label = paste0('italic(', species,')')), size = 2, offset = 1.75, parse = T, show.legend = F) +
  theme_tree(plot.margin = margin(0, 220, 0, 60)) +
  geom_tippoint(aes(x = p$data$x+0.3, size = background, color = color), alpha = 0.65, stroke = 0, show.legend = T) +
  geom_tippoint(aes(x = p$data$x+0.95, size = methane, color = color), alpha = 0.65, stroke = 0, show.legend = F) +
  scale_size_continuous(range = c(0.1, 6)) +
  geom_text(aes(x = 0, y = 0, label = ""), show.legend=F) +
  theme(legend.position = c(1.5,0.5)) +
  labs(size = "Relative abundance (%)", colour = "Phylum")
p2
dev.off()

##**********************
## Here I look at pathways of interest in the Ralstonia bins
# "maxbin.161.fasta.contigs" "maxbin.104.fasta.contigs"
Ralstonia.bins.sub <- subsetBins(full.sub5, Ralstonia.bins)
Ralstonia.bins.sub$taxa$species
Ralstonia.bins.sub$functions$KEGG$tpm       # 1088 x 12
dim(Ralstonia.bins.sub$taxa$species$abund)  # 15 x 12

## List of some pathways of interest
# "K03319","K11616","K24180","K03300" = citP
# "K01673","K01725","K11921","K15576","K15577","K15579" = cynT, cynS, cynR, cynA, cynB, cynD
# "K02575","K10850","K15576","K15577","K15578","K15579","K00370","K00371","K00374" = nitrate transport
# "K10944","K10945","K10946","K16157","K16158" = methane/ammonia monooxygenase

# Exploring to see what's there
Ralstonia.bins.sub$functions$KEGG$tpm[rownames(Ralstonia.bins.sub$functions$KEGG$tpm) %in% 
                                        c("K03319","K11616","K24180","K03300",                             # citP
                                          "K01725","K02575","K10850","K15576","K15577","K15578","K15579",  # Nitrogen metabolism (cynT, cynS, cynR, cynA, cynB, cynD)
                                          "K00370","K00371","K00374",                                      # Nitrate transport (narG, narZ, nxrA, narH, narY, nxrB, narI, narV)
                                          "K10944"),]                                                      # Methane metabolism
full.sub5$functions$KEGG$tpm[rownames(full.sub5$functions$KEGG$tpm) == "K10944",]
full.sub5$functions$KEGG$tpm[rownames(full.sub5$functions$KEGG$tpm) == "K10946",]
mmo.sub <- subsetFun(full.sub5, c("K10944|K10945|K10946|K16157|K16158"))
mmo.sub$functions$KEGG$tpm 
mmo.sub$contigs$tax
full.sub5$contigs$tax[full.sub5$contigs$tax]
full.sub5$contigs$tax[grep("Nitrosomonadales", full.sub5$contigs$tax[,4]),]
full.sub5$contigs$tax[grep("Nitrosomonadaceae", full.sub5$contigs$tax[,5]),]
full.sub5$contigs$tax[grep("Nitrosomonas", full.sub5$contigs$tax[,6]),]
desert$contigs$tax[grep("Nitrosomonas", desert$contigs$tax[,6]),]


## 81 KEGG pathways in total that we are interested in
# From Black and Black
# Methylotrophic bacteria can fix nitrogen when provided with methane, methanol, or hydrogen from various substrates.
ko.list <- c("K00024|K00050|K00058|K00093|K00121|K00123|K00124|K00125|K00127|K00148|K00261|K00265|K00266|K00300|K00362|K00363|K00368|K00370|K00371|K00372|K00374|K00376|K00459|K00600|K00830|K00831|K00865|K00925|K01007|K01070|K01079|K01499|K01595|K01624|K01647|K01673|K01689|K01690|K01673|K01725|K01834|K01895|K02567|K02568|K02575|K03300|K03319|K03385|K03841|K04561|K05299|K05979|K10713|K10714|K10850|K10944|K10945|K10946|K11529|K11616|K11921|K13812|K14028|K14029|K15022|K15371|K15576|K15577|K15578|K15579|K15633|K15634|K15635|K15876|K15918|K16157|K16158|K16370|K22516|K22982|K23995|K24180")

## There are 18 of the 81 pathways present in Ralstonia
Ralstonia.bins.sub.kegg <- subsetFun(Ralstonia.bins.sub, ko.list)
Ralstonia.kegg.t <- as.data.frame(t(Ralstonia.bins.sub.kegg$functions$KEGG$tpm[,c(1,4,8,12,2,6,7,11,3,5,9,10)]))  # Convert to df

# Make the grouping factor (first 4 rows are background, last 8 are methane)
groups1 <- as.factor(c(rep("background", 4), rep("methane", 8)))

###********
## Comparing background and methane genera using Poisson regressions
A <- round(Ralstonia.kegg.t+1,0)
B <- groups1 # Grouping factor
# Loop the regression
ch4.pu = matrix(nrow=ncol(A), ncol=3) 
count = 1 # Start at the beginning eh?
for( i in 1:ncol(A) ){
  # Poisson regression
  zip.pu <- glm(A[,i] ~ B, family = poisson)
  # Add OTU names, slope estimate, and p-values
  ch4.pu[count,] = c(colnames(A)[i],coef(summary(zip.pu))[2,1],coef(summary(zip.pu))[2,4])
  # Increment the count
  count = count + 1
}

ch4.pu <- as.data.frame(ch4.pu)                               # Convert the output to a data frame
colnames(ch4.pu) <- c("kegg","Slope","Pvalue")                # Rename the columns
ch4.pu$Slope <- as.numeric(as.character(ch4.pu$Slope))        # Convert strings to numeric
ch4.pu$Pvalue <- as.numeric(as.character(ch4.pu$Pvalue))      # Convert strings to numeric
ch4.pu$bonf <- p.adjust(ch4.pu$Pvalue, method = "bonferroni") # Bonferroni correction on the p-values
ch4.pu.sig <- subset(ch4.pu, bonf < 0.1)                      # Remove any rows that aren't significant (bonf)
ch4.pu.sig1 <- subset(ch4.pu, Pvalue < 0.1)                   # Remove any rows that aren't significant (uncorrected)

##*********************************
## Run a PCA to compare composition among the various samples

abu.df.glom0 <- contig.abu1.df3                 # Create an abundance copy to work with
rownames(abu.df.glom0) <- abu.df.glom0$species  # Assign species as row names
abu.df.glom0$species <- NULL                    # Remove the species column
tax.df.glom0 <- contig.tax1.df3                 #
rownames(tax.df.glom0) <- tax.df.glom0$species  # Create a taxonomy copy to work with

# Make our compositional dataset (taxa must be columns)
system.time(glom.df0.gbm <- cmultRepl(t(abu.df.glom0), 
                                      method = "GBM", 
                                      output = "p-counts")) # Replace zeroes using GBM
glom.df0.gbm.clr <- codaSeq.clr(glom.df0.gbm, 
                                samples.by.row = T)         # Log-centred ratios
identical(rownames(samp.df), rownames(glom.df0.gbm.clr))    # Check the sample and CLR df row names match
glom.df0.gbm.clr1 <- cbind.data.frame(samp.df, 
                                      glom.df0.gbm.clr)     # Merge samples and CLR

# Run the PCA
glom.df0.gbm.clr1.pca <- rda(glom.df0.gbm.clr1[,c(6:ncol(glom.df0.gbm.clr1))])

# Sum the total variance
glom.df0.gbm.clr1.pca.mvar <- sum(glom.df0.gbm.clr1.pca$CA$eig)

# Calculate the PC1 and PC2 variance
glom.df0.gbm.clr1.pca.PC1 <- paste("PC1: ", sprintf("%.3f",(sum(glom.df0.gbm.clr1.pca$CA$eig[1])/glom.df0.gbm.clr1.pca.mvar), 3))
glom.df0.gbm.clr1.pca.PC2 <- paste("PC2: ", sprintf("%.3f",(sum(glom.df0.gbm.clr1.pca$CA$eig[2])/glom.df0.gbm.clr1.pca.mvar), 3))

# Extract the plotting coordinates for the biplot
glom.df0.gbm.clr1.pca.scrs <- vegan::scores(glom.df0.gbm.clr1.pca, display = c("sites","species"), scaling = 3) # Symmetric scaling
glom.df0.gbm.clr1.pca.xlim <- with(glom.df0.gbm.clr1.pca.scrs, range(sites[,1]))
glom.df0.gbm.clr1.pca.ylim <- with(glom.df0.gbm.clr1.pca.scrs, range(sites[,2]))
rownames(glom.df0.gbm.clr1.pca.scrs$species)
identical(rownames(glom.df0.gbm.clr1.pca.scrs$species), tax.df.glom0$species)     # TRUE
tax.df.glom0 <- tax.df.glom0[rownames(glom.df0.gbm.clr1.pca.scrs$species),]
identical(rownames(glom.df0.gbm.clr1.pca.scrs$species), tax.df.glom0$species)    # TRUE

# Assign treatment types
glom.df0.gbm.clr1$sample.type1 <- c("bulk","methane","methane","bulk","methane","methane","methane","bulk","methane","methane","methane","bulk")

# Create centroids to use to plot phyla
tmp.crp$Acidobacteria$center
tmp.crp0 <- unlist(tmp.crp)
tmp.pc1 <- tmp.crp0[grep("center.PC1",names(tmp.crp0))]
tmp.pc2 <- tmp.crp0[grep("center.PC2",names(tmp.crp0))]
spp.centroids <- cbind.data.frame(tmp.pc1,tmp.pc2)
rownames(spp.centroids) <- gsub(".center.PC1","",rownames(spp.centroids))

# Create the PCA graphic to export to Photoshop. Export at 4.7 x 4.5
jpeg("PCA_tree1.jpg", width = 4.7, height = 4.5, units = "in", res = 600)
par(mar = c(3.25,3.25,0.2,0.2))
plot.new()
plot.window(xlim = glom.df0.gbm.clr1.pca.xlim, ylim = glom.df0.gbm.clr1.pca.ylim, asp = 1)
abline(h = 0, lty = "dotted", col = "gray75"); abline(v = 0, lty = "dotted", col = "gray75")
colvec1 <- wes_palette("Rushmore1")[c(4,3)]
vpch <- c(16,17)
with(glom.df0.gbm.clr1, 
     points(glom.df0.gbm.clr1.pca.scrs$sites, 
            pch = vpch[factor(glom.df0.gbm.clr1$diapir.presence)], 
            col = colvec1[factor(glom.df0.gbm.clr1$sample.type1)], 
            cex = 1))
text(spp.centroids$tmp.pc1, spp.centroids$tmp.pc2,
            labels = rownames(spp.centroids),
            col = col.pals0,
            cex = 0.6)
axis(side = 1, cex = 0.8); axis(side = 2, cex = 1)
mtext(side = 1, glom.df0.gbm.clr1.pca.PC1, line = 2.2, cex = 1)
mtext(side = 2, glom.df0.gbm.clr1.pca.PC2, line = 2.2, cex = 1)
box()
# legend("bottomright",
#        c("Bulk", expression(paste("CH"[4]))),
#        pch = 17, pt.cex = 1,
#        text.col = c("#35274A", "#0B775E"), text.font = 2, title = "G",
#        col = c("#35274A", "#0B775E"), bty = "n", inset = c(0,-0.0), x.intersp = 1, cex = 1)
# legend("bottomright",
#        rep("",2), pch = 16, pt.cex = 1,
#        col = c("#35274A", "#0B775E"), bty = "n", title = "D",
#        inset = c(0.34,-0.0), cex = 1)
dev.off()
  

##***************************************
## Plot completeness/contamination for the supp

## Read in the data
des.col <- read.csv("~/Dropbox/Cyclotron Research/Carbon_Oxidation_and_Fixation_Study/Metagenome/SteveM/Curated_Rcode/des.col.csv", row.names = 1)
cont.col <- read.csv("~/Dropbox/Cyclotron Research/Carbon_Oxidation_and_Fixation_Study/Metagenome/SteveM/Curated_Rcode/cont.col.csv", row.names = 1)
Completeness <- read.csv("~/Dropbox/Cyclotron Research/Carbon_Oxidation_and_Fixation_Study/Metagenome/SteveM/Curated_Rcode/Completeness.csv", row.names = 1)
Contamination <- read.csv("~/Dropbox/Cyclotron Research/Carbon_Oxidation_and_Fixation_Study/Metagenome/SteveM/Curated_Rcode/Contamination.csv", row.names = 1)
desert.df <- cbind.data.frame(Completeness = Completeness$x, Contamination = Contamination$x)
des.col <- des.col$x
des.col2 <- des.col
cont.col <- cont.col$x

# Create color indices
des.col2[des.col <= 50] <- "black"
des.col2[des.col > 90] <- wes_palette("Rushmore1")[1]
des.col2[des.col <= 90 & des.col > 70] <- wes_palette("Rushmore1")[3]
des.col2[des.col <= 70 & des.col > 50] <- wes_palette("Rushmore1")[5]

### Create the various completeness/contamination subsets to use (50-70%, 70-90%, >90%)

## Completeness
# 50-70%
comp.d1 <- des.col[des.col <= 70 & des.col > 50]
comp.d1a <- density(comp.d1)
comp.d1a <- cbind.data.frame(x = comp.d1a$x, y = comp.d1a$y)
comp.d1a <- rbind.data.frame(c(x = 0, y = 0), comp.d1a, c(x = 100, y = 0))
# 70-90%
comp.d2 <- des.col[des.col <= 90 & des.col > 70]
comp.d2a <- density(comp.d2)
comp.d2a <- cbind.data.frame(x = comp.d2a$x, y = comp.d2a$y)
comp.d2a <- rbind.data.frame(c(x = 0, y = 0), comp.d2a, c(x = 100, y = 0))
# >90%
comp.d3 <- des.col[des.col > 90]
comp.d3a <- density(comp.d3)
comp.d3a <- cbind.data.frame(x = comp.d3a$x, y = comp.d3a$y)
comp.d3a <- rbind.data.frame(c(x = 0, y = 0), comp.d3a, c(x = 100, y = 0))

## Contamination
# 50-70%
cont.d1 <- cont.col[des.col <= 70 & des.col > 50]
cont.d1a <- density(cont.d1)
cont.d1a <- cbind.data.frame(x = cont.d1a$x, y = cont.d1a$y)

# 70-90%
cont.d2 <- cont.col[des.col <= 90 & des.col > 70]
cont.d2a <- density(cont.d2)
cont.d2a <- cbind.data.frame(x = cont.d2a$x, y = cont.d2a$y)

# >90%
cont.d3 <- cont.col[des.col > 90]
cont.d3a <- density(cont.d3)
cont.d3a <- cbind.data.frame(x = cont.d3a$x, y = cont.d3a$y)

## Make a figure to compare completeness and contamination of the MAGs
jpeg("FigureS04.jpg", width = 6.2, height = 6, units = "in", res = 300)
m <- rbind(c(2, 2, 4), 
           c(3, 3, 1),
           c(3, 3, 1))
print(m)

layout(m)
par(mar = c(0.1, 0.1, 0.1, 0.1), oma = c(4,4,1,1))
plot(cont.d1a$y, cont.d1a$x, type = "n", ylim = c(0.00, 44.08), xlim = c(0,0.06), axes = F)
polygon(cont.d1a$y, cont.d1a$x, col = alpha(wes_palette("Rushmore1")[5],0.5), lty = 0)
polygon(cont.d2a$y, cont.d2a$x, col = alpha(wes_palette("Rushmore1")[3],0.5), lty = 0)
polygon(cont.d3a$y, cont.d3a$x, col = alpha(wes_palette("Rushmore1")[1],0.5), lty = 0)
axis(1)

plot(comp.d1a$x, comp.d1a$y, type = "n", xlim = c(0,100), ylim = c(0,0.15), axes = F)
polygon(comp.d1a$x, comp.d1a$y, col = alpha(wes_palette("Rushmore1")[5],0.5), lty = 0)
polygon(comp.d2a$x, comp.d2a$y, col = alpha(wes_palette("Rushmore1")[3],0.5), lty = 0)
polygon(comp.d3a$x, comp.d3a$y, col = alpha(wes_palette("Rushmore1")[1],0.5), lty = 0)
axis(2)

plot(desert.df$Completeness, desert.df$Contamination, pch = 16, cex = 1.5, col = alpha(des.col2,0.5), xlim = c(0,100), ylim = c(0.00, 44.08))
mtext("Completeness", side = 1, line = 2.2)
mtext("Contamination", side = 2, line = 2.2)
dev.off()

##***************************************
## Plot the individual histograms and scatterplots for the filter figure in the supp

## Zeroes

# KEGG TPM
sum(kegg.tpm.df0 == 0); print("of"); dim(kegg.df0)[1]*dim(kegg.df0)[2]  # 52257 of 168204
sum(kegg.tpm.df0 == 0)/(dim(kegg.df0)[1]*dim(kegg.df0)[2])              # 31%
sum(kegg.tpm.df1 == 0); print("of"); dim(kegg.df1)[1]*dim(kegg.df1)[2]  # 29477 of 134280
sum(kegg.tpm.df1 == 0)/(dim(kegg.df1)[1]*dim(kegg.df1)[2])              # 22%
sum(kegg.tpm.df2 == 0); print("of"); dim(kegg.df2)[1]*dim(kegg.df2)[2]  # 11427 of 80364
sum(kegg.tpm.df2 == 0)/(dim(kegg.df2)[1]*dim(kegg.df2)[2])              # 14%
sum(kegg.tpm.df3 == 0); print("of"); dim(kegg.df3)[1]*dim(kegg.df3)[2]  # 11427 of 80364
sum(kegg.tpm.df3 == 0)/(dim(kegg.df3)[1]*dim(kegg.df3)[2])              # 14%
sum(kegg.tpm.df4 == 0); print("of"); dim(kegg.df4)[1]*dim(kegg.df4)[2]  # 10465 of 79248
sum(kegg.tpm.df4 == 0)/(dim(kegg.df4)[1]*dim(kegg.df4)[2])              # 13%

# Taxa
sum(taxa.df0 == 0); print("of"); dim(taxa.df0)[1]*dim(taxa.df0)[2]  # 23732 of 53136
sum(taxa.df0 == 0)/(dim(taxa.df0)[1]*dim(taxa.df0)[2])              # 45%
sum(taxa.df1 == 0); print("of"); dim(taxa.df1)[1]*dim(taxa.df1)[2]  # 22100 of 49464
sum(taxa.df1 == 0)/(dim(taxa.df1)[1]*dim(taxa.df1)[2])              # 45%
sum(taxa.df2 == 0); print("of"); dim(taxa.df2)[1]*dim(taxa.df2)[2]  # 2361 of 7572
sum(taxa.df2 == 0)/(dim(taxa.df2)[1]*dim(taxa.df2)[2])              # 31%
sum(taxa.df3 == 0); print("of"); dim(taxa.df3)[1]*dim(taxa.df3)[2]  # 2361 of 7572
sum(taxa.df3 == 0)/(dim(taxa.df3)[1]*dim(taxa.df3)[2])              # 31%
sum(taxa.df4 == 0); print("of"); dim(taxa.df4)[1]*dim(taxa.df4)[2]  # 2132 of 7308
sum(taxa.df4 == 0)/(dim(taxa.df4)[1]*dim(taxa.df4)[2])              # 29%

## KEGG
# Hist 0
jpeg("hist0.jpg", width = 4.25, height = 4, units = "in", res = 300)
par(mar = c(2,2,0.5,0.5))
h0 <- hist(log1p(as.matrix(kegg.tpm.df0)), freq = TRUE, breaks = seq(0,20,0.5), 
           xlim = c(0,20), ylim = c(0,16000), main = "", ylab = "", xlab = "", col = wes_palette("FantasticFox1")[1])
h1 <- hist(log1p(as.matrix(kegg.tpm.df1)), freq = TRUE, breaks = seq(0,20,0.5), 
           xlim = c(0,20), ylim = c(0,16000), main = "", ylab = "", xlab = "", col = wes_palette("FantasticFox1")[2], add = T)
box()
dev.off()
# Hist 1
jpeg("hist1.jpg", width = 4.25, height = 4, units = "in", res = 300)
par(mar = c(2,2,0.5,0.5))
h1 <- hist(log1p(as.matrix(kegg.tpm.df1)), freq = TRUE, breaks = seq(0,20,0.5), 
     xlim = c(0,20), ylim = c(0,16000), main = "", ylab = "", xlab = "", col = wes_palette("FantasticFox1")[2])
h2 <- hist(log1p(as.matrix(kegg.tpm.df2)), freq = TRUE, breaks = seq(0,20,0.5), 
     xlim = c(0,20), ylim = c(0,16000), main = "", ylab = "", xlab = "", col = wes_palette("FantasticFox1")[3], add = T)
box()
dev.off()
# Hist 2
jpeg("hist2.jpg", width = 4.25, height = 4, units = "in", res = 300)
par(mar = c(2,2,0.5,0.5))
h2 <- hist(log1p(as.matrix(kegg.tpm.df2)), freq = TRUE, breaks = seq(0,20,0.5), 
     xlim = c(0,20), ylim = c(0,16000), main = "", ylab = "", xlab = "", col = wes_palette("FantasticFox1")[3])
h3 <- hist(log1p(as.matrix(kegg.tpm.df3)), freq = TRUE, breaks = seq(0,20,0.5), 
     xlim = c(0,20), ylim = c(0,16000), main = "", ylab = "", xlab = "", col = wes_palette("FantasticFox1")[5], add = T)
box()
dev.off()
# Hist 3
jpeg("hist3.jpg", width = 4.25, height = 4, units = "in", res = 300)
par(mar = c(2,2,0.5,0.5))
h3 <- hist(log1p(as.matrix(kegg.tpm.df3)), freq = TRUE, breaks = seq(0,20,0.5), 
     xlim = c(0,20), ylim = c(0,16000), main = "", ylab = "", xlab = "", col = wes_palette("FantasticFox1")[5])
h4 <- hist(log1p(as.matrix(kegg.tpm.df4)), freq = TRUE, breaks = seq(0,20,0.5), 
     xlim = c(0,20), ylim = c(0,16000), main = "", ylab = "", xlab = "", col = wes_palette("Rushmore1")[3], add = T)
box()
dev.off()

zero.kegg.bin <- rbind.data.frame(h0$counts[1],
                                  h1$counts[1],
                                  h2$counts[1],
                                  h3$counts[1],
                                  h4$counts[1])
names(zero.kegg.bin) <- "zero.bin"

## Taxa
# Hist 0
jpeg("hist0.tax.jpg", width = 4.25, height = 4, units = "in", res = 300)
par(mar = c(2,2,0.5,0.5))
h0.tax <- hist(log1p(as.matrix(taxa.df0)), freq = TRUE, breaks = seq(0,20,0.5), 
               xlim = c(0,20), ylim = c(0,5000), main = "", ylab = "", xlab = "", col = wes_palette("FantasticFox1")[1])
h1.tax <- hist(log1p(as.matrix(taxa.df1)), freq = TRUE, breaks = seq(0,20,0.5), 
               xlim = c(0,20), ylim = c(0,5000), main = "", ylab = "", xlab = "", col = wes_palette("FantasticFox1")[2], add = T)
box()
dev.off()
# Hist 1
jpeg("hist1.tax.jpg", width = 4.25, height = 4, units = "in", res = 300)
par(mar = c(2,2,0.5,0.5))
h1.tax <- hist(log1p(as.matrix(taxa.df1)), freq = TRUE, breaks = seq(0,20,0.5), 
     xlim = c(0,20), ylim = c(0,5000), main = "", ylab = "", xlab = "", col = wes_palette("FantasticFox1")[2])
h2.tax <- hist(log1p(as.matrix(taxa.df2)), freq = TRUE, breaks = seq(0,20,0.5), 
     xlim = c(0,20), ylim = c(0,5000), main = "", ylab = "", xlab = "", col = wes_palette("FantasticFox1")[3], add = T)
box()
dev.off()
# Hist 2
jpeg("hist2.tax.jpg", width = 4.25, height = 4, units = "in", res = 300)
par(mar = c(2,2,0.5,0.5))
h2.tax <- hist(log1p(as.matrix(taxa.df2)), freq = TRUE, breaks = seq(0,20,0.5), 
     xlim = c(0,20), ylim = c(0,5000), main = "", ylab = "", xlab = "", col = wes_palette("FantasticFox1")[3])
h3.tax <- hist(log1p(as.matrix(taxa.df3)), freq = TRUE, breaks = seq(0,20,0.5), 
     xlim = c(0,20), ylim = c(0,5000), main = "", ylab = "", xlab = "", col = wes_palette("FantasticFox1")[5], add = T)
box()
dev.off()
# Hist 3
jpeg("hist3.tax.jpg", width = 4.25, height = 4, units = "in", res = 300)
par(mar = c(2,2,0.5,0.5))
h3.tax <- hist(log1p(as.matrix(taxa.df3)), freq = TRUE, breaks = seq(0,20,0.5), 
     xlim = c(0,20), ylim = c(0,5000), main = "", ylab = "", xlab = "", col = wes_palette("FantasticFox1")[5])
h4.tax <- hist(log1p(as.matrix(taxa.df4)), freq = TRUE, breaks = seq(0,20,0.5), 
     xlim = c(0,20), ylim = c(0,5000), main = "", ylab = "", xlab = "", col = wes_palette("Rushmore1")[3], add = T)
box()
dev.off()

zero.taxa.bin <- rbind.data.frame(h0.tax$counts[1],h1.tax$counts[1],
                                  h2.tax$counts[1],h3.tax$counts[1],
                                  h4.tax$counts[1])
names(zero.taxa.bin) <- "zero.bin"

## Make the dfs for the scatterplots

# KEGG pathways
kegg.rel.df0 <- rowMeans(decostand(kegg.tpm.df0, MARGIN = 2, "total")*100)
kegg.rel.df1 <- kegg.rel.df0[rownames(kegg.tpm.df1)]
kegg.rel.df2 <- kegg.rel.df0[rownames(kegg.tpm.df2)]
kegg.rel.df3 <- kegg.rel.df0[rownames(kegg.tpm.df3)]
kegg.rel.df4 <- kegg.rel.df0[rownames(kegg.tpm.df4)]
kegg.prev.df0 <- rowSums(kegg.tpm.df0 > 0)
kegg.prev.df1 <- kegg.prev.df0[rownames(kegg.tpm.df1)]
kegg.prev.df2 <- kegg.prev.df0[rownames(kegg.tpm.df2)]
kegg.prev.df3 <- kegg.prev.df0[rownames(kegg.tpm.df3)]
kegg.prev.df4 <- kegg.prev.df0[rownames(kegg.tpm.df4)]
kegg.rp.df0 <- cbind.data.frame(prevalence = kegg.prev.df0, abundance = kegg.rel.df0)
kegg.rp.df1 <- cbind.data.frame(prevalence = kegg.prev.df1, abundance = kegg.rel.df1)
kegg.rp.df2 <- cbind.data.frame(prevalence = kegg.prev.df2, abundance = kegg.rel.df2)
kegg.rp.df3 <- cbind.data.frame(prevalence = kegg.prev.df3, abundance = kegg.rel.df3)
kegg.rp.df4 <- cbind.data.frame(prevalence = kegg.prev.df4, abundance = kegg.rel.df4)

# Taxonomy
taxa.rel.df0 <- rowMeans(decostand(taxa.df0, MARGIN = 2, "total")*100)
taxa.rel.df1 <- taxa.rel.df0[rownames(taxa.df1)]
taxa.rel.df2 <- taxa.rel.df0[rownames(taxa.df2)]
taxa.rel.df3 <- taxa.rel.df0[rownames(taxa.df3)]
taxa.rel.df4 <- taxa.rel.df0[rownames(taxa.df4)]
taxa.prev.df0 <- rowSums(taxa.df0 > 0)
taxa.prev.df1 <- taxa.prev.df0[rownames(taxa.df1)]
taxa.prev.df2 <- taxa.prev.df0[rownames(taxa.df2)]
taxa.prev.df3 <- taxa.prev.df0[rownames(taxa.df3)]
taxa.prev.df4 <- taxa.prev.df0[rownames(taxa.df4)]
tax.rp.df0 <- cbind.data.frame(prevalence = taxa.prev.df0, abundance = taxa.rel.df0)
tax.rp.df1 <- cbind.data.frame(prevalence = taxa.prev.df1, abundance = taxa.rel.df1)
tax.rp.df2 <- cbind.data.frame(prevalence = taxa.prev.df2, abundance = taxa.rel.df2)
tax.rp.df3 <- cbind.data.frame(prevalence = taxa.prev.df3, abundance = taxa.rel.df3)
tax.rp.df4 <- cbind.data.frame(prevalence = taxa.prev.df4, abundance = taxa.rel.df4)

# KEGG scatt0
jpeg("scatt0.kegg.jpg", width = 4.25, height = 4, units = "in", res = 300)
par(mar = c(2,2,0.5,0.5))
plot(abundance ~ jitter(prevalence, 0.5), kegg.rp.df0, log = "y", pch = 16, col = alpha(wes_palette("FantasticFox1")[1],0.9), ylim = c(0.0000001,100), cex = 0.7)
points(jitter(kegg.rp.df1$prevalence,0.5), kegg.rp.df1$abundance, pch = 21, bg = alpha(wes_palette("FantasticFox1")[2],1), lwd = 0.65, cex = 0.7)
dev.off()
# KEGG scatt1
jpeg("scatt1.kegg.jpg", width = 4.25, height = 4, units = "in", res = 300)
par(mar = c(2,2,0.5,0.5))
plot(abundance ~ jitter(prevalence, 0.5), kegg.rp.df1, log = "y", pch = 16, col = alpha(wes_palette("FantasticFox1")[2],0.9), ylim = c(0.0000001,100), cex = 0.7)
points(jitter(kegg.rp.df2$prevalence,0.5), kegg.rp.df2$abundance, pch = 21, bg = alpha(wes_palette("FantasticFox1")[3],1), lwd = 0.65, cex = 0.7)
dev.off()
# KEGG scatt2
jpeg("scatt2.kegg.jpg", width = 4.25, height = 4, units = "in", res = 300)
par(mar = c(2,2,0.5,0.5))
plot(abundance ~ jitter(prevalence, 0.25), kegg.rp.df2, log = "y", pch = 16, col = alpha(wes_palette("FantasticFox1")[3],0.9), ylim = c(0.0000001,100), cex = 0.7)
points(jitter(kegg.rp.df3$prevalence,0.25), kegg.rp.df3$abundance, pch = 21, bg = alpha(wes_palette("FantasticFox1")[5],1), lwd = 0.65, cex = 0.7)
dev.off()
# KEGG scatt3
jpeg("scatt3.kegg.jpg", width = 4.25, height = 4, units = "in", res = 300)
par(mar = c(2,2,0.5,0.5))
plot(abundance ~ jitter(prevalence, 0.25), kegg.rp.df3, log = "y", pch = 16, col = alpha(wes_palette("FantasticFox1")[5],0.9), ylim = c(0.0000001,100), cex = 0.7)
points(jitter(kegg.rp.df4$prevalence,0.25), kegg.rp.df4$abundance, pch = 21, bg = alpha(wes_palette("Rushmore1")[3],1), lwd = 0.65, cex = 0.7)
dev.off()

# KEGG taxa0
jpeg("scatt0.taxa.jpg", width = 4.25, height = 4, units = "in", res = 300)
par(mar = c(2,2,0.5,0.5))
plot(abundance ~ jitter(prevalence, 0.25), tax.rp.df0, log = "y", pch = 16, col = alpha(wes_palette("FantasticFox1")[1],0.9), xlim = c(1,12), ylim = c(0.0000004,15), cex = 0.7)
points(jitter(tax.rp.df1$prevalence,0.25), tax.rp.df1$abundance, pch = 21, bg = alpha(wes_palette("FantasticFox1")[2],1), lwd = 0.65, cex = 0.7)
dev.off()
# KEGG taxa1
jpeg("scatt1.taxa.jpg", width = 4.25, height = 4, units = "in", res = 300)
par(mar = c(2,2,0.5,0.5))
plot(abundance ~ jitter(prevalence, 0.25), tax.rp.df1, log = "y", pch = 16, col = alpha(wes_palette("FantasticFox1")[2],0.9), xlim = c(1,12), ylim = c(0.0000004,15), cex = 0.7)
points(jitter(tax.rp.df2$prevalence,0.25), tax.rp.df2$abundance, pch = 21, bg = alpha(wes_palette("FantasticFox1")[3],1), lwd = 0.65, cex = 0.7)
dev.off()
# KEGG taxa2
jpeg("scatt2.taxa.jpg", width = 4.25, height = 4, units = "in", res = 300)
par(mar = c(2,2,0.5,0.5))
plot(abundance ~ jitter(prevalence, 0.25), tax.rp.df2, log = "y", pch = 16, col = alpha(wes_palette("FantasticFox1")[3],0.9), xlim = c(1,12), ylim = c(0.0000004,15), cex = 0.7)
points(jitter(tax.rp.df3$prevalence,0.25), tax.rp.df3$abundance, pch = 21, bg = alpha(wes_palette("FantasticFox1")[5],1), lwd = 0.65, cex = 0.7)
dev.off()
# KEGG taxa3
jpeg("scatt3.taxa.jpg", width = 4.25, height = 4, units = "in", res = 300)
par(mar = c(2,2,0.5,0.5))
plot(abundance ~ jitter(prevalence, 0.25), tax.rp.df3, log = "y", pch = 16, col = alpha(wes_palette("FantasticFox1")[5],0.9), xlim = c(1,12), ylim = c(0.0000004,15), cex = 0.7)
points(jitter(tax.rp.df4$prevalence,0.25), tax.rp.df4$abundance, pch = 21, bg = alpha(wes_palette("Rushmore1")[3],1), lwd = 0.65, cex = 0.7)
dev.off()