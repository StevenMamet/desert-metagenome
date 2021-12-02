##*********************
## This first block of code is for importing, subsetting, and QA/QC of the MG data from SqueezeMeta on the barley server.

#### See here for documentation on jupyterlab 
# # MacOS:
# https://jupyter.org/install.html
# 
# # Windows:
# https://mycarta.wordpress.com/2019/07/09/from-zero-to-jupyterlab-pro-on-windows-10/
# 
# ## Install SqueezeMeta and SQMtools on MacOS or Linux

# # Step 1: Log onto the server if this is where you are running the MG sequences
# ssh nsid@barley.usask.ca # Your NSID
# 
# # Step 2: Install and update conda if necessary. In terminal:
# ## Linux
## This should all be done within the terminal in JupyterLab
# cd $HOME
# wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh
# bash ~/miniconda.sh -b -p $HOME/miniconda
# export PATH="$HOME/miniconda/bin:$PATH"
# conda update conda

# ## MacOSX
# cd $HOME
# wget https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.sh -O ~/miniconda.sh
# bash  ~/miniconda.sh -b -p $HOME/miniconda
# export PATH="$HOME/miniconda/bin:$PATH"
# conda update conda
#
# # Step 3: Install SqueezeMeta:
# conda create -n SqueezeMeta -c bioconda -c conda-forge -c fpusan squeezemeta
# 
# # Activating/deactivating SqueezeMeta:
# source activate SqueezeMeta # May also be: conda activate SqueezeMeta
# source deactivate SqueezeMeta # May also be: conda deactivate SqueezeMeta
# 
# # Step 4: Test installation without database integration
# miniconda/envs/SqueezeMeta/SqueezeMeta/utils/install_utils/test_install.pl
# 
# # Step 5: Direct SqueezeMeta to the database loaded on /u1 # Note this only works on Barley at present (8-Mar-2021)
# ~/miniconda/envs/SqueezeMeta/SqueezeMeta/utils/install_utils/configure_nodb.pl /u1/SqueezeMeta/SqueezeMeta_db/db
# 
# # Step 6: Retest the installation and there should be no error messages this time
# miniconda/envs/SqueezeMeta/SqueezeMeta/utils/install_utils/test_install.pl
#
# # Step 7: Install SQMtools in R. First, in R:
# lapply(names(sessionInfo()$loadedOnly), require, character.only = TRUE)
# invisible(lapply(paste0('package:', names(sessionInfo()$otherPkgs)), detach, character.only=TRUE, unload=TRUE, force=TRUE))
# 
# .rs.restartR()
# 
# sessionInfo()
# 
# 
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("pathview")
# # Then, in terminal find where SQMtools is after installing
# # SqueezeMeta and run this command, but substituting your directory:
# sudo R CMD INSTALL /Users/sdm231/miniconda3/pkgs/squeezemeta-1.3.0-pl526r36_0/SqueezeMeta/lib/SQMtools
# R CMD INSTALL /miniconda/pkgs/squeezemeta-1.3.1-pl526r36_0/SqueezeMeta/lib/SQMtools
#
# Now restart R studio and you should be good to go.

# install.packages("wesanderson") # If necessary

library(reticulate) # Bring in Python
library(SQMtools)   # Required to read in the SqueezeMeta database
library(tidyverse)  # Because it's amazing
library(scales)     # Transparency in plotting
library(vegan)      # Used for rarefaction curves
library(parallel)   # Used to parallelized the rarefaction
library(wesanderson)

# reticulate::use_condaenv("SqueezeMeta") # Configure which version of Python to use

## Read in the desert database (17-May-2021: Currently these data are on the barley server)
## I compressed the sample files here: /u1/SqueezeMeta/SqueezeMeta_Project/Combined_SDS_Analysis/data/sam/
## If there is an error loading, you may need to decompress these files

desert <- loadSQM("/u1/SqueezeMeta/SqueezeMeta_Project/Combined_SDS_Analysis")

#### 1. QA/QC

##*********************
### i. Filter first by taxonomy, then by completion
# To track down the sequences associated with the taxonomy of the final dataset,
# we need to link the taxonomy labels (megahit_XXXX) to a sequence file 

# Taxonomy is in the misc list of the SQM object, which lives here (just FYI):
# /u1/SqueezeMeta/SqueezeMeta_Project/Combined_SDS_Analysis/results/tables

# Subset the full object into only the three kingdoms we're interested in.
full.sub0 <- subsetTax(desert, "superkingdom", c("Bacteria","Archaea","Eukaryota"))

# Then subset to only high-quality bins and the contigs and ORFs contained within them
topBinNames <- rownames(full.sub0$bins$table[full.sub0$bins$table$Completeness > 50 & full.sub0$bins$table$Contamination < 10,]) # n = 49
full.sub5 <- subsetBins(full.sub0, topBinNames)

nrow(full.sub5$taxa$species$abund); print("of"); nrow(desert$taxa$species$abund) # 193 of 4428 taxa preserved
nlevels(as.factor(rownames(full.sub0$taxa$phylum$abund)))  # 148 phyla
nlevels(as.factor(rownames(full.sub5$taxa$phylum$abund)))  # 24 phyla
nlevels(as.factor(rownames(full.sub5$taxa$superkingdom$abund)))  # 3 kingdoms
nlevels(as.factor(rownames(full.sub5$taxa$species$abund)))  # 193 species

full.sub5$taxa$class$abund[grep("Gamma", rownames(full.sub5$taxa$class$abund)),]      # 1 Gamma
full.sub5$taxa$class$abund[grep("Alpha", rownames(full.sub5$taxa$class$abund)),]      # 1 Alpha
full.sub5$taxa$genus$abund[grep("Methylo", rownames(full.sub5$taxa$genus$abund)),]    # 1 Methylococcales
full.sub5$taxa$genus$abund[grep("Crenothrix", rownames(full.sub5$taxa$genus$abund)),]   # 0
full.sub5$taxa$genus$abund[grep("V4", rownames(full.sub5$taxa$genus$abund)),]           # 0
full.sub5$taxa$genus$abund[grep("kam1", rownames(full.sub5$taxa$genus$abund)),]         # 0
full.sub5$taxa$genus$abund[grep("SolV", rownames(full.sub5$taxa$genus$abund)),]         # 0

# Here is our taxonomy table we can link to the tree later
full.tax <- as.data.frame(full.sub5$contigs$tax, stringsAsFactors = F)

# Table that only includes species
contig.tax.df0 <- as.data.frame(full.sub5$contigs$tax[,7], stringsAsFactors = F)
contig.tax.df <- contig.tax.df0
names(contig.tax.df) <- "species"

# There are multiple megahit ids for the same species. Here I just select one for each species.
# 7089 megahit ids, 193 unique species
contig.tax.df2 <- contig.tax.df
contig.tax.df2$mega <- rownames(contig.tax.df)
contig.tax.df3 <- contig.tax.df2 %>%
  group_by(species) %>%
  sample_n(1)
mega.sub <- contig.tax.df3$mega

write.csv(full.tax, "/u1/SqueezeMeta/R_analyses_SDM/full.tax1.csv")
write.csv(mega.sub, "/u1/SqueezeMeta/R_analyses_SDM/mega.sub1.csv")

##*********************
#### 2. Filtering

## KEGG pathways

# Full dataset
desert.df <- as.data.frame(desert$functions$KEGG$abund)
desert.tpm.df <- as.data.frame(desert$functions$KEGG$tpm)
dim(desert.df) # 14017 x 12
# Taxonomy subset
full.df0 <- as.data.frame(full.sub0$functions$KEGG$abund)
full.tpm.df0 <- as.data.frame(full.sub0$functions$KEGG$tpm)
dim(full.df0) # 10541 x 12
# Completion filtered
full.df1 <- as.data.frame(full.sub5$functions$KEGG$abund)
full.tpm.df1 <- as.data.frame(full.sub5$functions$KEGG$tpm)
dim(full.df1) # 4551 x 12

## Step 1: Noise filter
# Remove any KEGG pathways with total # reads < 2
full.df2 <- full.df1[rowSums(full.df1) > 2,]
full.tpm.df2 <- full.tpm.df1[rowSums(full.tpm.df1) > 2,]
dim(full.df2)       # 4551 x 12
dim(full.tpm.df2)   # 4551 x 12
# Remove zero-sum samples 
full.df2 <- full.df2[,colSums(full.df2) > 2]
full.tpm.df2 <- full.tpm.df2[,colSums(full.tpm.df2) > 2]
dim(full.df2)       # 4551 x 12
dim(full.tpm.df2)   # 4551 x 12

## Step 2: Trace filter
# Remove any KEGG prevalence singltons/doubletons
full.df3 <- full.df2[rowSums(full.df2 > 0) > 2,]
full.tpm.df3 <- full.tpm.df2[rowSums(full.tpm.df2 > 0) > 2,]
dim(full.df3)       # 4408 x 12
dim(full.tpm.df3)   # 4408 x 12
# Remove zero-sum samples 
full.df3 <- full.df3[,colSums(full.df3 != 0) > 0]
full.tpm.df3 <- full.tpm.df3[,colSums(full.tpm.df3 != 0) > 0]
dim(full.df3)       # 4408 x 12
dim(full.tpm.df3)   # 4408 x 12

## Note to self:
# SQM.df$$taxa$species$abund is the same as summing SQM.df$contigs$abund by species

## Taxa abundance
# Full dataset
desert.tax.df <- as.data.frame(desert$contigs$tax)
desert.tax.df <- as.data.frame(mutate_if(desert.tax.df, is.factor, as.character))
desert.tax.df[rownames(desert.tax.df) %in% contig.tax.df3$mega,]
desert.tax.df0 <- desert.tax.df[desert.tax.df$species %in% rownames(desert$taxa$species$abund),]
desert.abu.df <- as.data.frame(desert$taxa$species$abund)
dim(desert.abu.df) # 4428 x 12 (3,623,269 x 12)

# Taxonomy subset
full.abu.df0 <- as.data.frame(full.sub0$taxa$species$abund)
full.abu.df0 <- full.abu.df0[rowSums(full.abu.df0) != 0,]
dim(full.abu.df0) # 3003 x 12 (925050 x 12)

# Completion filtered
full.abu.df1 <- as.data.frame(full.sub5$taxa$species$abund)
dim(full.abu.df1) # 193 x 12 (7089 x 12)

## Step 1: Noise filter
# Remove any taxa with total # reads < 2
full.abu.df2 <- full.abu.df1[rowSums(full.abu.df1) > 2,]
dim(full.abu.df2) # 193 x 12 (7089 x 12)

# Remove zero-sum samples 
full.abu.df2 <- full.abu.df2[,colSums(full.abu.df2) > 2]
dim(full.abu.df2) # 193 x 12 (7089 x 12)

## Step 2: Trace filter
# Remove any KEGG prevalence singltons/doubletons
full.abu.df3 <- full.abu.df2[rowSums(full.abu.df2 > 0) > 2,]
dim(full.abu.df3) # 183 x 12 (6805 x 12)

# Remove zero-sum samples 
full.abu.df3 <- full.abu.df3[,colSums(full.abu.df3 != 0) > 0]
dim(full.abu.df3) # 183 x 12 (6805 x 12)

filter.df <- cbind.data.frame(filter = c("SqueezeMeta","tax","completion","noise","trace"),
                              id = seq(1,5),
                              KEGG = c(dim(desert.tpm.df)[1],
                                       dim(full.tpm.df0)[1],
                                       dim(full.tpm.df1)[1],
                                       dim(full.tpm.df2)[1],
                                       dim(full.tpm.df3)[1]),
                              abundance = c(dim(desert.abu.df)[1],
                                            dim(full.abu.df0)[1],
                                            dim(full.abu.df1)[1],
                                            dim(full.abu.df2)[1],
                                            dim(full.abu.df3)[1]))
## Note: using TPM or just abundance doesn't change the prevalence

(filter.df$KEGG[2] / filter.df$KEGG[1])*0.9 # 0.677
(filter.df$KEGG[3] / filter.df$KEGG[1])*0.9 # 0.292
(filter.df$KEGG[4] / filter.df$KEGG[1])*0.9 # 0.292
(filter.df$KEGG[5] / filter.df$KEGG[1])*0.9 # 0.283
(filter.df$abundance[2] / filter.df$abundance[1])*0.9 # 0.610
(filter.df$abundance[3] / filter.df$abundance[1])*0.9 # 0.039
(filter.df$abundance[4] / filter.df$abundance[1])*0.9 # 0.039
(filter.df$abundance[5] / filter.df$abundance[1])*0.9 # 0.037

write.csv(filter.df, "/u1/SqueezeMeta/R_analyses_SDM/filter.barplot.csv")
write.csv(desert.df, "/u1/SqueezeMeta/R_analyses_SDM/kegg.full.csv")
write.csv(full.df0, "/u1/SqueezeMeta/R_analyses_SDM/kegg.filt0.csv")
write.csv(full.df1, "/u1/SqueezeMeta/R_analyses_SDM/kegg.filt1.csv")
write.csv(full.df2, "/u1/SqueezeMeta/R_analyses_SDM/kegg.filt2.csv")
write.csv(full.df3, "/u1/SqueezeMeta/R_analyses_SDM/kegg.filt3.csv")

write.csv(desert.tpm.df, "/u1/SqueezeMeta/R_analyses_SDM/kegg.full.tpm.csv")
write.csv(full.tpm.df0, "/u1/SqueezeMeta/R_analyses_SDM/kegg.filt0.tpm.csv")
write.csv(full.tpm.df1, "/u1/SqueezeMeta/R_analyses_SDM/kegg.filt1.tpm.csv")
write.csv(full.tpm.df2, "/u1/SqueezeMeta/R_analyses_SDM/kegg.filt2.tpm.csv")
write.csv(full.tpm.df3, "/u1/SqueezeMeta/R_analyses_SDM/kegg.filt3.tpm.csv")

write.csv(desert.tax.df, "/u1/SqueezeMeta/R_analyses_SDM/taxonomy.full.csv")
write.csv(desert.abu.df, "/u1/SqueezeMeta/R_analyses_SDM/taxa.full.csv")
write.csv(full.abu.df0, "/u1/SqueezeMeta/R_analyses_SDM/taxa.filt0.csv")
write.csv(full.abu.df1, "/u1/SqueezeMeta/R_analyses_SDM/taxa.filt1.csv")
write.csv(full.abu.df2, "/u1/SqueezeMeta/R_analyses_SDM/taxa.filt2.csv")
write.csv(full.abu.df3, "/u1/SqueezeMeta/R_analyses_SDM/taxa.filt3.csv")

# Rarefaction curve tells you about the rate at which new taxa are detected
# as you increase the number of sequences sampled. It does this by taking
# random subsamples from 1, up to the size of your sample, and computing the
# number of taxa present in each subsample. Ideally you would like
# you rarefaction curves to be relatively flat as this indicates that
# additional sampling would not likely yield further taxa.

# Here's a parallelized version of the 'rarecurve' function in vegan.
# Borrowed from here: https://dave-clark.github.io/post/speeding-up-rarefaction-curves-for-microbial-community-ecology/

## Note: even on the server, this will take a few days to run
quickRareCurve <- function (x, step = 1, sample, xlab = "Sample Size",
                            ylab = "Species", label = TRUE, col, lty, max.cores = T, nCores = 1, ...)
{
  require(parallel)
  x <- as.matrix(x)
  if (!identical(all.equal(x, round(x)), TRUE))
    stop("function accepts only integers (counts)")
  if (missing(col))
    col <- par("col")
  if (missing(lty))
    lty <- par("lty")
  tot <- rowSums(x) # calculates library sizes
  S <- specnumber(x) # calculates n species for each sample
  if (any(S <= 0)) {
    message("empty rows removed")
    x <- x[S > 0, , drop = FALSE]
    tot <- tot[S > 0]
    S <- S[S > 0]
  } # removes any empty rows
  nr <- nrow(x) # number of samples
  col <- rep(col, length.out = nr)
  lty <- rep(lty, length.out = nr)
  # parallel mclapply
  # set number of cores
  mc <- getOption("mc.cores", ifelse(max.cores, detectCores(), nCores))
  message(paste("Using ", mc, " cores"))
  out <- mclapply(seq_len(nr), mc.cores = mc, function(i) {
    n <- seq(1, tot[i], by = step)
    if (n[length(n)] != tot[i])
      n <- c(n, tot[i])
    drop(rarefy(x[i, ], n))
  })
  Nmax <- sapply(out, function(x) max(attr(x, "Subsample")))
  Smax <- sapply(out, max)
  plot(c(1, max(Nmax)), c(1, max(Smax)), xlab = xlab, ylab = ylab,
       type = "n", ...)
  if (!missing(sample)) {
    abline(v = sample)
    rare <- sapply(out, function(z) approx(x = attr(z, "Subsample"),
                                           y = z, xout = sample, rule = 1)$y)
    abline(h = rare, lwd = 0.5)
  }
  for (ln in seq_along(out)) {
    N <- attr(out[[ln]], "Subsample")
    lines(N, out[[ln]], col = col[ln], lty = lty[ln], ...)
  }
  if (label) {
    ordilabel(cbind(tot, S), labels = rownames(x), ...)
  }
  invisible(out)
}

# Look at rarefaction for the KEGG gene abundances
rare.df <- quickRareCurve(t(desert$functions$KEGG$abund), max.cores = F, nCores = 60, step = 5, sample = 100) # use 60 cores
abline(v=(min(rowSums(t(desert$functions$KEGG$abund)))))

# Look at rarefaction for the taxa gene abundances
rare.abund.df <- quickRareCurve(t(desert$taxa$species$abund), max.cores = F, nCores = 60, step = 5, sample = 100) # use 60 cores
abline(v=(min(rowSums(t(desert$taxa$species$abund)))))

col.palw <- wes_palette("Rushmore1", 12, type = "continuous")

jpeg("rarify.jpg", width = 8, height = 4, units = "in", res = 300)
par(mfrow = c(1,2), mar = c(1,2.5,0.1,0.1), oma = c(2.5, 2.5, 0.1, 0.1))
plot.new()
plot.window(xlim = c(0,2.1e7), ylim = c(0,12000))
lines(rare.df[[1]], col = col.palw[9])
lines(rare.df[[2]], col = col.palw[8])
lines(rare.df[[3]], col = col.palw[2])
lines(rare.df[[4]], col = col.palw[1])
lines(rare.df[[5]], col = col.palw[5])
lines(rare.df[[6]], col = col.palw[6])
lines(rare.df[[7]], col = col.palw[3])
lines(rare.df[[8]], col = col.palw[7])
lines(rare.df[[9]], col = col.palw[10])
lines(rare.df[[10]], col = col.palw[11])
lines(rare.df[[11]], col = col.palw[12])
lines(rare.df[[12]], col = col.palw[4])
box()
axis(side = 1, at = seq(0,3e7,0.5e7))
axis(side = 2, at = seq(0,12000,2000))

plot.new()
plot.window(xlim = c(0,2.1e7), ylim = c(0,3100))
lines(rare.abund.df[[1]], col = col.palw[9])
lines(rare.abund.df[[2]], col = col.palw[8])
lines(rare.abund.df[[3]], col = col.palw[2])
lines(rare.abund.df[[4]], col = col.palw[1])
lines(rare.abund.df[[5]], col = col.palw[5])
lines(rare.abund.df[[6]], col = col.palw[6])
lines(rare.abund.df[[7]], col = col.palw[3])
lines(rare.abund.df[[8]], col = col.palw[7])
lines(rare.abund.df[[9]], col = col.palw[10])
lines(rare.abund.df[[10]], col = col.palw[11])
lines(rare.abund.df[[11]], col = col.palw[12])
lines(rare.abund.df[[12]], col = col.palw[4])
box()
axis(side = 1, at = seq(0,3e7,0.5e7))
axis(side = 2, at = seq(0,4000,1000))
dev.off()