library(PLNmodels)
library(ggplot2)
# args = commandArgs(trailingOnly=TRUE)

rm(list = ls())

##### Load our dataset #########################################################

# abundance <- read.csv(file = args[1],row.names=1,header=TRUE)
abundance <- read.csv("~/Dropbox/Cyclotron Research/Carbon_Oxidation_and_Fixation_Study/Metagenome/SteveM/Curated_Rcode/abu.species.glom1.csv", row.names = 2)[-1]
taxa_names = list(row.names(abundance))
abundance = t(abundance)
covariates <- read.csv("~/Dropbox/Cyclotron Research/Carbon_Oxidation_and_Fixation_Study/Metagenome/SteveM/Network/Transformed_Covariates.csv",row.names=1,header=TRUE)
identical(rownames(abundance), rownames(covariates))
rownames(covariates) <- gsub(" ","", rownames(covariates))
abundance <- abundance[rownames(covariates),]
identical(rownames(abundance), rownames(covariates))

dataset <- prepare_data(abundance,covariates)

############################ Training the network on data ######################

network_models <- PLNnetwork(Abundance ~ 1 + offset(log(Offset)), data = dataset)

##################### Plotting the model ######################################

model_StARS <- getBestModel(network_models, "StARS")
plot(model_StARS, output = "igraph")

model_StARS.igraph <- as.matrix(plot(model_StARS, output = "corrplot"))
sum(model_StARS.igraph[rownames(model_StARS.igraph) == "Ralstonia pickettii",] != 0)
hist(sort(rowSums(model_StARS.igraph != 0)), breaks = seq(0,40,1))
model_StARS.degree <- data.frame(sort(rowSums(model_StARS.igraph != 0)))
names(model_StARS.degree) <- "degree"
write.csv(model_StARS.igraph, "~/Dropbox/Cyclotron Research/Carbon_Oxidation_and_Fixation_Study/Metagenome/SteveM/Curated_Rcode/model_StARS.igraph.csv")
write.csv(model_StARS.degree, "~/Dropbox/Cyclotron Research/Carbon_Oxidation_and_Fixation_Study/Metagenome/SteveM/Curated_Rcode/model_StARS.degree.csv")


## Coneccted to Ralstonia
# Acidobacterium capsulatum                    31
# Burkholderia gladioli                        13
# Corallococcus coralloides                     6
# Cutibacterium acnes                          16
# Cutibacterium granulosum                     15
# Lawsonella clevelandensis                    14
# Nocardia farcinica                           20
# Paraburkholderia caffeinilytica              21
# Pseudomonas yamanorum                        14
# Staphylococcus capitis                       15
# Staphylococcus epidermidis                   14
# Unclassified Acidovorax                      22
# Unclassified Bacillales                      19
# Unclassified Bacilli                         14
# Unclassified Bacteroidetes                    7
# Unclassified Corynebacterium                 18
# Unclassified Cutibacterium                   17
# Unclassified Verrucomicrobia                 23
# Unclassified Gemmatimonadetes                14


data.frame(
  fitted   = as.vector(fitted(model_StARS)),
  observed = as.vector(dataset$Abundance)
) %>% 
  ggplot(aes(x = observed, y = fitted)) + 
  geom_point(size = .5, alpha =.25 ) + 
  scale_x_log10(limits = c(1,1000)) + 
  scale_y_log10(limits = c(1,1000)) + 
  theme_bw() + annotation_logticks()
