## Grouping co-occurring OTUs into modules using the WGCNA package

# Subset OTU_main_LC_filtered to keep only the columns listed in Keep$SampleID
OTU_main_LC_filtered_ok <- OTU_main_LC_filtered[, colnames(OTU_main_LC_filtered) %in% Keep$SampleID]

# Keep only the columns in OTU_main_LC that match SampleID in dataset
OTU_main_LC_filtered <- OTU_main_LC[, colnames(OTU_main_LC) %in% dataset$SampleID]
# Filter TAX_table_bacteria to keep only rows where zOTU is in the row names of OTU_main_LC_filtered
TAX_fungi_filtered <- TAX_fungi[row.names(TAX_fungi) %in% rownames(OTU_main_LC_filtered), ]

# Convert all columns to numeric
seq_table[] <- lapply(seq_table, function(x) as.numeric(as.character(x)))

library(microeco)
library(dplyr)

# Let's create a microtable object with more information
mt <- microtable$new(otu_table = seq_table, tax_table = tax_table)
mt

mt$tidy_dataset()
mt

# first clone the data
mt_rarefied <- clone(mt)
# use sample_sums to check the sequence numbers in each sample
mt_rarefied$sample_sums() %>% range

# As an example, use 10000 sequences in each sample
mt_rarefied$rarefy_samples(sample.size = 1100)

#WGCNA The Bowman Lab tutorial

Taxa_table=as.data.frame(mt_rarefied$otu_table)

# Assuming Taxa_table is a data frame or matrix where rows are taxa and columns are samples
Taxa_table_rel_abund <- sweep(Taxa_table, 2, colSums(Taxa_table), FUN = "/")

# Filter OTUs that are present (relative abundance > 0) in at least 10% samples
datExpr0 <- Taxa_table_rel_abund[rowSums(Taxa_table_rel_abund > 0) >= 10, ]
datExpr0 = as.data.frame(datExpr0)

# Calculate mean relative abundance per OTU (row)
mean_abundance <- mean(colSums(datExpr0))

# Filter OTUs with mean relative abundance >= 0.0001
Taxa_table_final <- Taxa_table_filtered[mean_abundance >= 0.00001, ]
average_total_abundance <- mean(rowMeans(datExpr0))
write.csv(Taxa_table_final,file="taxa_table_filtered.csv")


#change working directory if needed
data<-read.table("MB.0.03.subsample.fn.txt",header=T,na.strings="NA")
#Remove all columns with no OTU data
data1 = data[-1][-1][-1]

#Use Hellinger Transformation to calculate (square root of) relative abundance
library(vegan)
HellingerData<-decostand(data1,method = "hellinger")

#Sort OTUs by abundance, since you have to limit the OTUs to the most frequent ones. You need to decide here!
lessdata <- data1[,colSums(data1) > 0.05]

#Making your relative abundance matrix
RelAbun1 = data.frame(data[2],HellingerData[1:750])
write.table(RelAbun1, file = "MontereyRelAbun.txt", sep="\t")

#Could I start from here if I already have my rel abu table from microeco?
library(WGCNA)
#OTUs<-read.table("MontereyRelAbun.txt",header=T,sep="\t")


datExpr0 = as.data.frame(datExpr0)
#rownames have to be sample names, colnames have to be OTU IDs

#Turn first column into row names so that only OTUs make up actual columns
#use this code or another
datExpr0 = as.data.frame((OTUs[,-c(1)]));
names(datExpr0) = names(OTUs)[-c(1)];
rownames(datExpr0) = OTUs$Group;

#check data for excessive missingness, you should get TRUE if you got rid of all OTUs that have many 0s
#gsg = goodSamplesGenes(datExpr0[-1], verbose = 3)
gsg = goodSamplesGenes(datExpr0, verbose = 3)
gsg$allOK #you should get TRUE

#Cluster the samples to see if there are any obvious outliers, this will cluster the samples to see if some sample is outlier
sampleTree = hclust(dist(datExpr0), method = "average");

sizeGrWindow(12,9)

par(cex = 0.6);
par(mar = c(0,4,2,0))

plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", 
     cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)

#now read in trait data
#First column has to be SampleID then the variables (multifunctionality?)
traitData <- dataset[, c(1, 9:16)]
#traitData = dataset

#create a data frame including OTU rel abu and environmental data
#OTU abundances are stored in datExpr0 and the environmental traits are stored in datTraits
OTUSamples = row.names(datExpr0);
traitRows = match(OTUSamples, traitData$SampleID);
datTraits = traitData[traitRows, -1];
datTraits = as.data.frame(datTraits)
rownames(datTraits) = traitData$SampleID
#rownames(datTraits) = traitData$SampleID
collectGarbage()

sampleTree2 = hclust(dist(datExpr0), method = "average")
#number of objects in sampleTree2 has to be number of SAMPLES! If you have number of OTUs here, you need to transpose datExpr0

traitColors = numbers2colors(datTraits[1:8], signed = FALSE)

#This is very preliminary, this is how samples relate to traits (no modules here...)
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(datTraits[1:8]),
                    main = "Sample dendrogram and trait heatmap")

#save(datExpr0, datTraits, file = "Monterey-dataInput.RData")

#Network analysis
options(stringsAsFactors = FALSE)
enableWGCNAThreads()
#lnames = load(file = "Monterey-dataInput.RData")

#Creating a weighted OTU co-expression network
#choosing a set of soft thresholding powers
powers = c(c(1:10), seq(from = 11, to=30, by=1))

#Call the network topology analysis function
library(doParallel)
registerDoSEQ()  # Forces sequential execution
sft <- pickSoftThreshold(datExpr0, powerVector = powers, verbose = 5, networkType = "signed", blockSize = 1000)
#This will show the Power (soft thresholding value), the R2 for the scale independence for each power (should be >0.8),
#the mean number of connections each node has at each power (mean.k), the median number of connections, and the max num of con.

#PLOT RESULTS
sizeGrWindow(9, 5)
par(mfrow = c(1,2))
cex1 = 0.9

plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red")
abline(h=0.8,col="red")

plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
#based on these plots, we proceed with a soft thresholding value of 10, because >0.8 and mean connectivity >10
#now let's calculate adjacencies, using soft thresholding power of 10
softPower = 8
adjacency = adjacency(datExpr0, power = softPower, type = "signed")

# Then we transform the adjacency matrix into a Topological Overlap Matrix (TOM) 
# and calculate corresponding dissimilarity
TOM = TOMsimilarity(adjacency, TOMType = "signed")
dissTOM = 1-TOM

#Create a dendogram using a hierarchical clustering tree and then call the hierarchical clustering function
TaxaTree = hclust(as.dist(dissTOM), method = "average")
sizeGrWindow(12,9)
plot(TaxaTree, xlab="", sub="", main = "Taxa clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04)

#you need to decide on minimum module size
set.seed(1)
minModuleSize = 15

#Module identification using dynamic tree cut
dynamicMods = cutreeDynamic(dendro = TaxaTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize);
table(dynamicMods)

dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
# below each module (now assigned a color) you can see the number of OTUs that belong to that module

#Plot the dendrogram with module colors underneath
sizeGrWindow(8,6)
plotDendroAndColors(TaxaTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Taxa dendrogram and module colors")

#Now we will quantify co-expression similarity of the entire modules 
#using eigengenes and cluster them based on their correlation:
#An eigengene is 1st principal component of a module expression matrix and represents a 
#suitably defined average OTU community.

#Calculate eigengenes:

MEList = moduleEigengenes(datExpr0, colors = dynamicColors)
MEs = MEList$eigengenes
MEDiss = 1-cor(MEs)

METree = hclust(as.dist(MEDiss), method = "average")

sizeGrWindow(7, 6)
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")

#Now we will see if any of the modules should be merged. 
#I chose a height cut of 0.30, corresponding to a similarity of 0.70 to merge
MEDissThres = 0.30
abline(h=MEDissThres, col = "red")
#This will not do anything if none of the modules should be merged
merge = mergeCloseModules(datExpr0, dynamicColors, cutHeight = MEDissThres, verbose = 3)
mergedColors = merge$colors
mergedMEs = merge$newMEs
sizeGrWindow(12, 9)

plotDendroAndColors(TaxaTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

moduleColors = mergedColors
#Construct numerical labels corresponding to the colors
colorOrder = c("grey", standardColors(50));
moduleLabels = match(moduleColors, colorOrder)-1;
MEs = mergedMEs;
save(MEs, moduleLabels, moduleColors, TaxaTree, file = "Monterey-networkConstruction-stepByStep.RData")

#Relating modules to external information and IDing important taxa
#we will use Eigengenes of each module since these represent a suitable defined average OTU community

#Defining numbers of OTUs and samples:
nTaxa = ncol(datExpr0)
nSamples = nrow(datExpr0)

#Recalculate MEs (module eigengenes)
MEs0 = moduleEigengenes(datExpr0, moduleColors)$eigengenes
MEs = orderMEs(MEs0)

moduleTraitCor = cor(MEs, datTraits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)

#Now we will visualize it
sizeGrWindow(10,6)

textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3))

labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))

#Now imagine that the red module is interesting, you can further explore individual taxa within this module
#First define the variable we are interested in from datTrait

EMF_average = as.data.frame(datTraits$EMF_average);
names(EMF_average) = "EMF_average"

modNames = substring(names(MEs), 3)
TaxaModuleMembership = as.data.frame(cor(datExpr0, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(TaxaModuleMembership), nSamples));
names(TaxaModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");
TaxaTraitSignificance = as.data.frame(cor(datExpr0, EMF_average, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(TaxaTraitSignificance), nSamples));
names(TaxaTraitSignificance) = paste("GS.", names(EMF_average), sep="");
names(GSPvalue) = paste("p.GS.", names(EMF_average), sep="");

module = "blue"
column = match(module, modNames)
moduleTaxa = moduleColors==module

sizeGrWindow(7, 7)
par(mfrow = c(1,1))

verboseScatterplot(abs(TaxaModuleMembership[moduleTaxa, column]),
                   abs(TaxaTraitSignificance[moduleTaxa, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Taxa significance for EMF_average",
                   main = paste("Module membership vs. Taxa significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)

#This graph shows you how each taxa (each red dot is an OTU that 
#belongs in the Red module) correlated with the Environmental trait 
#of interest and how important it is to the module.
#MORE INFO ON THE TUTORIAL

#Annotating the module
#Now lets get more info about the taxa that make up the Red module
#names(datExpr0)
#names(datExpr0)[moduleColors=="blue"]

annot = TAX_fungi_filtered
dim(annot)
names(annot)
probes = colnames(datExpr0)
probes2annot = match(probes, annot$OTU)

sum(is.na(probes2annot))

TaxaInfo0 = data.frame(Taxon = probes,
                       TaxaSymbol = annot$OTU[probes2annot],
                       LinkID = annot$TAX[probes2annot],
                       moduleColor = moduleColors)

modOrder = order(-abs(cor(MEs, CR, use = "p")))

for (mod in 1:ncol(TaxaModuleMembership))
{
  oldNames = names(TaxaInfo0)
  TaxaInfo0 = data.frame(TaxaInfo0, TaxaModuleMembership[, modOrder[mod]],
                         MMPvalue[, modOrder[mod]]);
  names(TaxaInfo0) = c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""),
                       paste("p.MM.", modNames[modOrder[mod]], sep=""))
}

TaxaOrder = order(TaxaInfo0$moduleColor, -abs(TaxaInfo0$GS.CR));
TaxaInfo = TaxaInfo0[TaxaOrder, ]

write.csv(TaxaInfo, file = "TaxaInfo.csv")

library(tidyr)

# Separate LinkID column into six new columns
TaxaInfo0_ok <- TaxaInfo0 %>%
  separate(LinkID, into = c("phylum", "class", "order", "family", "genus", "species"), sep = ";", fill = "right")


# Assume your data is in a data.frame called 'dataset'
# Functions are in columns 10 to 15
function_data <- dataset[, 10:15]


## Calculating ecosystem multifunctionality indices

# ----------------------------
# Multiple-threshold indices
# ----------------------------
max_values <- apply(function_data, 2, max, na.rm = TRUE)

thresholds <- sweep(matrix(rep(max_values, each = 3), nrow = 3, byrow = TRUE),
                    1,
                    c(0.25, 0.50, 0.75),
                    `*`)
rownames(thresholds) <- c("25", "50", "75")

dataset$EMF_25 <- rowSums(sweep(function_data, 2, thresholds["25",], `>=`), na.rm = TRUE)
dataset$EMF_50 <- rowSums(sweep(function_data, 2, thresholds["50",], `>=`), na.rm = TRUE)
dataset$EMF_75 <- rowSums(sweep(function_data, 2, thresholds["75",], `>=`), na.rm = TRUE)

# -----------------------------------------
# Z-score index (mean of standardized Zs)
# -----------------------------------------
z_mat <- scale(function_data)                          # column-wise z-standardization
dataset$EMF_zscore <- rowMeans(z_mat, na.rm = TRUE)

# ----------------------------------------------------
# 0–1 normalized average (min–max per function/column)
# ----------------------------------------------------
col_min <- apply(function_data, 2, min, na.rm = TRUE)
col_rng <- apply(function_data, 2, function(x) diff(range(x, na.rm = TRUE)))
col_rng[col_rng == 0] <- 1                              # avoid division by zero
mm_mat <- sweep(sweep(function_data, 2, col_min, `-`), 2, col_rng, `/`)
dataset$EMF_avg01 <- rowMeans(mm_mat, na.rm = TRUE)

# ----------------------------------------------------------
# Weighted index:
#  - enzymes share 0.25 total (i.e., ~0.0833 each)
#  - respiration = 0.25, aggregation = 0.25, NPP = 0.25
#   (works by names if available; otherwise falls back to positions)
# ----------------------------------------------------------
fn_names <- colnames(function_data)

# Try to detect columns by (common) names
enz_candidates <- c("NAG","NAG_activity","N-acetylglucosaminidase",
                    "Phosphatase","Phosphatase_activity",
                    "Xylosidase","Xylosidase_activity")
resp_candidates <- c("Basal_resp","BasalRespiration","Soil_basal_respiration","Respiration")
npp_candidates  <- c("NPP","NetPrimaryProductivity","PrimaryProductivity","Primary_Prod")
agg_candidates  <- c("Aggregate_stability","Soil_aggregation","Aggregation","Agg")

enz_idx <- which(fn_names %in% enz_candidates)
resp_idx <- which(fn_names %in% resp_candidates)
npp_idx  <- which(fn_names %in% npp_candidates)
agg_idx  <- which(fn_names %in% agg_candidates)

# Fallback to positional assumption if needed (cols 10:15 in this order):
# [1]=Resp, [2:4]=Enzymes, [5]=NPP, [6]=Aggregation
if (length(enz_idx) == 0 || length(resp_idx) == 0 || length(npp_idx) == 0 || length(agg_idx) == 0) {
  resp_idx <- 1
  enz_idx  <- 2:4
  npp_idx  <- 5
  agg_idx  <- 6
}

# Build weights vector aligned to columns
w <- rep(0, ncol(function_data))
# enzymes share 0.25 total
w[enz_idx] <- 0.25 / length(enz_idx)
# each singletons get 0.25
w[resp_idx] <- 0.25
w[npp_idx]  <- 0.25
w[agg_idx]  <- 0.25

# Use 0–1 normalized values for weighting (consistent scale)
weighted_vals <- as.matrix(mm_mat) %*% matrix(w, ncol = 1)
dataset$EMF_weighted <- as.numeric(weighted_vals)

## Variance partition analyses

dataset <- read_excel("working_dataset_140825.xlsx", 
                      sheet = "dataset140825")
#dataset <- subset(dataset, 
#                 dataset$LC_simpl_2018 == "Cropland")

# Function to compute unique and shared variance per land use category

set.seed(1)
compute_variance_partition <- function(data_subset, group_label) {
  Y <- data_subset$RESP_STD
  edaphic    <- data_subset[, 16:22]
  climate    <- data_subset[, 14:15]
  microbiome <- data_subset[, 23:45]
  
  lm_all   <- lm(Y ~ ., data = cbind(edaphic, climate, microbiome))
  lm_C_M   <- lm(Y ~ ., data = cbind(climate, microbiome))
  lm_E_M   <- lm(Y ~ ., data = cbind(edaphic, microbiome))
  lm_E_C   <- lm(Y ~ ., data = cbind(edaphic, climate))
  lm_EC_M  <- lm(Y ~ ., data = cbind(edaphic, climate, microbiome))
  
  unique_edaphic    <- summary(lm_EC_M)$adj.r.squared - summary(lm_C_M)$adj.r.squared
  unique_climate    <- summary(lm_EC_M)$adj.r.squared - summary(lm_E_M)$adj.r.squared
  unique_microbiome <- summary(lm_EC_M)$adj.r.squared - summary(lm_E_C)$adj.r.squared
  shared            <- summary(lm_all)$adj.r.squared - (unique_edaphic + unique_climate + unique_microbiome)
  shared            <- ifelse(shared < 0, 0, shared)
  
  data.frame(
    LandUse = group_label,
    Category = c("Edaphic", "Climate", "Microbiome", "Shared"),
    Variance_Explained = c(unique_edaphic, unique_climate, unique_microbiome, shared) * 100
  )
}

## Random Forest models

dataset <- read_excel("working_dataset_220825.xlsx", 
                      sheet = "dataset140825")

dataset <- subset(dataset, 
                  dataset$PH_Class == "Neutral")

dataset <- dataset[, -c(1:12)]


set.seed(1)
rfPermute = rfPermute(EMF_weighted ~ .,
                      data = dataset, ntree = 1000)
library(caret)
set.seed(1)
ctrl <- trainControl(method = "cv", number = 10)
fit_cv <- train(EMF_weighted ~ ., data = dataset,
                method = "rf",
                trControl = ctrl,
                ntree = 1000,
                metric = "Rsquared")
max(fit_cv$results$Rsquared)  # cross-validated R^2

# Extract final OOB MSE
oob_mse <- rfPermute$rf$mse[length(rfPermute$rf$mse)]
# Convert to numeric just in case
oob_mse <- as.numeric(oob_mse)
# Calculate RMSE
oob_rmse <- sqrt(oob_mse)
oob_rmse

#predictor importance
importance(rfPermute)

## Multi-group Structural Equation Models (SEMs)

#SEM

library(lavaan)
library(dplyr)

dataset <- read_excel("~/Desktop/MINOTAUR-EMF/Revisions_June_2025/Final_results_290625/working_dataset_300625.xlsx", 
                      sheet = "Random_Forest")
#dataset <- subset(dataset, 
#                 dataset$LC_simpl_2018 == "Cropland")
dataset$Total_nitrogen = NULL

set.seed(1)

# Ensure the group variable is a factor
dataset$PH_Class <- as.factor(dataset$PH_Class)

# Standardize all predictors (columns 7 to 37)
dataset[, 7:37] <- scale(dataset[, 7:37])

# Define variable groups
micro_vars   <- names(dataset)[7:29]   # Microbial composition (observed)
soil_vars    <- names(dataset)[30:35]  # Soil properties
climate_vars <- names(dataset)[36:37]  # Climate variables

# 1. Direct effects on EMF_weighted
direct_effects <- paste(
  "EMF_weighted ~", 
  paste(c(micro_vars, soil_vars, climate_vars), collapse = " + ")
)

# 2. Microbes predicted by Soil
indirect_micro_from_soil <- paste(
  sapply(micro_vars, function(mv) paste(mv, "~", paste(soil_vars, collapse = " + "))),
  collapse = "\n"
)

# 3. Microbes predicted by Climate
indirect_micro_from_climate <- paste(
  sapply(micro_vars, function(mv) paste(mv, "~", paste(climate_vars, collapse = " + "))),
  collapse = "\n"
)

# Combine into full model
model <- paste(direct_effects, indirect_micro_from_soil, indirect_micro_from_climate, sep = "\n")

# Fit the multi-group SEM
fit_group <- sem(model, 
                 data = dataset, 
                 group = "PH_Class", 
                 fixed.x = FALSE)

# Show results
summary(fit_group, standardized = TRUE, fit.measures = TRUE)

# Fit a constrained model where all regressions are equal across groups
fit_constrained <- sem(model, 
                       data = dataset, 
                       group = "PH_Class", 
                       group.equal = "regressions", 
                       fixed.x = FALSE)

# Compare models
anova(fit_group, fit_constrained)

# Extract parameter estimates
estimates <- parameterEstimates(fit_group, standardized = TRUE)
# Write to CSV
write.csv(estimates, "SEM_parameter_estimates.csv", row.names = FALSE)
# Extract fit measures
fit_stats <- fitMeasures(fit_group)
# Convert to data frame
fit_stats_df <- as.data.frame((fit_stats))
# Write to CSV
write.csv(fit_stats_df, "SEM_fit_measures.csv", row.names = TRUE)

# Load and extract parameter estimates from SEM
est <- parameterEstimates(fit_group, standardized = TRUE)

# Match lavaan group IDs to land-use labels
group_labels <- data.frame(
  group = 1:3,
  LandUse = levels(dataset$PH_Class)
)

# Define your variable groups
micro_vars   <- names(dataset)[7:29]
soil_vars    <- names(dataset)[30:35]
climate_vars <- names(dataset)[36:37]

direct_effects <- est %>%
  filter(lhs == "EMF_weighted", op == "~") %>%
  mutate(category = case_when(
    rhs %in% micro_vars   ~ "Microbial",
    rhs %in% soil_vars    ~ "Soil",
    rhs %in% climate_vars ~ "Climate",
    TRUE                  ~ NA_character_
  )) %>%
  filter(!is.na(category)) %>%
  group_by(group, category) %>%
  summarise(
    mean_effect = mean(abs(std.all), na.rm = TRUE),
    effect_type = "Direct",
    .groups = "drop"
  )

indirect_effects <- est %>%
  filter(rhs %in% c(soil_vars, climate_vars), lhs %in% micro_vars, op == "~") %>%
  mutate(category = case_when(
    rhs %in% soil_vars    ~ "Soil",
    rhs %in% climate_vars ~ "Climate",
    TRUE                  ~ NA_character_
  )) %>%
  filter(!is.na(category)) %>%
  group_by(group, category) %>%
  summarise(
    mean_effect = mean(abs(std.all), na.rm = TRUE),
    effect_type = "Indirect",
    .groups = "drop"
  )

# Combine direct and indirect effects
summary_effects <- bind_rows(direct_effects, indirect_effects)

# Add land use labels
summary_effects <- left_join(summary_effects, group_labels, by = "group") %>%
  select(LandUse, effect_type, category, mean_effect)

print(summary_effects)
