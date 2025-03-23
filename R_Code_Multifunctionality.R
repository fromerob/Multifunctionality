# Ecosystem Multifunctionality and Kruskal Wallis test
# Step 1: Load the dataset
data <- read.csv("dataset.csv")

# Step 2: Variables to be used for multifunctionality calculation are NAG, APA, XYL, NPP, RESP, WSA

# Step 3: Apply z-score transformation to the variables mentioned in Step 2
data$AP_z <- scale(data$AP, center = TRUE, scale = TRUE)
data$XYL_z <- scale(data$XYL, center = TRUE, scale = TRUE)
data$NAG_z <- scale(data$NAG, center = TRUE, scale = TRUE)
data$RESP_z <- scale(data$RESP, center = TRUE, scale = TRUE)
data$NPP_z <- scale(data$NPP, center = TRUE, scale = TRUE)
data$WSA_z <- scale(data$WSA, center = TRUE, scale = TRUE)

# Step 4: Aggregate (mean) the z-score transformed variables
data$Multifunctionality <- rowMeans(data[, 1:6])

library(dplyr)

filtered_data <- data %>% 
  filter(LC_simpl_2018 %in% c("Cropland", "Grassland", "Woodland"))

kruskal_test <- kruskal.test(Multifunctionality ~ LC_simpl_2018, data = filtered_data)

print(kruskal_test)

# Random Forest
# Load necessary libraries
library(dplyr)

# Assuming the dataset is already loaded and named 'dataset'
# Split the dataset into three subsets based on the 'LC_simpl_2018' column
# This piece of code has to be repeated for each land use separately
data_grassland <- filter(dataset, LC_simpl_2018 == "Grassland")
data_cropland <- filter(dataset, LC_simpl_2018 == "Cropland")
data_woodland <- filter(dataset, LC_simpl_2018 == "Woodland")
data_grassland$LC_simpl_2018 = NULL
data_cropland$LC_simpl_2018 = NULL
data_woodland$LC_simpl_2018 = NULL

library(rfPermute)
rfPermute_woodland = rfPermute(Multifunctionality ~ .,
                           data = data_woodland, ntree = 1000)

print(rfPermute_woodland$rf)
importance(rfPermute_woodland)
# Do not forget to generate rfPermute_cropland and rfPermute_grassland

#Code for Structural Equations Model (SEM)

library(nlme)
library(lme4)
library(piecewiseSEM)
library(QuantPsyc)

##Read data
SEM_woodlands=read.delim('clipboard',header=T)
SEM_grasslands=read.delim('clipboard',header=T)
SEM_croplands=read.delim('clipboard',header=T)

####Multiple regression based on lm mixed effects model

Function.list <- list(
  # Climate → Soil properties
  lm(MB ~ PRE + TEMP, data = SEM_data),
  lm(DEN ~ PRE + TEMP, data = SEM_data),
  lm(SIL ~ PRE + TEMP, data = SEM_data),
  lm(PH ~ PRE + TEMP, data = SEM_data),
  lm(SOC ~ PRE + TEMP, data = SEM_data),
  lm(K ~ PRE + TEMP, data = SEM_data),
  
  # Soil properties → Microbial groups
  lm(AP ~ PRE + TEMP + MB + DEN + SIL + PH + SOC + K, data = SEM_data),
  lm(CHY ~ PRE + TEMP + MB + DEN + SIL + PH + SOC + K, data = SEM_data),
  lm(ECT ~ PRE + TEMP + MB + DEN + SIL + PH + SOC + K, data = SEM_data),
  lm(MOR ~ PRE + TEMP + MB + DEN + SIL + PH + SOC + K, data = SEM_data),
  lm(SAC ~ PRE + TEMP + MB + DEN + SIL + PH + SOC + K, data = SEM_data),
  lm(UMB ~ PRE + TEMP + MB + DEN + SIL + PH + SOC + K, data = SEM_data),
  lm(YEA ~ PRE + TEMP + MB + DEN + SIL + PH + SOC + K, data = SEM_data),
  lm(AGP7 ~ PRE + TEMP + MB + DEN + SIL + PH + SOC + K, data = SEM_data),
  lm(AGP5 ~ PRE + TEMP + MB + DEN + SIL + PH + SOC + K, data = SEM_data),
  lm(AGP22 ~ PRE + TEMP + MB + DEN + SIL + PH + SOC + K, data = SEM_data),
  lm(ACT ~ PRE + TEMP + MB + DEN + SIL + PH + SOC + K, data = SEM_data),
  lm(CHT ~ PRE + TEMP + MB + DEN + SIL + PH + SOC + K, data = SEM_data),
  lm(DPRO ~ PRE + TEMP + MB + DEN + SIL + PH + SOC + K, data = SEM_data),
  lm(NIT ~ PRE + TEMP + MB + DEN + SIL + PH + SOC + K, data = SEM_data),
  lm(PHO ~ PRE + TEMP + MB + DEN + SIL + PH + SOC + K, data = SEM_data),
  
  # Microbial groups + soil properties → Multifunctionality
  lm(Multifunctionality ~ PRE + TEMP + MB + DEN + SIL + PH + SOC + K +
       AP + CHY + ECT + MOR + SAC + UMB + YEA +
       AGP7 + AGP5 + AGP22 + ACT + CHT + DPRO + NIT + PHO,
     data = SEM_data),
  
  # Correlated errors (example suggestions, adjust as needed)
  PRE %~~% TEMP,
  DEN %~~% MB,
  SIL %~~% PH,
  SOC %~~% K,
  MOR %~~% ECT,
  AP %~~% PHO,
  ACT %~~% CHY,
  YEA %~~% UMB,
  AGP7 %~~% AGP22,
  AGP5 %~~% AGP7
)

#AIC
AIC(Function.psem,AIC.type = "dsep", aicc = FALSE)
