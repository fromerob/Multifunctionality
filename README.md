README

This repository contains R code used to analyze the relationship between soil microbial community composition, soil properties, and ecosystem multifunctionality across European soils. The workflow integrates co-occurrence network analysis, multifunctionality metrics, variance partitioning, random forest models, and structural equation modeling. Below is a description of each major section of the code.

Grouping co-occurring OTUs into modules using the WGCNA package

Operational Taxonomic Units (OTUs) were filtered, transformed, and processed using the WGCNA framework to identify modules of co-occurring taxa. Module eigengenes (first principal component of each module) were calculated and used as predictors of ecosystem multifunctionality and individual functions.

Calculating ecosystem multifunctionality indices

Four complementary indices of multifunctionality were computed:
	1.	Multiple-threshold index (number of functions above 25%, 50%, 75% of maximum observed values).
	2.	Z-score index (average of standardized values).
	3.	Normalized average (0–1 min–max scaling).
	4.	Weighted index (enzymes together receive 0.25 total weight, while basal respiration, aggregation, and primary productivity each receive 0.25).

Variance partition analyses

Variance partitioning was conducted to determine the unique and shared contributions of soil properties, microbial community composition, and climate to ecosystem multifunctionality. Results were computed across different groupings (land use, soil texture, pH, and climate region).

Random Forest models

Random forest analyses were performed to assess the relative importance of soil properties, climate, and microbial modules as predictors of multifunctionality. Cross-validated R² and out-of-bag error are reported for each model, ensuring robust model evaluation.

Multi-group Structural Equation Models (SEMs)

Structural equation modeling (lavaan) was applied to explore direct and indirect effects of microbial composition, soil properties, and climate on multifunctionality. Multi-group SEMs were implemented to compare model structure across soil groupings (land use, texture, pH, and climate region).

How to run

1) Software & R packages
	•	R ≥ 4.2 (tested on 4.3.x)
	•	Install required packages:

pkgs <- c(
  "tidyverse","readxl","WGCNA","vegan","microeco","caret","randomForest","rfPermute",
  "lavaan","Matrix","reshape2","RColorBrewer","doParallel","ggplot2","dplyr","tidyr"
)
install.packages(setdiff(pkgs, rownames(installed.packages())))

Note: WGCNA needs allowWGCNAThreads() / enableWGCNAThreads(); avoid parallel conflicts.

2) Data inputs

Place these files in the working directory:
	•	working_dataset_220825.xlsx (sheet: dataset140825) – main dataset for RF.
	•	working_dataset_140825.xlsx / working_dataset_300625.xlsx as referenced in code (update paths if needed).
	•	OTU/ASV tables and taxonomy:
	•	seq_table (OTU abundance; samples as columns, OTUs as rows)
	•	tax_table (taxonomy for seq_table)
	•	TAX_fungi (fungal taxonomy table)
	•	OTU_main_LC, OTU_main_LC_filtered (if used)
	•	A sample metadata table with at least:
	•	SampleID
	•	Functional proxies (6 columns; respiration, 3 enzymes, NPP, aggregation)
	•	Soil properties (e.g., TN, TOC, microbial biomass)
	•	Climate variables
	•	Grouping variables: LC_simpl_2018 (land use), SOIL_TYPE_SIMPL (texture), PH_Class, Climate_zone_simpl

Column order for functions: the code assumes functions are in columns 10:15. If different, adjust indices or rename using the candidate-name detection in the script.

3) Execution order (modules of the script)
	1.	WGCNA modules
	•	Filter/transform OTU table
	•	Build co-occurrence network, detect modules, compute module eigengenes (MEs)
	•	Export module assignments / eigengenes if needed
	2.	Multifunctionality indices
	•	Compute: multiple-threshold (25/50/75%), z-score mean, 0–1 normalized average, weighted index (enzymes share 0.25 total; respiration/NPP/aggregation = 0.25 each)
	3.	Variance partition analyses
	•	Run per grouping (land use, climate region, pH, texture)
	•	Output unique/shared contributions (%)
	4.	Random Forest models
	•	Train RF with rfPermute and caret::train (10-fold CV)
	•	Record cross-validated R² and OOB RMSE
	•	Export variable importance
	5.	Multi-group SEMs (lavaan)
	•	Standardize predictors
	•	Fit unconstrained and constrained models by group (e.g., PH_Class)
	•	Compare via anova()
	•	Export parameter estimates and fit indices

4) Reproducibility notes
	•	Set seeds (already included): set.seed(1) before RF/SEM steps.
	•	Save session info for provenance:

writeLines(capture.output(sessionInfo()), "sessionInfo.txt")

	•	Keep file paths consistent; if you change filenames/sheets, update the corresponding read_excel() calls.

5) Outputs
	•	Tables/CSVs: module assignments, variance partition summaries, RF metrics (CV R², OOB RMSE, importance), SEM parameter estimates and fit measures.
	•	Figures: WGCNA diagnostics, dendrograms/heatmaps, multifunctionality comparisons, RF scatter + correlations, SEM summaries (as produced by plotting code).

6) Common pitfalls
	•	Mismatched SampleIDs: ensure SampleID keys align across OTU tables, taxonomy, and dataset.
	•	Sparse OTUs: WGCNA assumes filtered, non-excessively sparse data (retain OTUs present in ≥10% samples as in script).
	•	Scaling: SEM section scales predictors; avoid re-scaling responses.
	•	Function names/positions: if your function columns aren’t at 10:15, adjust indices or update candidate name lists in the weighted-index block.
