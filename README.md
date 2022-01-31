# COPD-Mortality-Prediction

This repository shows the list of R code files that were used to undertake the main analysis. This code, however, requires data which needs to be requested 
from the data provider (CPRD). 

Once the data is available, the analysis files can be run sequentially. The recommended sequence to follow is:
A: Load the raw data (analysis_load_data.R)
B: Convert the data into a counting process format (analysis_counting_format.R)
C: Process the covariates (analysis_counting_format.R)
D: Plot the Kaplan-Meir curves (plot_KM_curves)
E: Train and test the model using the processing pipeline developed (analysis_train_test.R)
