# CancerSEEKv1
R code for the algorithm used in Cohen et al. "Detection and localization of surgically resectable cancers with a multi-analyte blood test", Science 2018, 359(6378):926-930.

There are 4 steps to run the full code.

1. Run Analysis_1.sh for submitting the code to the cluster to run it in parallel. It will call for Analysis_1.R that is used to analyze all mutations. A 10-fold cross validation is performed each time Analysis_1.R is called. 

2. Run Analysis_2.R used to calculate the omega value for each mutation and find the maximum omega value for each sample to be used in Step 3. 

Note: you can skip step 1 and 2 and directly use the final results of those steps as found in the file “maxValuesPerSample_20171209_FORJosh.rda” and in the 10 files “allResults_fromLuMethod_*_20171209.rda”, where * ranges from 1 to 10.

3. Run Analysis_3.R to obtain the logistic regression scores needed for the classification calls (cancer vs normal), as well as for sensitivity analysis. It also generates the sensitivities when removing one feature at a time. See file “SummaryResults.txt” for the summary of these analyses.   

4. Run Analysis_4 to predict the tumor tissue type (localization) using only correctly classified cancer samples from Step 3. Random Forest is applied to protein data, gender, and maximum omega score per sample, as obtained in Step 2. The results are saved as “tissueRecog_results_2018-01-13 .xlsx”.
