# Survival Analysis with TCGA in R
This repository contains one R script (survival.R) that performs survival analysis on breast cancer patients with high and low expression of SCYL3 gene. 

## Study Data
TCGA-BRCA from GDC portal using TCGAbiolinks R package

Goal: To understand whether high expression of SCYL3 influence the prognosis in breast cancer patients.


## Requirements
- R version 4.2.1
- Package:
•	TCGAbiolinks (v2.25.0): To download the TCGA breast cancer data from TCG portal
•	survminer (v0.4.9): To perform survival analysis and plot survival curves
•	survival (v3.2-13): To perform survival analysis and plot survival curves
•	SummarizedExperiment (v1.24.0): To get gene expression data in a summarized experiment object
•	tidyverse (v1.3.1): To manipulate data
•	DESeq2 (v1.32.0): To obtain DeSeq data set
•	ggplot: For visual representation

##### 3 main things for survival analysis
1. Time-to-event: time till an event happens 
•	In this study, the event is death.
•	So, there is a requirement of the time until death of every patient occurs.
2. Status: censored data; which patients need to be considered in the analysis
•	Need to indicate which patients have to be censored and which will have to be considered for the analysis based on time to death available for the patients.
•	For the patients that have not experienced the event of death will have to be censored out.
3. Event: Strata or group they belong to
•	Divide the cohort into in to two groups: patients with high expression of SCYL3 genes and patients with low expression of the same.  



## Result
1. survfit: The fit object gives us information on the strata that we compared.
- ‘n’ is the number of individuals present in each of the strata.
- Three events were detected in high expressing group and four in low expressing group.
- For low expressing group, the median survival was 2417 days and for high expressing group it could not be computed, suggesting that these patients have not reached the median survival point. This lack of data prevents a direct comparison of survival probabilities between the high and low expressing groups.
- 0.95LCL is the lower bound of the 95% confidence interval.
- 0.95UCL is the upper bound of the 95% confidence interval.

The fit object result is visualized in the Kaplan-Meier curve. 
The risk table below the Kaplan-Meier curve represents the number of patients at risk of an event (death) at specific time points during the study. 
It provides information on the number of patients remaining in each group at different time intervals, indicating how many patients are still being monitored for survival outcomes.
This table helps track the progression of patients in the study and provides insights into the sample size at various time points.
So, after 1000 days 7 patients out of 10 with high expressing SCYL3 gene experienced death and only 3 were being monitored for study and 2 patients out of 10 died and 8 patients were still alive and being monitored. 
After 3000 days, still 1 patient with high expressing SCYL3 gene was alive and being monitored but all the patients with low expressing gene experienced death.

2. survdiff: The fit2 object also gives us information on the strata that we compared.
‘N’ is the number of individuals present in each of the strata.
Observed events: Three events were detected in high expressing group and four in low expressing group.
Expected events (considering that there are no differences in the survival probabilities in both the groups): 2.69 in high expressing group and 4.31 in low expressing group.
The last two columns give us information on the contribution of each group to the chi-square statistics and the variance of the contribution of each group to the chi-square statistics.



## Conclusion
There is significant difference in survival probabilities in patients having high expression of SCYL3 and low expression of SCYL3.
We see a poor prognosis for patients with high SCYL3 expression till 1000 days. 
But after that there is not much difference between the prognosis of the two groups. This could be due to the fact that this study was done only on 20 patients and not the entire cohort.
