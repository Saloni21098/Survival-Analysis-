# 29-03-2024
# survival analysis using TCGA data
# setwd("C:/Users/salon/OneDrive/Desktop/Survival_Analysis_TCGA")
# =============================================================================================

# Load the libraries
install.packages('fastmap')
install.packages('BiocManager')
library(BiocManager)
BiocManager::install('Biostrings')
library(TCGAbiolinks) # to download the TCGA breast cancer data from TCG portal
install.packages('survminer')
library(survminer)  # to perform survival analysis and plot survival curves
library(survival)  # to perform survival analysis and plot survival curves
library(SummarizedExperiment) # to get gene expression data in a summarized experiment object
install.packages('tidyverse')
library(tidyverse) # to manipulate the data
library(DESeq2)
library(ggplot2)
# ==============================================================================================

# Get clinical data for TCGA-BRCA cohort 
clinical_data_brca <- GDCquery_clinic('TCGA-BRCA')

any(colnames(clinical_data_brca) %in% c('vital_status', 'days_to_last_follow_up', 'days_to_death'))
which(colnames(clinical_data_brca) %in% c('vital_status', 'days_to_last_follow_up', 'days_to_death'))

# Looking at some variables associated with survival
clinical_data_brca[,c(9,39,45)]
table(clinical_data_brca$vital_status)

# days_to_death, that is the number of days passed from the initial diagnosis to the patient’s death (clearly, this is only relevant for dead patients)
# days_to_last_follow_up that is the number of days passed from the initial diagnosis to the last visit.

# Change the vital_status values to get the status (Alive = False, Dead = True)
clinical_data_brca$status <- ifelse(clinical_data_brca$vital_status == 'Alive', FALSE, TRUE)

# Create an 'overall survival' variable that is equal to days_to_death for dead patients, 
# and to days_to_last_follow_up for patients who are still alive
clinical_data_brca$overall_survival <- ifelse(clinical_data_brca$vital_status == 'Alive',
                                              clinical_data_brca$days_to_last_follow_up,
                                              clinical_data_brca$days_to_death)
# =====================================================================================================================================================

# Get gene expression data for entire cohort
# https://bioconductor.org/packages/release/bioc/vignettes/TCGAbiolinks/inst/doc/query.html
gdc_projects <- getGDCprojects()
gdc_project_summary <- getProjectSummary('TCGA-BRCA')

query_TCGA_brca <- GDCquery(project = 'TCGA-BRCA', data.category = 'Transcriptome Profiling',
                            data.type = 'Gene Expression Quantification',
                            workflow.type = 'STAR - Counts', access = 'open',
                            experimental.strategy = 'RNA-Seq', sample.type = 'Primary Tumor')
output_brca <- getResults(query_TCGA_brca)

# Get 20 primary tissue sample barcodes
tumor <- output_brca$cases[1:20]

# Get gene expression data of 20 primary tumors
query_TCGA_brca <- GDCquery(project = 'TCGA-BRCA', data.category = 'Transcriptome Profiling',
                            data.type = 'Gene Expression Quantification',
                            workflow.type = 'STAR - Counts', access = 'open',
                            experimental.strategy = 'RNA-Seq', 
                            sample.type = 'Primary Tumor', barcode = tumor)
# ============================================================================================

# Download data
GDCdownload(query_TCGA_brca)
# ============================================================================================

# Prepare data
# Get counts
tcga_brca_data <- GDCprepare(query_TCGA_brca, summarizedExperiment = TRUE)
brca_matrix <- assay(tcga_brca_data, 'unstranded')
brca_matrix[1:10, 1:10]

# Extract gene and the sample metadata from summarizedExperiment object
gene_metadata <- as.data.frame(rowData(tcga_brca_data))
col_data <- as.data.frame(colData(tcga_brca_data))
# ============================================================================================

# Variance stabilization transformation (vst) counts to be used in survival analysis
# Set up countData object
dds <- DESeqDataSetFromMatrix(countData = brca_matrix, colData = col_data, design = ~ 1)

# Removing genes with sum total of 10 reads across all samples
dds10 <- dds[rowSums(counts(dds) >= 10),]

# vst
vsd <- vst(dds10, blind = F) # returns another summarizedExperiment
brca_matrix_vst <- assay(vsd)
brca_matrix_vst[1:10, 1:10]

# Get data for SCYL3 gene and add gene_metadata information to it
brca_scyl3 <- brca_matrix_vst %>% 
  as.data.frame() %>% 
  rownames_to_column(var = 'gene_id') %>%
  gather(key = 'case_id', value = 'counts', -gene_id) %>% 
  left_join(., gene_metadata, by = 'gene_id') %>% 
  filter(gene_name == 'SCYL3') 

# Get median value
median_value <- median(brca_scyl3$counts)

# Denote which cases have higher or lower expression of SCYL3 than median count
brca_scyl3$strata <- ifelse(brca_scyl3$counts >= median_value, "HIGH", 'LOW')

# Add clinical information to brca_rad52
brca_scyl3$case_id <- gsub("-01.*", '', brca_scyl3$case_id)
brca_scyl3 <- merge(brca_scyl3, clinical_data_brca, by.x = 'case_id', by.y = 'submitter_id')
# ===========================================================================================

# Fitting survival curve (for visual analysis)
fit <- survfit(Surv(overall_survival, status) ~ strata, data = brca_scyl3)
fit

ggsurvplot(fit, data = brca_scyl3, pval = T, risk.table = T)+
  labs(x = "Time (in days)")

# Fitting survival curve (for programmatic analysis)
fit2 <- survdiff(Surv(overall_survival, status) ~ strata, data = brca_rad52)
fit2
# If the observed and the expected value are different, then it can be concluded that there is a difference in survival pattern of
# two groups of cancer patients having high and low SCYL3 gene expression.
  

  
  
  
  
