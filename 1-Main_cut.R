library(rstudioapi) 
library(stringr)
library(dplyr)
library(progeny)
library(survival)
library(glmnet)
library(tidyverse)
library(survminer)

# Set Dynamic Working Directory
directory <- getSourceEditorContext()$path
directory <- str_replace(directory, "/[^/]*$", "")
directory <- str_replace(directory, "/[^/]*$", "")
setwd(directory)

# Process mRNA seq Data ##############

# Read data_mrna_seq_v2_rsem.txt
data_mrna_seq <- read.table(paste(directory, "data/lihc_tcga_pan_can_atlas_2018/data_mrna_seq_v2_rsem.txt", sep = "/"), header = TRUE, stringsAsFactors = FALSE, sep = "\t")
# Drop genes that do not have Hugo_Symbol
data_mrna_seq <- data_mrna_seq[data_mrna_seq$Hugo_Symbol != "" & !is.na(data_mrna_seq$Hugo_Symbol), ]
# Replace NA values with 0
data_mrna_seq[is.na(data_mrna_seq)] <- 0
# Save Row names
rownames(data_mrna_seq) <- make.unique(data_mrna_seq$Hugo_Symbol)
# Drop Entrez_Gene_Id and Hugo_Symbol
data_mrna_seq <- subset(data_mrna_seq, select = -c(Entrez_Gene_Id, Hugo_Symbol))
# Rename columns to match samples
colnames(data_mrna_seq) <- gsub("\\.", "-", colnames(data_mrna_seq))
# Transpose the data_mrna_seq data to match the patient/sample data
data_mrna_seq <- as.data.frame(t(data_mrna_seq))

# Process Patient & Sample Data ################################################

# Read data_clinical_patient.txt and ata_clinical_sample.txt
data_clinical_patient <- read.table(paste(directory, "data/lihc_tcga_pan_can_atlas_2018/data_clinical_patient.txt", sep = "/"), header = TRUE, stringsAsFactors = FALSE, sep = "\t")
data_clinical_sample <- read.csv(paste(directory, "data/lihc_tcga_pan_can_atlas_2018/data_clinical_sample.txt", sep = "/"), header = TRUE, stringsAsFactors = FALSE, sep = "\t")
# As we are reading the file via .csv because read.table does not work for some reason, we need to rename coulmns and drop additional rows
colnames(data_clinical_sample) <- data_clinical_sample[4, ]
data_clinical_sample <- data_clinical_sample[-c(1:4), ]
# Select the columns we might need
#data_clinical_patient <- subset(data_clinical_patient, select = c(PATIENT_ID, AGE, SEX, SMOKING_PACK_YEARS, STAGE, OS_MONTHS, OS_STATUS))
data_clinical_patient <- subset(data_clinical_patient, select = c(PATIENT_ID, OS_MONTHS, OS_STATUS))
data_clinical_sample <- subset(data_clinical_sample, select = c(PATIENT_ID, SAMPLE_ID, CANCER_TYPE_DETAILED, TMB_NONSYNONYMOUS))
# Merge Two data frmaes
data_clinical <- merge(data_clinical_sample, data_clinical_patient, by = "PATIENT_ID", all.x = TRUE)
# Drop patients that do not have OS_STATUS or OS_MONTHS
data_clinical <- data_clinical[data_clinical$OS_STATUS != "" & !is.na(data_clinical$OS_STATUS), ]
data_clinical <- data_clinical[data_clinical$OS_MONTHS != "" & !is.na(data_clinical$OS_MONTHS), ]
data_clinical <- data_clinical[((data_clinical$OS_MONTHS > 1 & data_clinical$OS_STATUS == "0:LIVING") | data_clinical$OS_STATUS == "1:DECEASED"), ]
# Surv() function in the {survival} package accepts by default TRUE/FALSE, where TRUE is event and FALSE is censored; 1/0 where 1 is event and 0 is censored
data_clinical$SURV_STATUS <- ifelse(data_clinical$OS_STATUS == "1:DECEASED", 1, 0)

# Select Pateints who have mRNA seq data
data_clinical <- data_clinical[data_clinical$SAMPLE_ID %in% rownames(data_mrna_seq), ]
data_mrna_seq <- data_mrna_seq[rownames(data_mrna_seq) %in% data_clinical$SAMPLE_ID,]

# Process Data Mutations File ##################################################

# Read Data Mutations file
data_mutations <- read.table(paste(directory, "data/lihc_tcga_pan_can_atlas_2018/data_mutations.txt", sep = "/"), header = TRUE, stringsAsFactors = FALSE, sep = "\t")
# Subset the columns we might need
data_mutations <- subset(data_mutations, select = c(Hugo_Symbol, Consequence, Tumor_Sample_Barcode))
# Merge mutations with their Consequence
data_mutations$Hugo_Symbol <- paste(data_mutations$Hugo_Symbol, data_mutations$Consequence, sep = "_")
# Drop Consequence column
data_mutations <- subset(data_mutations, select = -c(Consequence))
# Drop Samples from mutations file that are not in clinical data
data_mutations <- data_mutations[data_mutations$Tumor_Sample_Barcode %in% data_clinical$SAMPLE_ID, ]
# Spread the data set to create a matrix, where columns will be samples, rows mutations and value would indicate if a mutation is present in sample
data_mutations <- data_mutations %>%
  distinct(Hugo_Symbol, Tumor_Sample_Barcode) %>%
  mutate(value = 1) %>%
  spread(Tumor_Sample_Barcode, value, fill = 0)
# Name the row names
rownames(data_mutations) <- data_mutations$Hugo_Symbol
# Drop the Gene names from column
data_mutations <- subset(data_mutations, select = -c(Hugo_Symbol))
# Transpose the data_mutations data to match the patient/sample data
data_mutations <- as.data.frame(t(data_mutations))

# Process mRNA for PROGENY #####################################################

# Create a data frame for progeny
progeny_data_mrna_seq <- data.frame("Gene" = colnames(data_mrna_seq))
# Add expression values to it
progeny_data_mrna_seq <- cbind(progeny_data_mrna_seq, t(data_mrna_seq))
# Run Progeny and save the output result
progeny_data <- as.data.frame(progeny(as.matrix(progeny_data_mrna_seq[,-1]), scale=TRUE, organism="Human", top = 100, perm = 1))

# Process Methylation data #########################################################

# Read Data Methylation file
data_methylation <- read.table(paste(directory, "data/lihc_tcga_pan_can_atlas_2018/data_methylation_hm27_hm450_merged.txt", sep = "/"), header = TRUE, stringsAsFactors = FALSE, sep = "\t")
# Drop rows with NA in the "NAME" column
data_methylation <- data_methylation[complete.cases(data_methylation$NAME), ]
rownames(data_methylation) <- make.unique(data_methylation$NAME)
data_methylation <- subset(data_methylation, select = -c(ENTITY_STABLE_ID, NAME, DESCRIPTION, TRANSCRIPT_ID))
# Replace NA values with 0
data_methylation[is.na(data_methylation)] <- 0
# Rename columns to match samples
colnames(data_methylation) <- gsub("\\.", "-", colnames(data_methylation))
data_methylation <- as.data.frame(t(data_methylation))

################################################################################
# Feature extraction ###########################################################
################################################################################

# We also need to supply the survival object to the glmnet function
data_clinical_surv <- Surv(time = data_clinical$OS_MONTHS, event = data_clinical$SURV_STATUS)

# RNA Seq Feature extraction ###################################################

# Run glmnet()
fit_rna <- glmnet(x = data_mrna_seq, data_clinical_surv, family = "cox")
# Select the genes to explain as much variation as possible
cfs_rna <- as.data.frame(as.matrix((coef(fit_rna, s = 0.001789)))) # 84%
# Create Gene column
cfs_rna$Genes <- rownames(cfs_rna)
# Rename rownames and column names
rownames(cfs_rna) <- NULL
colnames(cfs_rna) <- c("Coefficient", "Gene")
# Replace . in Gene column with -
cfs_rna$Gene <- gsub("\\.", "-", cfs_rna$Gene)
# Filter Genes by Coefficient
cfs_rna <- cfs_rna %>% filter(Coefficient > 0.1 | Coefficient < -0.1)
# Create Column with absolute values of Genes Coefficients
cfs_rna$Abs_Coefficient <- abs(cfs_rna$Coefficient)
# Reorder genes based on absolute value
cfs_rna <- cfs_rna[order(-cfs_rna$Abs_Coefficient),]
# Print formula for Cox regression
print(paste(cfs_rna$Gene[which(cfs_rna$Abs_Coefficient > 1)], collapse = " + "))

# Run Cox Regression (Drop OR9K2 gene because it has infinite coefficient)
fit.coxph.mrna <- coxph(data_clinical_surv ~ OR13C8 + OR10J1 + OR2M4 + OR5AP2 + PPBPP1 + OR11H4 + OR2J2 + HPVC1 + `SNAR-G1` + TAS2R39 + C4orf11 + OR11L1 + LINC00917 + ACSM4 + C4orf35 + C1orf68 + THEG + PRSS37 + SNORA71D + `PHACTR2-AS1` + SPANXN2 + CSN1S1 + PRDM14 + SPDYE4 + OR6F1 + TAS2R42 + GJD4 + CIB4 + PRSS42 + CXorf50B + MORF4 + `KRTAP3-2` + RTP2,
                        data = data_mrna_seq)
ggforest(fit.coxph.mrna, data = data_mrna_seq)
# Drop OR13C8, OR13C8, OR10J1, `SNAR-G1`, TAS2R39, OR2M4, PPBPP1, C4orf11, OR11L1, LINC00917 Genes due to too high Hazard Ratio values
fit.coxph.mrna <- coxph(data_clinical_surv ~ OR5AP2 + OR11H4 + OR2J2 + HPVC1 + ACSM4 + C4orf35 + C1orf68 + THEG + PRSS37 + SNORA71D + `PHACTR2-AS1` + SPANXN2 + CSN1S1 + PRDM14 + SPDYE4 + OR6F1 + TAS2R42 + GJD4 + CIB4 + PRSS42 + CXorf50B + MORF4 + `KRTAP3-2` + RTP2,
                        data = data_mrna_seq)
ggforest(fit.coxph.mrna, data = data_mrna_seq)

# Cox Regression on genes from literature: TERT, CTNNB1, AXIN1, APC, P53, RB1, CCNA2, CCNE1, PTEN, ARID1A, ARID2, RPS6KA3, NFE2L2, CCND1, FGF19, VEGFA, MYC, MET, HFE, SERPINA1, G6PC and SLC37A4, HMBS and UROD, FAH, ATP7B, GPC3, SPP1, GOLM1, EGF, PDGFB, IGF1, and MTOR 
fit.coxph.mrna.literature <- coxph(data_clinical_surv ~ TERT + CTNNB1 + AXIN1 + APC + TP53 + RB1 + CCNA2 + CCNE1 + PTEN + ARID1A + ARID2 + RPS6KA3 + NFE2L2 + CCND1 + FGF19 + VEGFA + MYC + MET + HFE + SERPINA1 + G6PC + SLC37A4 + HMBS + UROD + FAH + ATP7B + GPC3 + SPP1 + GOLM1 + EGF + PDGFB + IGF1 + MTOR,
                        data = data_mrna_seq)
ggforest(fit.coxph.mrna.literature, data = data_mrna_seq)
# Drop MTOR, PDGFB, ATP7B, RB1, NFE2L2, MYC, CCND1, RPS6KA3, FGF19, PTEN, UROD, FAH, AXIN1, TP53, CCNE1, SERPINA1, SLC37A4, CCNA2
fit.coxph.mrna.literature <- coxph(data_clinical_surv ~ TERT + CTNNB1 + APC + ARID1A + ARID2 + VEGFA + MET + HFE + G6PC + HMBS + GPC3 + SPP1 + GOLM1 + EGF + IGF1,
                        data = data_mrna_seq)
ggforest(fit.coxph.mrna.literature, data = data_mrna_seq)

# We need to find a best way to separate the data frame
best_split_function <- function(values, survival) {
  # Define base threshold
  threshold <- min(values)
  # Find min and max values of mRNA
  min_value <- min(values)
  max_value <- max(values)
  # Create a data frame for mRNA and Survival + Output Dataframe for storing thresholds
  values_df <- data.frame("Gene" = values, "Surv" = survival)
  threshold_df <- data.frame(Threshold = numeric(),AUC_Difference = numeric(), stringsAsFactors = FALSE)
  # For each threshold with 0.01 step:
  while(threshold<max(values)){
    # Put here to skip the min separation
    threshold = threshold+0.01
    # Copy mRNA values and find there min and max values
    gg_values_df <- values_df
    min_category <- paste0(min_value,"-",threshold)
    max_category <- paste0(threshold,"-",max_value)
    # Replace mRNA values with two levels
    gg_values_df$Gene <- ifelse(values<threshold, min_category, max_category)
    # Add survival data
    gg_values_df$Surv <- survival
    # Check if in the separation we have two groups
    if(length(unique(gg_values_df$Gene))>1){
      # Get labels of two groups
      group_1 <- unique(gg_values_df$Gene)[1]
      group_2 <- unique(gg_values_df$Gene)[2]
      # Separate two groups by labels
      gg_values_1 <- gg_values_df[gg_values_df$Gene==group_1,]
      gg_values_2 <- gg_values_df[gg_values_df$Gene==group_2,]
      # Check if both groups have at least 2 members
      if(nrow(gg_values_1)>1 & nrow(gg_values_2)>1)
      {
        # Compute a survival curve for both groups
        fit_1 <- survfit(gg_values_1$Surv~1, data = gg_values_1)
        fit_2 <- survfit(gg_values_2$Surv~1, data = gg_values_2)
        # Obtain data for Integration
        surv_prob_1 <- fit_1$surv
        time_points_1 <- fit_1$time
        surv_prob_2 <- fit_2$surv
        time_points_2 <- fit_2$time
        # Perform Integration via Trapezoidal Rule Formula to obtain Area Under the Curve for two graphs
        auc_1 <- sum(diff(time_points_1) * (surv_prob_1[-1] + surv_prob_1[-length(surv_prob_1)]) / 2)
        auc_2 <- sum(diff(time_points_2) * (surv_prob_2[-1] + surv_prob_2[-length(surv_prob_2)]) / 2)
        # Find difference in Area under the Curve of the two graphs
        auc_dif <- abs(auc_1-auc_2)
        # Save data to repeat the process 
        threshold_df <- rbind(threshold_df, data.frame(Threshold = threshold, AUC_Difference = auc_dif)) 
      }
    }
  }
  # Return the coefficient granting AUC difference closest to mean AUC difference
  return(threshold_df$Threshold[which.min(abs(threshold_df$AUC_Difference - mean(threshold_df$AUC_Difference)))])
  # Return the coefficient granting max AUC
  #return(threshold_df$Threshold[which.max(abs(threshold_df$AUC_Difference))])
}

# Copy mRNA data
gg_data_mrna <- data_mrna_seq
# Cut the values of important genes for Hazard Ratio into two groups
gg_data_mrna$OR5AP2 <- cut(gg_data_mrna$OR5AP2, breaks = c(min(gg_data_mrna$OR5AP2), best_split_function(data_mrna_seq$OR5AP2, data_clinical_surv), max(gg_data_mrna$OR5AP2)), include.lowest = TRUE)
gg_data_mrna$OR11H4 <- cut(gg_data_mrna$OR11H4, breaks = c(min(gg_data_mrna$OR11H4), best_split_function(data_mrna_seq$OR11H4, data_clinical_surv), max(gg_data_mrna$OR11H4)), include.lowest = TRUE)
gg_data_mrna$ACSM4 <- cut(gg_data_mrna$ACSM4, breaks = c(min(gg_data_mrna$ACSM4), best_split_function(data_mrna_seq$ACSM4, data_clinical_surv), max(gg_data_mrna$ACSM4)), include.lowest = TRUE)
gg_data_mrna$C4orf35 <- cut(gg_data_mrna$C4orf35, breaks = c(min(gg_data_mrna$C4orf35), best_split_function(data_mrna_seq$C4orf35, data_clinical_surv), max(gg_data_mrna$C4orf35)), include.lowest = TRUE)
gg_data_mrna$C1orf68 <- cut(gg_data_mrna$C1orf68, breaks = c(min(gg_data_mrna$C1orf68), best_split_function(data_mrna_seq$C1orf68, data_clinical_surv), max(gg_data_mrna$C1orf68)), include.lowest = TRUE)
gg_data_mrna$PRSS37 <- cut(gg_data_mrna$PRSS37, breaks = c(min(gg_data_mrna$PRSS37), best_split_function(data_mrna_seq$PRSS37, data_clinical_surv), max(gg_data_mrna$PRSS37)), include.lowest = TRUE)
gg_data_mrna$PRDM14 <- cut(gg_data_mrna$PRDM14, breaks = c(min(gg_data_mrna$PRDM14), best_split_function(data_mrna_seq$PRDM14, data_clinical_surv), max(gg_data_mrna$PRDM14)), include.lowest = TRUE)
gg_data_mrna$SPDYE4 <- cut(gg_data_mrna$SPDYE4, breaks = c(min(gg_data_mrna$SPDYE4), best_split_function(data_mrna_seq$SPDYE4, data_clinical_surv), max(gg_data_mrna$SPDYE4)), include.lowest = TRUE)
gg_data_mrna$CIB4 <- cut(gg_data_mrna$CIB4, breaks = c(min(gg_data_mrna$CIB4), best_split_function(data_mrna_seq$CIB4, data_clinical_surv), max(gg_data_mrna$CIB4)), include.lowest = TRUE)

# Build survival curves for each gene
fit_OR5AP2 <- survfit(data_clinical_surv ~ OR5AP2, data = gg_data_mrna)
fit_OR11H4 <- survfit(data_clinical_surv ~ OR11H4, data = gg_data_mrna)
fit_ACSM4 <- survfit(data_clinical_surv ~ ACSM4, data = gg_data_mrna)
fit_C4orf35 <- survfit(data_clinical_surv ~ C4orf35, data = gg_data_mrna)
fit_C1orf68 <- survfit(data_clinical_surv ~ C1orf68, data = gg_data_mrna)
fit_PRSS37 <- survfit(data_clinical_surv ~ PRSS37, data = gg_data_mrna)
fit_PRDM14 <- survfit(data_clinical_surv ~ PRDM14, data = gg_data_mrna)
fit_SPDYE4 <- survfit(data_clinical_surv ~ SPDYE4, data = gg_data_mrna)
fit_CIB4 <- survfit(data_clinical_surv ~ CIB4, data = gg_data_mrna)

# Draw plots for each gene
ggsurvplot(fit_OR5AP2, data = gg_data_mrna, pval = T) # https://www.frontiersin.org/articles/10.3389/fonc.2018.00033/full / https://www.frontiersin.org/articles/10.3389/fphys.2020.574082/full
ggsurvplot(fit_OR11H4, data = gg_data_mrna, pval = T) 
ggsurvplot(fit_ACSM4, data = gg_data_mrna, pval = T) #  https://pubmed.ncbi.nlm.nih.gov/33340617/
ggsurvplot(fit_C4orf35, data = gg_data_mrna, pval = T)
ggsurvplot(fit_C1orf68, data = gg_data_mrna, pval = T) # PVAL>0.05!
ggsurvplot(fit_PRSS37, data = gg_data_mrna, pval = T)
ggsurvplot(fit_PRDM14, data = gg_data_mrna, pval = T)
ggsurvplot(fit_SPDYE4, data = gg_data_mrna, pval = T)
ggsurvplot(fit_CIB4, data = gg_data_mrna, pval = T)

# Progeny Feature Extraction ###################################################

# Run glmnet()
fit_progeny <- glmnet(x = progeny_data, data_clinical_surv, family = "cox")
# Select genes
cfs_progeny <- as.data.frame(as.matrix((coef(fit_progeny, s = 0.001482)))) # 3.10% 
cfs_progeny$Pathway <- rownames(cfs_progeny)
rownames(cfs_progeny) <- NULL
colnames(cfs_progeny) <- c("Coefficient", "Pathway")
cfs_progeny <- cfs_progeny %>% filter(Coefficient > 0.1 | Coefficient < -0.1)
# Create Column with absolute values of Genes Coefficients
cfs_progeny$Abs_Coefficient <- abs(cfs_progeny$Coefficient)
# Reorder genes based on absolute value
cfs_progeny <- cfs_progeny[order(-cfs_progeny$Abs_Coefficient),]
# Print formula for Cox regression
print(paste(cfs_progeny$Pathway[which(cfs_progeny$Abs_Coefficient > 0.1)], collapse = " + "))

fit.coxph.progeny <- coxph(data_clinical_surv ~ Trail + TNFa + p53 + Hypoxia + Androgen + `JAK-STAT` + Estrogen + VEGF, 
                           data = progeny_data)
ggforest(fit.coxph.progeny, data = progeny_data)

# Hypoxia Genes: HIF1A, HIF3A, PDGFB, IGF1, INS, NRG1, SMAD3, SMAD4, ETS1 https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6908932/
# p53 Genes: MDM4, MDM2, STEAP3, PTEN, DDB2, TP53, CDKN1A, SFN, RPRM, GADD45G, GTSE1, ZNF385A, FAS, LRDD, BAX, SIVA1, EI24, SHISA5, AIFM2, IGFBP3, SERPINE1
"SERPINE1" %in% colnames(data_mrna_seq)

# Methylation Feature Extraction ###############################################

# We also need to supply the survival object to the glmnet function
data_clinical_methylation <- data_clinical[data_clinical$SAMPLE_ID %in% rownames(data_methylation),]
data_methylation <- data_methylation[rownames(data_methylation) %in% data_clinical$SAMPLE_ID,]
data_clinical_methylation_surv <- Surv(time = data_clinical_methylation$OS_MONTHS, event = data_clinical_methylation$SURV_STATUS)

# Run glmnet()
fit_methylation <- glmnet(x = data_methylation, data_clinical_methylation_surv, family = "cox")
# Select the genes to explain as much variation as possible
cfs_methylation <- as.data.frame(as.matrix((coef(fit_methylation, s = 0.001822))))
# Create Gene column
cfs_methylation$Gene <- rownames(cfs_methylation)
# Rename rownames and column names
rownames(cfs_methylation) <- NULL
colnames(cfs_methylation) <- c("Coefficient", "Gene")
# Filter Genes by Coefficient
cfs_methylation <- cfs_methylation %>% filter(Coefficient > 0.1 | Coefficient < -0.1)
# Create Column with absolute values of Genes Coefficients
cfs_methylation$Abs_Coefficient <- abs(cfs_methylation$Coefficient)
# Reorder genes based on absolute value
cfs_methylation <- cfs_methylation[order(-cfs_methylation$Abs_Coefficient),]
# Print formula for Cox regression
print(paste(cfs_methylation$Gene[which(cfs_methylation$Abs_Coefficient > 0.5 & cfs_methylation$Abs_Coefficient < 2.5)], collapse = " + "))
print(paste(cfs_methylation$Gene[which(cfs_methylation$Abs_Coefficient > 1)], collapse = " + "))

# Run Cox Regression - methylation data
fit.coxph.methylation <- coxph(data_clinical_methylation_surv ~ PVRL1 + `MAP4K2;MEN1` + C12orf48 + CALM1.1 + TCEA1.1 + C11orf65.1 + CDC23 + CHMP2B + RHOG + ZMAT2.1 + CH25H + WAPAL.1 + `OTOP2;USH1G.1` + ALG13.1 + CCDC19 + `C16orf42;GNPTG` + TLN2 + TRIOBP + APC2 + DNAH8 + CD44.1 + GRB10.1 + `PARD6A;ACD.2` + GATA4.3 + GJA4 + PRKD3 + KCNG1 + BBS5.1 + LOX + MDFI.1 + `SLC12A4;LCAT.1` + CCDC141 + SCN4B + ZNF273 + ZMPSTE24.1 + TMEM22 + `YIF1B;C19orf33` + PDE6B + EIF2C2 + ZNF33B + NQO1 + MS4A3 + LRRC20 + COMMD2 + CLTC + CRABP2 + GMIP.1 + RDH11.1 + ZRANB1 + ZFP2, 
                               data = data_methylation)
ggforest(fit.coxph.methylation, data = data_methylation)
# Drop PVRL1, CALM1.1, COMMD2, EIF2C2, ZNF273, TCEA1.1, ZMPSTE24.1, CDC23, CHMP2B, C16orf42;GNPTG, PARD6A;ACD.2, DNAH8, PRKD3, GMIP.1, ALG13.1, TLN2, ZMAT2.1, CD44.1, LOX, ZRANB1, CRABP2, CLTC, LRRC20, ZNF33B, `YIF1B;C19orf33`
fit.coxph.methylation <- coxph(data_clinical_methylation_surv ~ `MAP4K2;MEN1` + C12orf48 + C11orf65.1 + RHOG + CH25H + WAPAL.1 + `OTOP2;USH1G.1` + CCDC19 + TRIOBP + APC2 + GRB10.1 + GATA4.3 + GJA4 + KCNG1 + BBS5.1 + MDFI.1 + `SLC12A4;LCAT.1` + CCDC141 + SCN4B + TMEM22 + PDE6B + NQO1 + MS4A3 + RDH11.1 + ZFP2, 
                               data = data_methylation)
ggforest(fit.coxph.methylation, data = data_methylation)

# Copy Methylation data
gg_data_methylation <- data_methylation
# Cut the values of important genes for Hazard Ratio into two groups
gg_data_methylation$MAP4K2_MEN1 <- cut(gg_data_methylation$`MAP4K2;MEN1`, breaks = c(min(gg_data_methylation$`MAP4K2;MEN1`), best_split_function(data_methylation$`MAP4K2;MEN1`, data_clinical_methylation_surv), max(gg_data_methylation$`MAP4K2;MEN1`)), include.lowest = TRUE)
gg_data_methylation$C12orf48 <- cut(gg_data_methylation$C12orf48, breaks = c(min(gg_data_methylation$C12orf48), best_split_function(data_methylation$C12orf48, data_clinical_methylation_surv), max(gg_data_methylation$C12orf48)), include.lowest = TRUE)
gg_data_methylation$C11orf65.1 <- cut(gg_data_methylation$C11orf65.1, breaks = c(min(gg_data_methylation$C11orf65.1), best_split_function(data_methylation$C11orf65.1, data_clinical_methylation_surv), max(gg_data_methylation$C11orf65.1)), include.lowest = TRUE)
gg_data_methylation$WAPAL.1 <- cut(gg_data_methylation$WAPAL.1, breaks = c(min(gg_data_methylation$WAPAL.1), best_split_function(data_methylation$WAPAL.1, data_clinical_methylation_surv), max(gg_data_methylation$WAPAL.1)), include.lowest = TRUE)
gg_data_methylation$OTOP2_USH1G.1 <- cut(gg_data_methylation$`OTOP2;USH1G.1`, breaks = c(min(gg_data_methylation$`OTOP2;USH1G.1`), best_split_function(data_methylation$`OTOP2;USH1G.1`, data_clinical_methylation_surv), max(gg_data_methylation$`OTOP2;USH1G.1`)), include.lowest = TRUE)
gg_data_methylation$APC2 <- cut(gg_data_methylation$APC2, breaks = c(min(gg_data_methylation$APC2), best_split_function(data_methylation$APC2, data_clinical_methylation_surv), max(gg_data_methylation$APC2)), include.lowest = TRUE)
gg_data_methylation$GRB10.1 <- cut(gg_data_methylation$GRB10.1, breaks = c(min(gg_data_methylation$GRB10.1), best_split_function(data_methylation$GRB10.1, data_clinical_methylation_surv), max(gg_data_methylation$GRB10.1)), include.lowest = TRUE)
gg_data_methylation$GATA4.3 <- cut(gg_data_methylation$GATA4.3, breaks = c(min(gg_data_methylation$GATA4.3), best_split_function(data_methylation$GATA4.3, data_clinical_methylation_surv), max(gg_data_methylation$GATA4.3)), include.lowest = TRUE)
gg_data_methylation$GJA4 <- cut(gg_data_methylation$GJA4, breaks = c(min(gg_data_methylation$GJA4), best_split_function(data_methylation$GJA4, data_clinical_methylation_surv), max(gg_data_methylation$GJA4)), include.lowest = TRUE)
gg_data_methylation$KCNG1 <- cut(gg_data_methylation$KCNG1, breaks = c(min(gg_data_methylation$KCNG1), best_split_function(data_methylation$KCNG1, data_clinical_methylation_surv), max(gg_data_methylation$KCNG1)), include.lowest = TRUE)
gg_data_methylation$BBS5.1 <- cut(gg_data_methylation$BBS5.1, breaks = c(min(gg_data_methylation$BBS5.1), best_split_function(data_methylation$BBS5.1, data_clinical_methylation_surv), max(gg_data_methylation$BBS5.1)), include.lowest = TRUE)
gg_data_methylation$SLC12A4_LCAT.1 <- cut(gg_data_methylation$`SLC12A4;LCAT.1`, breaks = c(min(gg_data_methylation$`SLC12A4;LCAT.1`), best_split_function(data_methylation$`SLC12A4;LCAT.1`, data_clinical_methylation_surv), max(gg_data_methylation$`SLC12A4;LCAT.1`)), include.lowest = TRUE)
gg_data_methylation$SCN4B <- cut(gg_data_methylation$SCN4B, breaks = c(min(gg_data_methylation$SCN4B), best_split_function(data_methylation$SCN4B, data_clinical_methylation_surv), max(gg_data_methylation$SCN4B)), include.lowest = TRUE)
gg_data_methylation$TMEM22 <- cut(gg_data_methylation$TMEM22, breaks = c(min(gg_data_methylation$TMEM22), best_split_function(data_methylation$TMEM22, data_clinical_methylation_surv), max(gg_data_methylation$TMEM22)), include.lowest = TRUE)
gg_data_methylation$PDE6B <- cut(gg_data_methylation$PDE6B, breaks = c(min(gg_data_methylation$PDE6B), best_split_function(data_methylation$PDE6B, data_clinical_methylation_surv), max(gg_data_methylation$PDE6B)), include.lowest = TRUE)
gg_data_methylation$RDH11.1 <- cut(gg_data_methylation$RDH11.1, breaks = c(min(gg_data_methylation$RDH11.1), best_split_function(data_methylation$RDH11.1, data_clinical_methylation_surv), max(gg_data_methylation$RDH11.1)), include.lowest = TRUE)
gg_data_methylation$ZFP2 <- cut(gg_data_methylation$ZFP2, breaks = c(min(gg_data_methylation$ZFP2), best_split_function(data_methylation$ZFP2, data_clinical_methylation_surv), max(gg_data_methylation$ZFP2)), include.lowest = TRUE)

# Build survival curves for each gene
fit_MAP4K2_MEN1 <- survfit(data_clinical_methylation_surv ~ MAP4K2_MEN1, data = gg_data_methylation)
fit_C12orf48 <- survfit(data_clinical_methylation_surv ~ C12orf48, data = gg_data_methylation)
fit_C11orf65.1 <- survfit(data_clinical_methylation_surv ~ C11orf65.1, data = gg_data_methylation)
fit_WAPAL.1 <- survfit(data_clinical_methylation_surv ~ WAPAL.1, data = gg_data_methylation)
fit_OTOP2_USH1G.1 <- survfit(data_clinical_methylation_surv ~ OTOP2_USH1G.1, data = gg_data_methylation)
fit_APC2 <- survfit(data_clinical_methylation_surv ~ APC2, data = gg_data_methylation)
fit_GRB10.1 <- survfit(data_clinical_methylation_surv ~ GRB10.1, data = gg_data_methylation)
fit_GATA4.3 <- survfit(data_clinical_methylation_surv ~ GATA4.3, data = gg_data_methylation)
fit_GJA4 <- survfit(data_clinical_methylation_surv ~ GJA4, data = gg_data_methylation)
fit_KCNG1 <- survfit(data_clinical_methylation_surv ~ KCNG1, data = gg_data_methylation)
fit_BBS5.1 <- survfit(data_clinical_methylation_surv ~ BBS5.1, data = gg_data_methylation)
fit_SLC12A4_LCAT.1 <- survfit(data_clinical_methylation_surv ~ SLC12A4_LCAT.1, data = gg_data_methylation)
fit_SCN4B <- survfit(data_clinical_methylation_surv ~ SCN4B, data = gg_data_methylation)
fit_TMEM22 <- survfit(data_clinical_methylation_surv ~ TMEM22, data = gg_data_methylation)
fit_PDE6B <- survfit(data_clinical_methylation_surv ~ PDE6B, data = gg_data_methylation)
fit_RDH11.1 <- survfit(data_clinical_methylation_surv ~ RDH11.1, data = gg_data_methylation)
fit_ZFP2 <- survfit(data_clinical_methylation_surv ~ ZFP2, data = gg_data_methylation)

# Draw plots for each gene
ggsurvplot(fit_MAP4K2_MEN1, data = gg_data_methylation, pval = T)
ggsurvplot(fit_C12orf48, data = gg_data_methylation, pval = T)
ggsurvplot(fit_C11orf65.1, data = gg_data_methylation, pval = T)
ggsurvplot(fit_WAPAL.1, data = gg_data_methylation, pval = T)
ggsurvplot(fit_OTOP2_USH1G.1, data = gg_data_methylation, pval = T)
ggsurvplot(fit_APC2, data = gg_data_methylation, pval = T)
ggsurvplot(fit_GRB10.1, data = gg_data_methylation, pval = T)
ggsurvplot(fit_GATA4.3, data = gg_data_methylation, pval = T)
ggsurvplot(fit_GJA4, data = gg_data_methylation, pval = T)
ggsurvplot(fit_KCNG1, data = gg_data_methylation, pval = T)
ggsurvplot(fit_BBS5.1, data = gg_data_methylation, pval = T)
ggsurvplot(fit_SLC12A4_LCAT.1, data = gg_data_methylation, pval = T)
ggsurvplot(fit_SCN4B, data = gg_data_methylation, pval = T)
ggsurvplot(fit_TMEM22, data = gg_data_methylation, pval = T)
ggsurvplot(fit_PDE6B, data = gg_data_methylation, pval = T)
ggsurvplot(fit_RDH11.1, data = gg_data_methylation, pval = T)
ggsurvplot(fit_ZFP2, data = gg_data_methylation, pval = T)

# Cox Regression on MAP4K2 and MEN1 genes based on mrna data
fit.coxph.mrna.methylation <- coxph(data_clinical_surv ~ MAP4K2 + MEN1 + C12orf48 + C11orf65 + GATA4 + BBS5 + SCN4B + TMEM22 + PDE6B + ZFP2,
                                        data = data_mrna_seq)
ggforest(fit.coxph.mrna.methylation, data = data_mrna_seq) # Nothing

# Run Cox regression for literature genes based on methylation data
print(grep("MTOR", colnames(data_methylation), value = TRUE))
fit.coxph.methylation.lit <- coxph(data_clinical_methylation_surv ~ CTNNB1 + AXIN1 + APC + APC.1 + APC.2 + RB1 + RB1.1 + RB1.2 + RB1.3 + RB1.4 + RB1.5 + RB1.6 + RB1.7 + RB1.8 + RB1.9 +
                                 CCNA2 + CCNA2.1 + `PTEN;KILLIN` + `PTEN;KILLIN.1` + `PTEN;KILLIN.2` + ARID1A + RPS6KA3 + NFE2L2 + CCND1 + CCND1.1 + CCND1.2 + CCND1.3 + CCND1.4 + CCND1.5 + CCND1.6 +
                                 FGF19 + VEGFA + MET + MET.1 + SERPINA1 + SERPINA1.1 + SLC37A4 + HMBS + `UROD;HECTD3` + `ALG11;ATP7B` + GOLM1 + PDGFB + IGF1 + `ANGPTL7;MTOR` + MTOR + `ANGPTL7;MTOR.1`, 
                               data = data_methylation)
# Drop Genes: NFE2L2, CTNNB1, APC, PTEN;KILLIN, PTEN;KILLIN.1, PTEN;KILLIN.2`, CCNA2, CCNA2.1, CCND1, CCND1.2, CCND1.3, UROD;HECTD3, CCND1.6, HMBS, ANGPTL7;MTOR.1, VEGFA, CCND1.5, RB1.4, MET, IGF1, RB1.6
fit.coxph.methylation.lit <- coxph(data_clinical_methylation_surv ~ AXIN1 + APC.1 + APC.2 + RB1 + RB1.1 + RB1.2 + RB1.3 + RB1.5 + RB1.7 + RB1.8 + RB1.9 + ARID1A + RPS6KA3 + CCND1.1 + CCND1.4 + FGF19 + MET.1 + SERPINA1 + SERPINA1.1 + SLC37A4 + `ALG11;ATP7B` + GOLM1 + PDGFB + `ANGPTL7;MTOR` + MTOR, 
                                   data = data_methylation)
ggforest(fit.coxph.methylation.lit, data = data_methylation)

# Copy Methylation data
gg_data_methylation <- data_methylation
# Cut the values of important genes for Hazard Ratio into two groups
gg_data_methylation$RB1 <- cut(gg_data_methylation$RB1, breaks = c(min(gg_data_methylation$RB1), best_split_function(data_methylation$RB1, data_clinical_methylation_surv), max(gg_data_methylation$RB1)), include.lowest = TRUE)
gg_data_methylation$APC.2 <- cut(gg_data_methylation$APC.2, breaks = c(min(gg_data_methylation$APC.2), best_split_function(data_methylation$APC.2, data_clinical_methylation_surv), max(gg_data_methylation$APC.2)), include.lowest = TRUE)
gg_data_methylation$RB1.2 <- cut(gg_data_methylation$RB1.2, breaks = c(min(gg_data_methylation$RB1.2), best_split_function(data_methylation$RB1.2, data_clinical_methylation_surv), max(gg_data_methylation$RB1.2)), include.lowest = TRUE)
gg_data_methylation$ALG11_ATP7B <- cut(gg_data_methylation$`ALG11;ATP7B`, breaks = c(min(gg_data_methylation$`ALG11;ATP7B`), best_split_function(data_methylation$`ALG11;ATP7B`, data_clinical_methylation_surv), max(gg_data_methylation$`ALG11;ATP7B`)), include.lowest = TRUE)
gg_data_methylation$MTOR <- cut(gg_data_methylation$MTOR, breaks = c(min(gg_data_methylation$MTOR), best_split_function(data_methylation$MTOR, data_clinical_methylation_surv), max(gg_data_methylation$MTOR)), include.lowest = TRUE)

# Build survival curves for each gene
fit_RB1 <- survfit(data_clinical_methylation_surv ~ RB1, data = gg_data_methylation)
fit_APC.2 <- survfit(data_clinical_methylation_surv ~ APC.2, data = gg_data_methylation)
fit_RB1.2 <- survfit(data_clinical_methylation_surv ~ RB1.2, data = gg_data_methylation)
fit_ALG11_ATP7B <- survfit(data_clinical_methylation_surv ~ ALG11_ATP7B, data = gg_data_methylation)
fit_MTOR <- survfit(data_clinical_methylation_surv ~ MTOR, data = gg_data_methylation)

# Draw plots for each gene
ggsurvplot(fit_RB1, data = gg_data_methylation, pval = T)
ggsurvplot(fit_APC.2, data = gg_data_methylation, pval = T)
ggsurvplot(fit_RB1.2, data = gg_data_methylation, pval = T)
ggsurvplot(fit_ALG11_ATP7B, data = gg_data_methylation, pval = T)
ggsurvplot(fit_MTOR, data = gg_data_methylation, pval = T)

# Run Cox Regression based methylation data for mRNA genes
print(grep("CIB4", colnames(data_methylation), value = TRUE))
fit.coxph.methylation <- coxph(data_clinical_methylation_surv ~ C4orf35 + C4orf35.1, 
                               data = data_methylation)
ggforest(fit.coxph.methylation, data = data_methylation)

# Mutations Feature Extraction #################################################

# We also need to supply the survival object to the glmnet function
data_clinical_mutations <- data_clinical[data_clinical$SAMPLE_ID %in% rownames(data_mutations), ]
data_clinical_mutations_surv <- Surv(time = data_clinical_mutations$OS_MONTHS, event = data_clinical_mutations$SURV_STATUS)

# Run glmnet()
fit_mut <- glmnet(x = data_mutations, data_clinical_mutations_surv, family = "cox")
# Select the genes to explain as much variation as possible
cfs_mutations = as.data.frame(as.matrix((coef(fit_mut, s = 0.001062)))) # 89.51%
# Create Mutations column
cfs_mutations$Genes <- rownames(cfs_mutations)
# Rename rownames and column names
rownames(cfs_mutations) <- NULL
colnames(cfs_mutations) <- c("Coefficient", "Mutation")
# Filter Mutations by Coefficient
cfs_mutations <- cfs_mutations %>% filter(Coefficient > 0.1 | Coefficient < -0.1)
cfs_mutations$Abs_Coefficient <- abs(cfs_mutations$Coefficient)
cfs_mutations <- cfs_mutations[order(-cfs_mutations$Abs_Coefficient),]
# Print formula for Cox regression
print(paste(cfs_mutations$Mutation[which(cfs_mutations$Abs_Coefficient > 1)], collapse = " + "))

# Cox Regression on mutations data
fit.coxph.mutations <- coxph(data_clinical_mutations_surv ~ KIF16B_synonymous_variant + C6orf223_synonymous_variant + PRKACB_5_prime_UTR_variant + DPYD_missense_variant + HORMAD2_missense_variant + MUC6_missense_variant + LCE6A_synonymous_variant + UGT2A2_missense_variant + PLCXD2_5_prime_UTR_variant + ENGASE_missense_variant + OR5AL1_downstream_gene_variant + JAK3_missense_variant + MUSTN1_synonymous_variant + STK39_synonymous_variant + HTT_missense_variant + F10_synonymous_variant + BID_3_prime_UTR_variant + ABCG1_stop_gained + HPCAL1_3_prime_UTR_variant + MYH9_intron_variant + ABCC8_synonymous_variant + CNOT4_missense_variant + GPR179_missense_variant,
                              data = data_mutations)
ggforest(fit.coxph.mutations, data = data_mutations)

# Survival curves for mutation genes
fit_PRKACB_5_prime_UTR_variant <- survfit(data_clinical_mutations_surv ~ PRKACB_5_prime_UTR_variant, data = data_mutations)
fit_DPYD_missense_variant <- survfit(data_clinical_mutations_surv ~ DPYD_missense_variant, data = data_mutations)
fit_MUC6_missense_variant <- survfit(data_clinical_mutations_surv ~ MUC6_missense_variant, data = data_mutations)
fit_OR5AL1_downstream_gene_variant <- survfit(data_clinical_mutations_surv ~ OR5AL1_downstream_gene_variant, data = data_mutations)
fit_STK39_synonymous_variant <- survfit(data_clinical_mutations_surv ~ STK39_synonymous_variant, data = data_mutations)
fit_F10_synonymous_variant <- survfit(data_clinical_mutations_surv ~ F10_synonymous_variant, data = data_mutations)
fit_BID_3_prime_UTR_variant <- survfit(data_clinical_mutations_surv ~ BID_3_prime_UTR_variant, data = data_mutations)
fit_CNOT4_missense_variant <- survfit(data_clinical_mutations_surv ~ CNOT4_missense_variant, data = data_mutations)
# Draw plots for mutation genes
ggsurvplot(fit_PRKACB_5_prime_UTR_variant, data = data_mutations, pval = T)
ggsurvplot(fit_DPYD_missense_variant, data = data_mutations, pval = T)
ggsurvplot(fit_MUC6_missense_variant, data = data_mutations, pval = T)
ggsurvplot(fit_OR5AL1_downstream_gene_variant, data = data_mutations, pval = T)
ggsurvplot(fit_STK39_synonymous_variant, data = data_mutations, pval = T)
ggsurvplot(fit_F10_synonymous_variant, data = data_mutations, pval = T)
ggsurvplot(fit_BID_3_prime_UTR_variant, data = data_mutations, pval = T)
ggsurvplot(fit_CNOT4_missense_variant, data = data_mutations, pval = T)

# Test discoreved genes from mutations on mRNA data
print(grep("OR5AL1", colnames(data_mrna_seq), value = TRUE))
fit.coxph.mrna.mutations <- coxph(data_clinical_surv ~ PRKACB + MUC6 + STK39 + F10 + BID + CNOT4, 
                        data = data_mrna_seq)
ggforest(fit.coxph.mrna.mutations, data = data_mrna_seq)
# Test discoreved genes from mutations on methylation data
print(grep("CNOT4", colnames(data_methylation), value = TRUE))
fit.coxph.methylation.mutations <- coxph(data_clinical_methylation_surv ~ PRKACB + PRKACB.1 + STK39 + CNOT4 + CNOT4.1, 
                                               data = data_methylation)
ggforest(fit.coxph.methylation.mutations, data = data_methylation)

# mRNA confirmed genes
print(grep("SPDYE4", colnames(data_mutations), value = TRUE))
fit.coxph.mutations.mrna <- coxph(data_clinical_mutations_surv ~ OR11H4_missense_variant + PRDM14_synonymous_variant, 
                                  data = data_mutations)
ggforest(fit.coxph.mutations.mrna, data = data_mutations)
# Build survival curves for OR11H4_missense_variant
fit_OR11H4_missense_variant <- survfit(data_clinical_mutations_surv ~ OR11H4_missense_variant, data = data_mutations)
# Draw plots for OR11H4_missense_variant
ggsurvplot(fit_OR11H4_missense_variant, data = data_mutations, pval = T)

# Cox Regression on mutations data - literature genes
print(grep("TP53_splice_region_variant", colnames(data_mutations), value = TRUE))
fit.coxph.mutations.lit <- coxph(data_clinical_mutations_surv ~ CTNNB1_missense_variant + AXIN1_frameshift_variant + AXIN1_missense_variant + APC_missense_variant + RB1_frameshift_variant + RPS6KA3_missense_variant + NFE2L2_missense_variant + TP53_frameshift_variant + TP53_inframe_deletion + TP53_missense_variant + TP53_stop_gained + ARID1A_frameshift_variant + ARID1A_missense_variant + ARID1A_stop_gained + ARID2_frameshift_variant + ARID2_stop_gained + MET_missense_variant + G6PC_missense_variant + SLC37A4_missense_variant + MTOR_missense_variant,
                             data = data_mutations)
ggforest(fit.coxph.mutations.lit, data = data_mutations)

# Survival curves for mutation genes (literature)
fit_TP53_inframe_deletion <- survfit(data_clinical_mutations_surv ~ TP53_inframe_deletion, data = data_mutations)
fit_TP53_stop_gained <- survfit(data_clinical_mutations_surv ~ TP53_stop_gained, data = data_mutations)
fit_ARID2_stop_gained <- survfit(data_clinical_mutations_surv ~ ARID2_stop_gained, data = data_mutations)
# Draw plots for mutation genes literature
ggsurvplot(fit_TP53_inframe_deletion, data = data_mutations, pval = T)
ggsurvplot(fit_TP53_stop_gained, data = data_mutations, pval = T)
ggsurvplot(fit_ARID2_stop_gained, data = data_mutations, pval = T)

# Cox Regression of TP53 gene based on mrna data
fit.coxph.mrna.p53 <- coxph(data_clinical_surv ~ TP53,
                                        data = data_mrna_seq)
ggforest(fit.coxph.mrna.p53, data = data_mrna_seq) # Nothing
# Cox Regression of TP53 gene based on methylation data
print(grep("TP53", colnames(data_methylation), value = TRUE))
fit.coxph.progeny.methylation.hypoxia <- coxph(data_clinical_methylation_surv ~ `TP53;WRAP53`, 
                                               data = data_methylation)
ggforest(fit.coxph.progeny.methylation.hypoxia, data = data_methylation)

# Copy Methylation data
gg_data_methylation <- data_methylation
# Cut the values of TP53 gene
gg_data_methylation$TP53_WRAP53 <- cut(gg_data_methylation$`TP53;WRAP53`, breaks = c(min(gg_data_methylation$`TP53;WRAP53`), best_split_function(data_methylation$`TP53;WRAP53`, data_clinical_methylation_surv), max(gg_data_methylation$`TP53;WRAP53`)), include.lowest = TRUE)
# Build survival curves for TP53
fit_TP53_WRAP53 <- survfit(data_clinical_methylation_surv ~ TP53_WRAP53, data = gg_data_methylation)
# Draw plot
ggsurvplot(fit_TP53_WRAP53, data = gg_data_methylation, pval = T)

################################################################################

# Additional Validation of Progeny Genes

# Cox Regression on Hypoxia genes on mrna data
fit.coxph.mrna.progeny.hypoxia <- coxph(data_clinical_surv ~ HIF1A + HIF3A + PDGFB + IGF1 + INS + NRG1 + SMAD3 + SMAD4 + ETS1,
                                        data = data_mrna_seq)
ggforest(fit.coxph.mrna.progeny.hypoxia, data = data_mrna_seq) # Nothing

# Cox Regression on p53 genes on mrna data
fit.coxph.mrna.progeny.p53 <- coxph(data_clinical_surv ~ MDM4 + MDM2 + STEAP3 + PTEN + DDB2 + TP53 + CDKN1A + SFN + RPRM + GADD45G + GTSE1 + ZNF385A + FAS + LRDD + BAX + SIVA1 + EI24 + SHISA5 + AIFM2 + IGFBP3 + SERPINE1,
                                    data = data_mrna_seq)
ggforest(fit.coxph.mrna.progeny.p53, data = data_mrna_seq) # Nothing

# Cox Regression on Hypoxia genes on methylation data
print(grep("SHISA5", colnames(data_methylation), value = TRUE))
fit.coxph.progeny.methylation.hypoxia <- coxph(data_clinical_methylation_surv ~ HIF1A + PDGFB + IGF1 + SMAD4 + ETS1 + `IGF2;INS-IGF2` + `IGF2;INS-IGF2.1` + `INS-IGF2;INS` + `IGF2;IGF2AS;INS-IGF2` + `IGF2;IGF2AS;INS-IGF2.1` + `INS-IGF2;INS.1` + `IGF2;IGF2AS;INS-IGF2.2` + `IGF2;IGF2AS;INS-IGF2.3` + `IGF2;INS-IGF2.2` + `IGF2;IGF2AS;INS-IGF2.4`, 
                               data = data_methylation)
ggforest(fit.coxph.progeny.methylation.hypoxia, data = data_methylation) # Not important for Hazard Ratio
# Drop Genes: IGF2;INS-IGF2, HIF1A, SMAD4, IGF2;IGF2AS;INS-IGF2.4, ETS1, IGF2;IGF2AS;INS-IGF2
fit.coxph.progeny.methylation.hypoxia <- coxph(data_clinical_methylation_surv ~ + PDGFB + IGF1 + `IGF2;INS-IGF2.1` + `INS-IGF2;INS` + `IGF2;IGF2AS;INS-IGF2.1` + `INS-IGF2;INS.1` + `IGF2;IGF2AS;INS-IGF2.2` + `IGF2;IGF2AS;INS-IGF2.3` + `IGF2;INS-IGF2.2`, 
                                               data = data_methylation)
ggforest(fit.coxph.progeny.methylation.hypoxia, data = data_methylation) # Not important for Hazard Ratio

# Cox Regression on p53 genes on methylation data
fit.coxph.progeny.methylation.p53 <- coxph(data_clinical_methylation_surv ~  MDM4 + MDM2 + STEAP3 + `PTEN;KILLIN` + `PTEN;KILLIN.1` + `PTEN;KILLIN.2` + DDB2 + CDKN1A + SFN + RPRM + GADD45G + `GTSE1;CN5H6.4` + `ACTA2;FAS` + LRDD + SIVA1 + EI24 + AIFM2 + IGFBP3 + SERPINE1,
                                           data = data_methylation)
ggforest(fit.coxph.progeny.methylation.p53, data = data_methylation) # Not important for Hazard Ratio
# Drop STEAP3, PTEN;KILLIN, PTEN;KILLIN.1, PTEN;KILLIN.2, SERPINE1, EI24, GADD45G, MDM2, RPRM
fit.coxph.progeny.methylation.p53 <- coxph(data_clinical_methylation_surv ~ MDM4 + DDB2 + CDKN1A + SFN + `GTSE1;CN5H6.4` + `ACTA2;FAS` + LRDD + SIVA1 + AIFM2 + IGFBP3, 
                               data = data_methylation)
ggforest(fit.coxph.progeny.methylation.p53, data = data_methylation) # Not important for Hazard Ratio

# Copy Methylation data
gg_data_methylation <- data_methylation
# Cut the values of important genes for Hazard Ratio into two groups
gg_data_methylation$INS_IGF2.2 <- cut(gg_data_methylation$`IGF2;IGF2AS;INS-IGF2.2`, breaks = c(min(gg_data_methylation$`IGF2;IGF2AS;INS-IGF2.2`), best_split_function(data_methylation$`IGF2;IGF2AS;INS-IGF2.2`, data_clinical_methylation_surv), max(gg_data_methylation$`IGF2;IGF2AS;INS-IGF2.2`)), include.lowest = TRUE)
gg_data_methylation$ACTA2_FAS <- cut(gg_data_methylation$`ACTA2;FAS`, breaks = c(min(gg_data_methylation$`ACTA2;FAS`), best_split_function(data_methylation$`ACTA2;FAS`, data_clinical_methylation_surv), max(gg_data_methylation$`ACTA2;FAS`)), include.lowest = TRUE)
# Build survival curves for each gene
fit_INS_IGF2.2 <- survfit(data_clinical_methylation_surv ~ INS_IGF2.2, data = gg_data_methylation)
fit_ACTA2_FAS <- survfit(data_clinical_methylation_surv ~ ACTA2_FAS, data = gg_data_methylation)
# Draw plots for each gene
ggsurvplot(fit_INS_IGF2.2, data = gg_data_methylation, pval = T)
ggsurvplot(fit_ACTA2_FAS, data = gg_data_methylation, pval = T)

# Cox Regression for INS_IGF2.2 based on mRNA data
fit.coxph.mrna.ins.igf2 <- coxph(data_clinical_surv ~ INS + IGF2,
                        data = data_mrna_seq)
ggforest(fit.coxph.mrna.ins.igf2, data = data_mrna_seq)

# Cox Regression for INS_IGF2.2 based on mutations data
print(grep("IGF2", colnames(data_mutations), value = TRUE))
# Cannot be fit
#fit.coxph.mutations.ins.igf2 <- coxph(data_clinical_mutations_surv ~ INS_missense_variant + IGF2_frameshift_variant,
#                                 data = data_mutations)
#ggforest(fit.coxph.mutations.ins.igf2, data = data_mutations)

# Cox Regression on Hypoxia genes on mutations data
print(grep("SERPINE1", colnames(data_mutations), value = TRUE))
fit.coxph.progeny.mutations.hypoxia <- coxph(data_clinical_mutations_surv ~ HIF1A_missense_variant + HIF3A_missense_variant + IGF1R_missense_variant + IGF1R_synonymous_variant + NRG1_intron_variant + NRG1_missense_variant + NRG1_synonymous_variant + SMAD4_missense_variant, 
                                  data = data_mutations)
ggforest(fit.coxph.progeny.mutations.hypoxia, data = data_mutations)

# Cox Regression on p53 genes on mutations data
fit.coxph.progeny.mutations.p53 <- coxph(data_clinical_mutations_surv ~ STEAP3_frameshift_variant + STEAP3_synonymous_variant + PTEN_5_prime_UTR_variant + PTEN_splice_donor_variant + TP53_frameshift_variant + TP53_inframe_deletion + TP53_intron_variant + TP53_missense_variant + TP53_splice_acceptor_variant + TP53_stop_gained + CDKN1A_3_prime_UTR_variant + CDKN1A_start_lost + EI24_missense_variant, 
                                            data = data_mutations)
ggforest(fit.coxph.progeny.mutations.p53, data = data_mutations)

# Build survival curves for progeny pathway genes based on mutation data
fit_HIF1A_missense_variant <- survfit(data_clinical_mutations_surv ~ HIF1A_missense_variant, data = data_mutations)
fit_PTEN_splice_donor_variant <- survfit(data_clinical_mutations_surv ~ PTEN_splice_donor_variant, data = data_mutations)
fit_TP53_inframe_deletion <- survfit(data_clinical_mutations_surv ~ TP53_inframe_deletion, data = data_mutations)
fit_TP53_stop_gained <- survfit(data_clinical_mutations_surv ~ TP53_stop_gained, data = data_mutations)
fit_EI24_missense_variant <- survfit(data_clinical_mutations_surv ~ EI24_missense_variant, data = data_mutations)
# Draw plots
ggsurvplot(fit_HIF1A_missense_variant, data = data_mutations, pval = T)
ggsurvplot(fit_PTEN_splice_donor_variant, data = data_mutations, pval = T)
ggsurvplot(fit_TP53_inframe_deletion, data = data_mutations, pval = T)
ggsurvplot(fit_TP53_stop_gained, data = data_mutations, pval = T)
ggsurvplot(fit_EI24_missense_variant, data = data_mutations, pval = T)

#Additional Validation of Methylation data

# Cox Regression on MAP4K2 and MEN1 genes based on mutations data
print(grep("ZFP2", colnames(data_mutations), value = TRUE))
fit.coxph.mutations.methylation <- coxph(data_clinical_mutations_surv ~ MAP4K2_missense_variant + `BBS5_missense_variant,splice_region_variant` + BBS5_splice_acceptor_variant + PDE6B_missense_variant + ZFP2_missense_variant, 
                                         data = data_mutations)
ggforest(fit.coxph.mutations.methylation, data = data_mutations)

# Build survival curves for MAP4K2 and MEN genes based on mutation data
print(grep("MEN1", colnames(data_mutations), value = TRUE))
fit_MAP4K2_missense_variant <- survfit(data_clinical_mutations_surv ~ MAP4K2_missense_variant, data = data_mutations)
fit_MAP4K2_missense_variant <- survfit(data_clinical_mutations_surv ~ MAP4K2_missense_variant, data = data_mutations)
fit_PTEN_splice_donor_variant <- survfit(data_clinical_mutations_surv ~ PTEN_splice_donor_variant, data = data_mutations)
# Draw plots
ggsurvplot(fit_HIF1A_missense_variant, data = data_mutations, pval = T)
ggsurvplot(fit_PTEN_splice_donor_variant, data = data_mutations, pval = T)