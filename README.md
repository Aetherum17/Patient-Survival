### HCC Survivor Biomarkers

# Background

Liver cancer is one of the most common types of cancer, the number of which is constantly rising. It is 
expected that the number of cases of this disease will exceed the threshold of one million patients per year 
by the year 2025.[1] Herewith, about 90% of cases of liver cancer belong to hepatocellular carcinoma (HCC), 
a malignancy derived from hepatocytes, from which about 781 thousand people died in 2018 worldwide.[1,3]
The contribution to the development of HCC is made by both environmental factors, including excessive 
alcohol consumption, infection with hepatitis B and C viruses, as well as the development of liver cirrhosis 
associated with them, and genetic factors, among which mutations in the promoter of the TERT gene 
encoding one of the components of telomerase, mutations in the CTNNB1, AXIN1, and APC genes leading 
to activation of the Wnt/β-catenin signalling pathway, as well as mutations in the P53, RB1, CCNA2, CCNE1, 
PTEN, ARID1A, ARID2, RPS6KA3, and NFE2L2 genes that control the cell cycle can be distinguished. Other 
genes important for the occurrence of HCC are CCND1, FGF19, VEGFA, MYC, and MET, the overexpression 
of which will lead to the activation of various oncogenic signalling pathways.[1] In addition, a number of 
diseases have been identified, the presence of which also increases the chances of developing HCC. Among 
them are hemochromatosis (HFE gene), alpha 1-antitrypsin deficiency (SERPINA1 gene), glycogen storage 
diseases (G6PC and SLC37A4 genes), porphyria (HMBS and UROD genes), tyrosinemia (FAH gene) and
Wilson's disease (ATP7B gene).[3] Moreover, various genes involved in epigenetic regulation, oxidative 
stress, and AKT-mTOR and MAPK signalling pathways are also involved in the development of HCC.[1]
The life expectancy of patients with HCC depends on many factors, but the main one is the stage of cancer 
at the time of diagnosis. In advanced stages, life expectancy is usually a few months, but with early 
detection of the disease, a patient's life can be extended up to five years. Thus, early detection of HCC is 
one of the key factors in increasing the survival of people with this disease. Therefore, an active search for 
biomarkers for diagnosing HCC is currently underway, during which GPC3, OPN, GP73 proteins and VEGF, 
EGF, PDGF, IGF, and mTOR genes have already been identified, the increased level of 
biosynthesis/expression of which is associated with the development of this type of cancer.[2]
However, there is nearly no information in the articles on biomarkers associated with life expectancy in 
patients who have already been diagnosed with HCC. Therefore, as the main goal of this project,
the search for biomarkers potentially associated with the overall survival of people with HCC was chosen.

# Materials and Methods
 
To accomplish the setup goal, the dataset "Liver Hepatocellular Carcinoma" by TCGA, PanCancer Atlas 
(2018) downloaded from the cBioPortal website was used.[4] This dataset was processed using a script 
written in R 5 programming language in RStudio IDE6 using rstudioapi[7], string[8], dplyr[9], progeny[10], 
survival[11], glmnet[12], tidyverse[13] and survminer[14] libraries. This involved extracting patient and sample data, 
processing mRNA Expression data, mutation data from whole exome sequencing and methylation data, 
performing Cox Proportional Hazards Regression and Progeny Analysis and building a prediction model of 
patient overall survival with features discovered during the performed analysis. 
Extracting patient and sample data involved selecting columns from their tables corresponding to the Patient 
ID, Sample ID and Overall Survival time in months. Alive patients who were censored before 1st month of 
observations were removed from the analysis so as not to bias it. Additionally, some patients/samples were 
removed while performing Cox Proportional Hazards Regression analysis with mRNA, Methylation, and 
Mutations data if their material was not used for those examinations.
Cox Proportional Hazards Regression analysis itself was performed on processed data, which involved its 
transformation in a form suitable for fitting a generalized linear model of the Cox family via the glmnet() 
function. The most important features of the data under study (mRNA expression normalized 
values/methylation normalized values/presence of mutations) were selected based on their coefficient 
value from the generalized linear model (GLM) with the highest lambda value to explain the most variation 
in the data. The exact thresholds of coefficient value were selected empirically in each test separately in 
order to exclude extremely high or low Hazard ratios, though the absolute coefficient values were usually 
between 0.1 and 2.5. However, as this selection procedure does not always remove all extreme Hazard 
Ratio values, some features were dropped manually if their Hazard Ratios were more than 100 or less than 
0.01, as those features might be present in too few numbers of samples.
The remaining data features were selected for plotting survival curves to prove their importance by having 
Hazard Ratio different from 1 and a p-value less than 0.05. For continuous variables, such as mRNA 
expression or Methylation, the split threshold was determined by looping through each possible threshold 
with the step of 0.1 and calculating the difference in the Area Under the Curve (AUC) sized of two graphs. 
The threshold with the closest AUC difference to the mean AUC difference of all thresholds was selected 
for plotting two survival curves. Additionally, the p-value for survival curves was calculated, and only the 
features with a p-value less than 0.05 were selected as important for the patient’s survival. AUC area for 
survival curves was calculated via Integration via the Trapezoidal Rule Formula.
