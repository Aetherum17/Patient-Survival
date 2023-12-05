# HCC Survivor Biomarkers

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
of which will lead to the activation of various oncogenic signalling pathways.[1] 

In addition, a number of diseases have been identified, the presence of which also increases the chances of developing HCC. Among 
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

# Results

First of all, Cox Proportional Hazards Regression Analysis (CPHRA) was performed on genes from mRNA 
data that, in total, explain 84% variation of the survival data and have an absolute GLM coefficient greater 
than 1. After dropping OR13C8, OR13C8, OR10J1, SNAR-G1, TAS2R39, OR2M4, PPBPP1, C4orf11, OR11L1 
and LINC00917 genes, which had extremely high Hazard Ratio values, a plot of Cox Proportional Hazards 
Model, showed on Image 1 was obtained. Next OR5AP2, OR11H4, ACSM4, C4orf35, C1orf68, PRSS37, 
PRDM14, SPDYE4 and CIB4 genes from this plot with a p-value less than 0.05 were selected for plotting 
survival curves, as shown in Image 2. The negative effect on the patient’s overall survival in cases of 
increased expression was confirmed for OR5AP2, OR11H4, ACSM4, C4orf35, PRSS37, PRDM14, SPDYE4 and 
CIB4 genes.

Additionally, CPHRA was performed on TERT, CTNNB1, AXIN1, APC, P53, RB1, CCNA2, CCNE1, PTEN, 
ARID1A, ARID2, RPS6KA3, NFE2L2, CCND1, FGF19, VEGFA, MYC, MET, HFE, SERPINA1, G6PC, SLC37A4, 
HMBS, UROD, FAH, ATP7B, GPC3, SPP1, GOLM1, EGF, PDGFB, IGF1, and MTOR genes are mentioned in the 
literature in regard to HCC.[1,2,3] Its plot is shown in Image 3 with MTOR, PDGFB, ATP7B, RB1, NFE2L2, MYC, 
CCND1, RPS6KA3, FGF19, PTEN, UROD, FAH, AXIN1, TP53, CCNE1, SERPINA1, SLC37A4, and CCNA2 genes are 
dropped to increase the clarity of the graph as their p-values were more than 0.05. However, all their 
Hazard ratios based on mRNA data were equal to 1, meaning that they are not important for the patient's 
survival.

Though as OR5AP2, OR11H4, ACSM4, C4orf35, PRDM14 and CIB4 genes were still defined to be important 
for patient survival, their Hazard Ratio was also checked based on methylation and mutation data. The 
results for methylation are presented in Image 4, featuring only the C4orf35 gene, as methylation of others 
was not measured, and as Hazard Ratios of C4orf35 are close to 1 with a p-value greater than 0.05, we can 
conclude that methylation level of this gene is not important for patients survival. However, the results of 
CPHRA on mutations data, despite containing only information about OR11H4 and PRDM14 gene, showed 
that a Missense Mutation of OR11H4 can be associated with a negative effect on survival, as it has a 
Hazard Ratio of 13.8 with a p-value of 0.01 (Img. 5).

![image](https://github.com/Aetherum17/Patient-Survival/assets/46795020/a5019d4b-f9f7-4852-aee2-eb710215d993)

![image](https://github.com/Aetherum17/Patient-Survival/assets/46795020/d2b357d9-2751-451d-b967-7a8a38940588)

![image](https://github.com/Aetherum17/Patient-Survival/assets/46795020/7e4c4b72-d9bd-48cb-b02b-4a919a3dfee0)

![image](https://github.com/Aetherum17/Patient-Survival/assets/46795020/9a8fd8b7-6c34-4548-9e62-4e5c7d24998c)

![image](https://github.com/Aetherum17/Patient-Survival/assets/46795020/9527516a-dba6-4adf-86df-872feda10620)

Afterwards, CPHRA was performed on the results of Progeny Analysis, highlighting two pathways -
“Hypoxia” and “p53”, activity of which might be related to the decreased overall survival of patients with 
HCC as the Hazard ratio of the “Hypoxia” pathway is 1.3 with a p-value lower than 0.05 and a Hazard ratio of 
“p53” pathway is 1.23 with a p-value also lower than 0.05 (Img. 6). As the genes HIF1A, HIF3A, PDGFB, IGF1, 
INS, NRG1, SMAD3, SMAD4, ETS1 play a major role in “Hypoxia” pathway15 and genes MDM4, MDM2, 
STEAP3, PTEN, DDB2, TP53, CDKN1A, SFN, RPRM, GADD45G, GTSE1, ZNF385A, FAS, LRDD, BAX, SIVA1, EI24, 
SHISA5, AIFM2, IGFBP3, SERPINE1 play a major role in “p53” pathway16 connection between their 
expression, methylation and mutations with patient survival were also checked. 

Results of CPHRA for those genes based on mRNA data are presented in image 7, showing no association 
with HCC as for all genes, either Hazard Ratio is equal to 1 or the p-value is greater than 0.05. 
Results of CPHRA for progeny genes based on methylation data (Img. 8) were more promising. In the 
“Hypoxia” pathway, INS-IGF2.2 genes look to be associated with survival, as their Hazard Ratio was 0.21 
with a p-value of 0.004. Further inspection of those genes using survival plots confirmed the association, as 
its p-value was less than 0.05 (Img. 9a). In the “p53” pathway, peculiar results were obtained for ACTA2;FAS 
genes, as their Hazard ratio was 3.67 with a p-value of 0.015. However, its p-value of 0.21 at the survival 
plot also does not allow us to confirm its association with HCC disease (Img. 9a). Moreover, analysis of INS
and IGF2 genes did not show an association with patient survival based on mRNA data.
Similar to methylation data, interesting results were obtained by performing CPHRA based on mutation 
data for the “Hypoxia” pathway, as there was a Missense Mutation of the HIF1A gene with a Hazard Ratio of 11.9 
(p-value<0.001), and for the “p53” pathway, as there four mutations with p-values less than 0.05 were 
discovered: PTEN (splice donor variant), TP53 (Inframe Deletion), TP53 (Stop Gained), EI24 (Missense 
Variant) (Img. 10). As can be seen in Image 11, showing survival plots of those genes, the presence of all those 
mutations can potentially be associated with HCC, as their p-values are less than 0.05.

![image](https://github.com/Aetherum17/Patient-Survival/assets/46795020/f623c3c8-770d-499e-9a9c-a950ab2ea232)

![image](https://github.com/Aetherum17/Patient-Survival/assets/46795020/f5b45b94-b98b-4531-9cf1-6128e1961249)

![image](https://github.com/Aetherum17/Patient-Survival/assets/46795020/9b8a8f41-ce19-4c0b-ab1e-90ba87641662)

![image](https://github.com/Aetherum17/Patient-Survival/assets/46795020/b4fe626a-1349-4c0a-b4b8-a6a8898c784c)

![image](https://github.com/Aetherum17/Patient-Survival/assets/46795020/f962a3a3-af89-406f-94a3-dd79fb0fe0fc)

![image](https://github.com/Aetherum17/Patient-Survival/assets/46795020/fbf91529-0fc6-4e6d-b03d-ee8fd5efe482)

![image](https://github.com/Aetherum17/Patient-Survival/assets/46795020/08716ef0-e522-46d8-9ac2-81f86f0a0733)

Next CPHRA was performed on MAP4K2;MEN1, C12orf48, C11orf65.1, RHOG, CH25H, WAPAL.1, 
OTOP2;USH1G.1, CCDC19, TRIOBP, APC2, GRB10.1, GATA4.3, GJA4, KCNG1, BBS5.1, MDFI.1, 
SLC12A4;LCAT.1, CCDC141, SCN4B, TMEM22, PDE6B, NQO1, MS4A3, RDH11.1 and ZFP2 genes, selected 
from methylation data based on their relatively high GKM coefficient for lambda explaining 83,84% of the 
variance, not too high and not too low Hazard Ratio and not too high p-value with additional criteria of not 
preventing GLM model to fit. The result of this CPHRA is presented in image 12, showing MAP4K2-MEN1, 
C12orf48, C11orf65.1, and WAPAL.1, OTOP2-USH1G.1, APC2, GRB10.1, GATA4.3, GJA4, KCNG1, BBS5.1, 
SLC12A4-LCAT.1, SCN4B, TMEM22, PDE6B, RDH11.1 and ZFP2 genes are potentially associated with HCC 
judging by their Hazard Ratio and p-values. To confirm the connection of the methylation of those genes 
with patients' survival, 17 survival curves were plotted in image 13, showing that increased methylation of 
MAP4K2-MEN1, GATA4, TMEM22 ZFP2 are associated with better survival, while increased methylation 
levels of C12orf48, C11orf65.1, BBS5, SCN4B, PDE6B are associated with worse survival. However, neither 
mRNA expression nor mutations in those genes cannot be connected with HCC survival due to having a p-value greater than 0.05, a Hazard ratio equal to 1 or a Hazard ratio being too extreme (Img. 14).

Genes associated with HCC based on literature data were also assessed using CPHRA based on 
methylation data. As can be seen in image 15, a number of them, namely RB1, APC.2, RB1.2, 
“ALG11_ATP7B” and MTOR could potentially be associated with improving or decreasing a patient’s survival. 
Plotted afterwards, survival curves showed an association of increased methylation levels of the APC gene with 
worsened survival, while increased methylation levels of the MTOR gene might be associated with better 
prognosis (Img. 16).

![image](https://github.com/Aetherum17/Patient-Survival/assets/46795020/89d61c0d-2f09-42c6-8536-412bf0b0998d)

![image](https://github.com/Aetherum17/Patient-Survival/assets/46795020/0736208d-7a2d-4f68-80fa-5358130e07f0)
![image](https://github.com/Aetherum17/Patient-Survival/assets/46795020/8a0a50b2-7613-4ddc-84f8-ac036c3e384f)
![image](https://github.com/Aetherum17/Patient-Survival/assets/46795020/c4e278b7-7cad-4908-9873-72538577961e)
![image](https://github.com/Aetherum17/Patient-Survival/assets/46795020/e76f991a-73cd-4fde-a7ee-ef47429f77e7)

![image](https://github.com/Aetherum17/Patient-Survival/assets/46795020/cb23acb0-2c44-4b71-94c3-9cd9c972c011)

![image](https://github.com/Aetherum17/Patient-Survival/assets/46795020/2b3375a9-2dd5-48dc-b983-2015f96bc4a9)

![image](https://github.com/Aetherum17/Patient-Survival/assets/46795020/b32dec59-5fed-4a4c-aade-78a6ed4e4dea)

























