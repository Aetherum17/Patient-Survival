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

Finally, CPHRA was performed on genes with their mutations, selected from mutation data based on their 
relatively high GKM coefficient for lambda, explaining 89.51% of the variance, not too high and not too low 
Hazard Ratio and not too high p-value with additional criteria of not preventing GLM model from fitting. The 
result of this analysis is presented in image 17, highlighting the 5’-prime UTR mutation of the PRKACB gene, 
missense mutation of the DPYD gene, missense mutation of the MUC6 gene, downstream gene mutation of 
OR5AL1 gene, synonymous mutation of STK39 gene, synonymous mutation of F10 gene, 3’-prime UTR 
mutation of the BID gene and missense variant of the CNOT4 gene. Built survival curves on image 18 confirm the 
association of the presence of mutations in PRKACB, MUC6, OR5AL1, STK39, F10, BID and CNOT4 genes 
with decreased patient survival.

The mRNA expression and methylation levels of those genes were also tested via CPHRA, but as can be 
seen from image 19, either their Hazard ratio was equal to 1 or the p-value was larger than 0.05, which 
does not allow to connect those genes to mRNA/methylation data.

Additionally, CPHRA using mutation data was performed on genes, mentioned in HCC literature. Its results 
are presented in image 20, highlighting three potential mutations – Inframe deletion in the TP53 gene, 
Stop gained in the TP53 gene, and Stop gained in the ARID2 gene. However, plotted survival curves were 
only able to confirm the association of the TP53 gene, already found in the analysis of Progeny results (Img. 
21). 

Performed next CPHRA for the TP53 gene based on mRNA expression data showed no association with 
patient’s survival, but the results based on methylation data were more promising, as its hazard ratio value 
was 0.28 with a p-value equal to 0.003 (Img. 22). Survival plot was able to confirm the association of 
methylation of TP53 gene with HCC survival (Img. 23).

![image](https://github.com/Aetherum17/Patient-Survival/assets/46795020/b5c5b7dc-098f-4c69-9956-ab3d33950010)

![image](https://github.com/Aetherum17/Patient-Survival/assets/46795020/a83df64f-2ce0-4dc4-9b3a-b1de87c8a5e0)

![image](https://github.com/Aetherum17/Patient-Survival/assets/46795020/a6e3c92f-dc09-4505-9128-dc9482b877e1)

![image](https://github.com/Aetherum17/Patient-Survival/assets/46795020/23684edc-dcae-40ca-8fbd-348649ea671e)

![image](https://github.com/Aetherum17/Patient-Survival/assets/46795020/8026d41c-9f8a-4264-adb6-ecf3ab4cb02a)

![image](https://github.com/Aetherum17/Patient-Survival/assets/46795020/78058ecb-be48-4456-95a9-846feb4672a3)

![image](https://github.com/Aetherum17/Patient-Survival/assets/46795020/3f9996bf-37f4-4f2c-8df4-afb7955ec5cf)

# Discussion

The performed dataset analysis has uncovered genes, expression, methylation or mutation which might 
affect the overall survival of patients with hepatocellular carcinoma. 

Increased expression of Olfactory receptors, to which genes OR5AP2 and OR11H4 belong, was detected in 
breast carcinoma tissues and cell lines, as well as in pancreatic carcinoma tissues, which allows using of 
them as potential cancer biomarkers and therapeutic targets.[17, 18] Increased expression of long-chain fatty 
acyl CoA synthetases, to which gene ACSM4 belong, was observed in hepatocellular carcinoma tissues, 
compared to the normal liver, giving potential for this gene to be used as a cancer biomarker.[19] Increased 
expression of the PRDM14 gene, a transcriptional regulator which plays a key role in maintaining 
pluripotency of embryonic cells, was also observed in cancer stem cells, making this gene a candidate 
oncogene and a potential therapeutic target.[20] Serine Proteases could also potentially be associated with 
cancer, as it was shown for secreted pseudo-serine protease PRSS35, which acts as a tumour suppressor.[21]
However, in the current analysis, increased expression of a similar protease, PRSS37, was associated with 
poor prognosis. Increased expression of other genes, namely C4orf35, SPDYE4 and CIB4, which respectively 
code calcium-binding protein[22], Speedy E4 protein[23] and calcium & integrin-binding protein[24], was also 
associated with poor patient survival, but no association with cancer was found for those genes in other 
scientific papers, apart from C4orf35 (CABS1) gene, associated with Hereditary Wilms' Tumor.[22]

The performed progeny analysis of signalling pathways, which might be associated with poor survival for 
patients with hepatocellular carcinoma have highlighted “Hypoxia” and “p53” pathways. Indeed, hypoxia-inducible factors play a key role in the development of HCC, as tumour growth prevents diffusion of oxygen 
from the blood to cancer cells, thus creating a hypoxic microenvironment to which malignant cells are 
adapting via the “Hypoxia” pathway, consequently making them a potential target for precise treatment of 
HCC.25 “p53” pathway also plays an important role in the occurrence of HCC, as genes in this pathway 
perform a crucial role in preventing normal cells from malignant transformation via control of their cell 
cycle, apoptosis and metabolism.[26]

In the methylation data a number of genes were also uncovered, increased or decreased methylation 
levels of which might be associated with better survival of patients with HCC. For instance, despite high 
expression levels of the MAP4K2 gene, which are associated with a better prognosis for patients with 
Acute Myeloid Leukemia[27] is suppressed by increased methylation of the MAP4K2 gene, at the same 
time, it also suppresses the expression of the MEN1 gene, which promotes hepatocellular carcinogenesis.[28]
Additionally, increased methylation of SCN4B, leading to its suppression, seems to be an important 
potential biomarker for poor outcomes of HCC disease, as this gene acts as a suppressor of metastasis in 
breast cancer as it prevents cell migration.[29] Suppression of TMEM22 via its methylation also can be 
confirmed as a positive factor for patient survival, as it was shown that overexpression of this gene is 
associated with the growth of renal cell carcinoma cells, and TMEM22 could be a potential therapeutic 
target for this type of cancer.[30] Methylation of INS and IGF2 genes, found as the result of progeny analysis, 
was also checked in literature for the association with cancer, and despite having no significant result for 
the insulin gene, overexpression of the IGF2 gene was observed in many cancers, and it is associated with 
poor patient survival. So, indeed, suppression of the IGF2 gene via its methylation can be associated with 
better outcomes.[31]

However, contradicting results to the published data were also obtained. For example, the C12orf48 gene, 
encoding a PARP1 binding protein, is overexpressed in several cancers, including gastric one, and is said to 
be a potential target for poly(ADP-Ribose) polymerase inhibitors, but in current results, its suppression via 
methylation is associated with poor prognosis.[32, 33] Another contradicting result comes for the PDE6 gene, 
which is expressed in breast cancer tissues but not in normal tissues with only the exception of 
photoreceptors, thus potentially playing a role in breast cancer development, but in current data, its 
suppression via methylation is associated with reduced patient survival.[34] The last contradicting result is 
associated with the GATA4 gene, which has been shown to suppress the development of Hepatocellular Carcinoma cells 
by organising the assembly of a tumour suppressor enhancing module, which inhibits β-catenin 
transcription, but in current data, the suppression of GATA4 via its methylation is associated with a better 
prognosis.[35] 

Methylation of a couple of genes, such as C11orf65, a mitochondrial fission factor 
interaction[36], BBS5, a Bardet-Biedl Syndrome 5 Protein[37] and ZFP2, a Zinc Finger Protein38, was also 
associated with the patient's survival. Suppression of C11orf65 and BBS5 via their methylation was linked 
to poor prognosis, while suppression of the ZFP2 gene was associated with a better prognosis. However, any connections between those genes with cancer in the literature apart from 
C11orf65, associated with Mantle Cell Lymphoma[36] were found.

Additionally, methylation of TP53 and WRAP53 gene, an increase of which is associated with better patient 
survival was also checked in literature sources, and it looks like such an effect, despite the suppression of 
TP53, a tumour-suppressor protein, can be explained by the suppression of the WRAP53 gene, whose function is 
important for telomerase localization and overexpression, which is considered a potential biomarker 
for various cancer types.[39]

Next, the methylation of APC and MTOR genes, known to be associated with the development of 
Hepatocellular carcinoma was also attempted to be confirmed in the literature. Suppression of MTOR via 
its methylation was indeed confirmed to be associated with better patient survivors, as the mTOR pathway 
was shown to be upregulated in HCC tissues.[40] However, the obtained results for the APC gene, telling that 
its suppression via methylation was associated with a worse prognosis did not find confirmation in the 
literature, as methylation of this gene increased in normal tissues and decreased in HCC tissues.[41]

Finally, several genes associated with the overall survival of patients with hepatocellular carcinoma were 
identified in mutation data. To note, there were no mutations detected, that would increase the patient's 
lifespan. Thus, missense mutation of the OR11H4 gene, 5’-prime UTR mutation of the PRKACB gene, 
missense mutation of the MUC6 gene, downstream gene mutation of OR5AL1 gene, synonymous mutation 
of STK39 gene, synonymous mutation of F10 gene, 3’-prime UTR mutation of BID gene, missense variant of 
CNOT4 genes were all associated with reduced patient survival. However, it is hard to speculate about the 
consequence of those mutations and how, in reality, they would affect the lifespan of other patients due to 
the extremely low number of those mutation counts in samples, which can also be seen on survival plots, 
having very steep survival curves with a handful of points on the graph. Additionally, not only mutations of 
genes, but their polymorphisms also play an important role in the development of hepatocellular 
carcinoma, as was shown for the MUC6 gene.42 Indeed, hepatocellular carcinoma is known for its wide 
mutation spectrum, making it even more difficult to find a cure for this disease.[26]

However, if we still would like to confirm the effect of discovered mutations, for example, for Genes, 
discovered from Progeny Analysis, like missense mutation of HIF1A, splice donor variant of PTEN gene, 
Inframe deletion of TP53 gene, stop gained mutation of TP53 gene or missense variant of EI24 gene, we 
could still do this for PTEN and TP53 genes, as the first one is frequently mutated or even deleted in various 
tumours[43], and mutations in the second one are just the most common mutations in HCC that affect its 
progression, prognosis, and result in the reduced immune response to malignant cells, which consequently 
negatively influences the survival of patients with this disease[44].

Additionally, during the performed analysis, separate tests were conducted to determine if mRNA 
expression, methylation or mutations of TERT, CTNNB1, AXIN1, APC, P53, RB1, CCNA2, CCNE1, PTEN, 
ARID1A, ARID2, RPS6KA3, NFE2L2, CCND1, FGF19, VEGFA, MYC, MET, HFE, SERPINA1, G6PC, SLC37A4, 
HMBS, UROD, FAH, ATP7B, GPC3, SPP1, GOLM1, EGF, PDGFB, IGF1, and MTOR genes, mentioned in the 
literature in association with HCC[1,2,3], is important for patient survival. This might be true only for 
methylation of APC and MTOR genes, and mutations of TP53, discussed earlier. To note, the same TP53 
mutations were also obtained from the results of progeny analysis, which could add additional validity to 
those findings.

# Conclusion

As a result of this work, several features from mRNA expression, methylation and mutation data were 
uncovered that might influence the overall survival of patients with Hepatocellular carcinoma. Most of 
those features were already mentioned in the literature in association with HCC, which supports their 
finding. Unfortunately, results for some of the obtained features contradict already known data 
about them, but this potentially could be caused by limitations of the data set analyzed – its small patient 
size, disregard for the molecular subtypes of the disease and missing data for some of the patients.

# References: 

1. Llovet, J.M., Kelley, R.K., Villanueva, A. et al. Hepatocellular carcinoma. Nat Rev Dis Primers 7, 6 
(2021). https://doi.org/10.1038/s41572-020-00240-3
2. Tunissiolli NM, Castanhole-Nunes MMU, Biselli-Chicote PM, Pavarino EC, da Silva RF, da Silva RC, 
Goloni-Bertollo EM. Hepatocellular Carcinoma: a Comprehensive Review of Biomarkers, Clinical 
Aspects, and Therapy. Asian Pac J Cancer Prev. 2017 Apr 1;18(4):863-872. doi: 
10.22034/APJCP.2017.18.4.863. PMID: 28545181; PMCID: PMC5494234.
3. N. Méndez-Sánchez, A. Valencia-Rodríguez, C. Coronel-Castillo, X. Qi. Narrative review of 
hepatocellular carcinoma: from molecular bases to therapeutic approach. Digestive Medicine 
Research. 2021;4:15 doi: http://dx.doi.org/10.21037/dmr-20-116
4. Cerami E, Gao J, Dogrusoz U, Gross BE, Sumer SO, Aksoy BA, Jacobsen A, Byrne CJ, Heuer ML, 
Larsson E, Antipin Y, Reva B, Goldberg AP, Sander C, Schultz N. The cBio cancer genomics portal: an 
open platform for exploring multidimensional cancer genomics data. Cancer Discov. 2012 
May;2(5):401-4. doi: 10.1158/2159-8290.CD-12-0095. Erratum in: Cancer Discov. 2012 
Oct;2(10):960. PMID: 22588877; PMCID: PMC3956037. Gao J, Aksoy BA, Dogrusoz U, Dresdner G, 
Gross B, Sumer SO, Sun Y, Jacobsen A, Sinha R, Larsson E, Cerami E, Sander C, Schultz N. Integrative 
analysis of complex cancer genomics and clinical profiles using the cBioPortal. Sci Signal. 2013 Apr 
2;6(269):pl1. doi: 10.1126/scisignal.2004088. PMID: 23550210; PMCID: PMC4160307. [Internet] 
cBioPortal Dataset “Liver Hepatocellular Carcinoma (TCGA, PanCancer Atlas)”, [cited 2023 Jun 15]. 
Available from: https://www.cbioportal.org/study/summary?id=lihc_tcga_pan_can_atlas_2018
5. R: A Language and Environment for Statistical Computing. [Internet] [cited 2023 Apr 25]. Available 
from: https://www.R-project.org
6. Posit. RStudio. [Internet] [cited 2023 Apr 25]. Available from: https://posit.co/download/rstudiodesktop/
7. rstudioapi: Safely Access the RStudio API [Internet] [cited 2023 Apr 25]. Available from:
https://cran.r-project.org/web/packages/rstudioapi/index.html
8. stringr: Simple, Consistent Wrappers for Common String Operations [Internet] [cited 2023 Apr 25]. 
Available from: https://cran.r-project.org/web/packages/stringr/index.html
9. dplyr: A Grammar of Data Manipulation [Internet] [cited 2023 Apr 25]. Available from: 
https://cran.r-project.org/web/packages/dplyr/index.html
10. progeny: Pathway RespOnsive GENes for activity inference from gene expression [Internet] [cited 2023 Apr 25]. Available from: https://bioconductor.org/packages/release/bioc/html/progeny.html
11. survival: Survival Analysis [Internet] [cited 2023 Apr 25]. Available from: https://cran.rproject.org/web/packages/survival/index.html
12. glmnet: Lasso and Elastic-Net Regularized Generalized Linear Models [Internet] [cited 2023 Apr 25]. 
Available from: https://cran.r-project.org/web/packages/glmnet/index.html
13. tidyverse: Easily Install and Load the 'Tidyverse' [cited 2023 Apr 25]. Available from: https://cran.rproject.org/web/packages/tidyverse/index.html
14. survminer: Drawing Survival Curves using 'ggplot2' [cited 2023 Apr 25]. Available from:
https://cran.r-project.org/web/packages/survminer/index.html
15. Guo Y, Xiao Z, Yang L, Gao Y, Zhu Q, Hu L, Huang D, Xu Q. Hypoxia‑inducible factors in 
hepatocellular carcinoma (Review). Oncol Rep. 2020 Jan;43(1):3-15. doi: 10.3892/or.2019.7397. 
Epub 2019 Nov 1. PMID: 31746396; PMCID: PMC6908932.
16. Kanehisa, M. and Goto, S.; KEGG: Kyoto Encyclopedia of Genes and Genomes. Nucleic Acids Res. 28, 
27-30 (2000). Kanehisa, M; Toward understanding the origin and evolution of cellular organisms. 
Protein Sci. 28, 1947-1951 (2019) Kanehisa, M., Furumichi, M., Sato, Y., Kawashima, M. and 
Ishiguro-Watanabe, M.; KEGG for taxonomy-based analysis of pathways and genomes. Nucleic 
Acids Res. 51, D587-D592 (2023). p53 signaling pathway - Homo sapiens (human) [cited 2023 Apr 
25]. Available from: https://www.genome.jp/pathway/hsa04115
17. Weber Lea, Maßberg Désirée, Becker Christian, Altmüller Janine, Ubrig Burkhard, Bonatz Gabriele, 
Wölk Gerhard, Philippou Stathis, Tannapfel Andrea, Hatt Hanns, Gisselmann Günter. Olfactory 
Receptors as Biomarkers in Human Breast Carcinoma Tissues. Frontiers in Oncology. 2018; vol. 8. 
doi: 10.3389/fonc.2018.00033
18. Chung C, Cho HJ, Lee C, Koo J. Odorant receptors in cancer. BMB Rep. 2022 Feb;55(2):72-80. doi: 
10.5483/BMBRep.2022.55.2.010. PMID: 35168702; PMCID: PMC8891625.
19. Ndiaye H, Liu JY, Hall A, Minogue S, Morgan MY, Waugh MG. Immunohistochemical staining reveals 
differential expression of ACSL3 and ACSL4 in hepatocellular carcinoma and hepatic gastrointestinal 
metastases. Biosci Rep. 2020 Apr 30;40(4):BSR20200219. doi: 10.1042/BSR20200219. PMID: 
32286604; PMCID: PMC7198044.
20. Tracey LJ, Justice MJ. Off to a Bad Start: Cancer Initiation by Pluripotency Regulator PRDM14. 
Trends Genet. 2019 Jul;35(7):489-500. doi: 10.1016/j.tig.2019.04.004. Epub 2019 May 23. PMID: 
31130394; PMCID: PMC6760870.
21. Wang T, Zhou Y, Zhou Z, Zhang P, Yan R, Sun L, Ma W, Zhang T, Shen S, Liu H, Lu H, Ye L, Feng J, 
Chen Z, Zhong X, Wu G, Cai Y, Jia W, Gao P, Zhang H. Secreted protease PRSS35 suppresses 
hepatocellular carcinoma by disabling CXCL2-mediated neutrophil extracellular traps. Nat Commun. 
2023 Mar 18;14(1):1513. doi: 10.1038/s41467-023-37227-z. PMID: 36934105; PMCID: 
PMC10024721.
22. The GeneCards Suite: From Gene Data Mining to Disease Genome Sequence Analyses (PMID: 27322403) Stelzer G, Rosen R, Plaschkes I, Zimmerman S, Twik M, Fishilevich S, Iny Stein T, Nudel R, 
Lieder I, Mazor Y, Kaplan S, Dahary, D, Warshawsky D, Guan - Golan Y, Kohn A, Rappaport N, Safran 
M, and Lancet D Current Protocols in Bioinformatics(2016), 54:1.30.1 - 1.30.33.doi: 10.1002 / cpbi.5 
CABS1 Gene [Internet] [cited 2023 Jun 19]. Available from: https://www.genecards.org/cgibin/carddisp.pl?gene=CABS1
23. The GeneCards Suite: From Gene Data Mining to Disease Genome Sequence Analyses (PMID: 27322403) Stelzer G, Rosen R, Plaschkes I, Zimmerman S, Twik M, Fishilevich S, Iny Stein T, Nudel R, 
Lieder I, Mazor Y, Kaplan S, Dahary, D, Warshawsky D, Guan - Golan Y, Kohn A, Rappaport N, Safran 
M, and Lancet D Current Protocols in Bioinformatics(2016), 54:1.30.1 - 1.30.33.doi: 10.1002 / cpbi.5 
SPDYE4 Gene [Internet] [cited 2023 Jun 19]. Available from: https://www.genecards.org/cgibin/carddisp.pl?gene=SPDYE4 
24. The GeneCards Suite: From Gene Data Mining to Disease Genome Sequence Analyses (PMID: 
27322403) Stelzer G, Rosen R, Plaschkes I, Zimmerman S, Twik M, Fishilevich S, Iny Stein T, Nudel R, 
Lieder I, Mazor Y, Kaplan S, Dahary, D, Warshawsky D, Guan - Golan Y, Kohn A, Rappaport N, Safran 
M, and Lancet D Current Protocols in Bioinformatics(2016), 54:1.30.1 - 1.30.33.doi: 10.1002 / cpbi.5 
CIB4 Gene [Internet] [cited 2023 Jun 19]. Available from: https://www.genecards.org/cgibin/carddisp.pl?gene=CIB4
25. Guo Y, Xiao Z, Yang L, Gao Y, Zhu Q, Hu L, Huang D, Xu Q. Hypoxia‑inducible factors in 
hepatocellular carcinoma (Review). Oncol Rep. 2020 Jan;43(1):3-15. doi: 10.3892/or.2019.7397. 
Epub 2019 Nov 1. PMID: 31746396; PMCID: PMC6908932.
26. Link T, Iwakuma T. Roles of p53 in extrinsic factor-induced liver carcinogenesis. Hepatoma Res. 
2017;3:95-104. doi: 10.20517/2394-5079.2017.07. Epub 2017 Jun 6. PMID: 30123836; PMCID: 
PMC6097626.
27. Bai Z, Yao Q, Sun Z, Xu F, Zhou J. Prognostic Value of mRNA Expression of MAP4K Family in Acute 
Myeloid Leukemia. Technol Cancer Res Treat. 2019 Jan 1;18:1533033819873927. doi: 
10.1177/1533033819873927. PMID: 31522654; PMCID: PMC6747867.
28. Xu B, Li SH, Zheng R, Gao SB, Ding LH, Yin ZY, Lin X, Feng ZJ, Zhang S, Wang XM, Jin GH. Menin 
promotes hepatocellular carcinogenesis and epigenetically up-regulates Yap1 transcription. Proc 
Natl Acad Sci U S A. 2013 Oct 22;110(43):17480-5. doi: 10.1073/pnas.1312022110. Epub 2013 Oct 
7. PMID: 24101467; PMCID: PMC3808599.
29. Bon E, Driffort V, Gradek F, Martinez-Caceres C, Anchelin M, Pelegrin P, Cayuela ML, MarionneauLambot S, Oullier T, Guibon R, Fromont G, Gutierrez-Pajares JL, Domingo I, Piver E, Moreau A, 
Burlaud-Gaillard J, Frank PG, Chevalier S, Besson P, Roger S. SCN4B acts as a metastasis-suppressor 
gene preventing hyperactivation of cell migration in breast cancer. Nat Commun. 2016 Dec 
5;7:13648. doi: 10.1038/ncomms13648. PMID: 27917859; PMCID: PMC5150224.
30. Dobashi S, Katagiri T, Hirota E, Ashida S, Daigo Y, Shuin T, Fujioka T, Miki T, Nakamura Y. 
Involvement of TMEM22 overexpression in the growth of renal cell carcinoma cells. Oncol Rep. 2009 Feb;21(2):305-12. PMID: 19148500.
31. Livingstone C. IGF2 and cancer. Endocr Relat Cancer. 2013 Oct 24;20(6):R321-39. doi: 10.1530/ERC13-0231. PMID: 24080445.
32. Lin L, Li H, Shi D, Liu Z, Wei Y, Wang W, Wu D, Li B, Guo Q. Depletion of C12orf48 inhibits gastric 
cancer growth and metastasis via up-regulating Poly r(C)-Binding Protein (PCBP) 1. BMC Cancer. 
2022 Jan 31;22(1):123. doi: 10.1186/s12885-022-09220-0. PMID: 35100974; PMCID: PMC8802463.
33. Paturel A, Hall J, Chemin I. Poly(ADP-Ribose) Polymerase Inhibition as a Promising Approach for 
Hepatocellular Carcinoma Therapy. Cancers (Basel). 2022 Aug 5;14(15):3806. doi: 
10.3390/cancers14153806. PMID: 35954469; PMCID: PMC9367559.
34. Dong H, Claffey KP, Brocke S, Epstein PM. Expression of phosphodiesterase 6 (PDE6) in human 
breast cancer cells. Springerplus. 2013 Dec 18;2:680. doi: 10.1186/2193-1801-2-680. PMID: 24683528; PMCID: PMC3967736.
35. Lu F, Zhou Q, Liu L, Zeng G, Ci W, Liu W, Zhang G, Zhang Z, Wang P, Zhang A, Gao Y, Yu L, He Q, Chen 
L. A tumor suppressor enhancing module orchestrated by GATA4 denotes a therapeutic 
opportunity for GATA4 deficient HCC patients. Theranostics. 2020 Jan 1;10(2):484-497. doi: 
10.7150/thno.38060. PMID: 31903133; PMCID: PMC6929984.
36. The GeneCards Suite: From Gene Data Mining to Disease Genome Sequence Analyses (PMID: 27322403) Stelzer G, Rosen R, Plaschkes I, Zimmerman S, Twik M, Fishilevich S, Iny Stein T, Nudel R, 
Lieder I, Mazor Y, Kaplan S, Dahary, D, Warshawsky D, Guan - Golan Y, Kohn A, Rappaport N, Safran 
M, and Lancet D Current Protocols in Bioinformatics(2016), 54:1.30.1 - 1.30.33.doi: 10.1002 / cpbi.5 
C11orf65 Gene [Internet] [cited 2023 Jun 19]. Available from: https://www.genecards.org/cgibin/carddisp.pl?gene=C11orf65
37. The GeneCards Suite: From Gene Data Mining to Disease Genome Sequence Analyses (PMID: 
27322403) Stelzer G, Rosen R, Plaschkes I, Zimmerman S, Twik M, Fishilevich S, Iny Stein T, Nudel R, 
Lieder I, Mazor Y, Kaplan S, Dahary, D, Warshawsky D, Guan - Golan Y, Kohn A, Rappaport N, Safran 
M, and Lancet D Current Protocols in Bioinformatics(2016), 54:1.30.1 - 1.30.33.doi: 10.1002 / cpbi.5 
BBS5 Gene [Internet] [cited 2023 Jun 19]. Available from: https://www.genecards.org/cgibin/carddisp.pl?gene=BBS5
38. The GeneCards Suite: From Gene Data Mining to Disease Genome Sequence Analyses (PMID: 27322403) Stelzer G, Rosen R, Plaschkes I, Zimmerman S, Twik M, Fishilevich S, Iny Stein T, Nudel R, 
Lieder I, Mazor Y, Kaplan S, Dahary, D, Warshawsky D, Guan - Golan Y, Kohn A, Rappaport N, Safran 
M, and Lancet D Current Protocols in Bioinformatics(2016), 54:1.30.1 - 1.30.33.doi: 10.1002 / cpbi.5 
ZFP2 Gene [Internet] [cited 2023 Jun 19]. Available from: https://www.genecards.org/cgibin/carddisp.pl?gene=ZFP2
39. Gadelha, R.B.; Machado, C.B.; Pessoa, F.M.C.d.P.; Pantoja, L.d.C.; Barreto, I.V.; Ribeiro, R.M.; de 
Moraes Filho, M.O.; de Moraes, M.E.A.; Khayat, A.S.; Moreira-Nunes, C.A. The Role of WRAP53 in 
Cell Homeostasis and Carcinogenesis Onset. Curr. Issues Mol. Biol. 2022, 44, 5498-5515. 
https://doi.org/10.3390/cimb44110372
40. Ferrín G, Guerrero M, Amado V, Rodríguez-Perálvarez M, De la Mata M. Activation of mTOR 
Signaling Pathway in Hepatocellular Carcinoma. Int J Mol Sci. 2020 Feb 13;21(4):1266. doi: 
10.3390/ijms21041266. PMID: 32070029; PMCID: PMC7072933.
41. Csepregi A, Röcken C, Hoffmann J, Gu P, Saliger S, Müller O, Schneider-Stock R, Kutzner N, Roessner 
A, Malfertheiner P, Ebert MP. APC promoter methylation and protein expression in hepatocellular 
carcinoma. J Cancer Res Clin Oncol. 2008 May;134(5):579-89. doi: 10.1007/s00432-007-0321-y. 
Epub 2007 Nov 1. PMID: 17973119; PMCID: PMC2757596.
42. Lee HL, Chien YC, Wang HL, Hua CH, Liu LC, Wu GW, Bai LY, Yang SF, Yu YL. Analysis of MUC6 
Genetic Variants on the Clinicopathologic Characteristics of Patients with Hepatocellular Carcinoma. 
J Cancer. 2022 Sep 6;13(11):3251-3257. doi: 10.7150/jca.75754. PMID: 36118520; PMCID: 
PMC9475359.
43. Shearn CT, Petersen DR. Understanding the tumor suppressor PTEN in chronic alcoholism and 
hepatocellular carcinoma. Adv Exp Med Biol. 2015;815:173-84. doi: 10.1007/978-3-319-09614-
8_10. PMID: 25427907.
44. Long J, Wang A, Bai Y, Lin J, Yang X, Wang D, Yang X, Jiang Y, Zhao H. Development and validation of 
a TP53-associated immune prognostic model for hepatocellular carcinoma. EBioMedicine. 2019 
Apr;42:363-374. doi: 10.1016/j.ebiom.2019.03.022. Epub 2019 Mar 16. PMID: 30885723; PMCID: 
PMC6491941.































