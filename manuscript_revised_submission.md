# Plasma 1H-NMR metabolomics identifies amino acid and carbohydrate pathway alterations associated with obesity

**Karine Cavalcanti Mauricio Sena-Evangelista**^1,2^, **Folurunsho Bright Omage**^3^, **Erik Sobrinho Braga**^3^, **Lucas Gelain Martins**^3^, **Paula Emilia Nunes Ribeiro Bellot**^2^, **Adriano Carlos de Souza Junior**^2^, **Lucia Fatima Campos Pedrosa**^1,2^, **Dirce Maria Lobo Marchioni**^4^, **Fernando Barbosa Jr.**^5^, **Severina Carla Vieira Cunha Lima**^1,2^, **Clelia Oliveira Lyra**^1,2^, **Ljubica Tasic**^3^

1. Department of Nutrition, Center for Health Sciences, Federal University of Rio Grande do Norte, Natal, Rio Grande do Norte, Brazil
2. Graduate Program in Nutrition, Center for Health Sciences, Federal University of Rio Grande do Norte, Natal, Rio Grande do Norte, Brazil
3. Biological Chemistry Laboratory, Department of Organic Chemistry, Institute of Chemistry, Universidade Estadual de Campinas (UNICAMP), Campinas, Sao Paulo, Brazil
4. Department of Nutrition, School of Public Health, University of Sao Paulo, Sao Paulo, Sao Paulo, Brazil
5. Department of Clinical Analyses, Toxicology and Food Sciences, School of Pharmaceutical Sciences of Ribeirao Preto of the University of Sao Paulo, Ribeirao Preto, Sao Paulo, Brazil

*Corresponding author:
E-mail: karine.sena@ufrn.br (KCMSE)

---

## Abstract

**Scope:** Obesity perturbs intermediary metabolism and elevates cardiometabolic risk. We profiled circulating metabolites to distinguish adults with BMI >= 30 kg/m^2^ from those with BMI < 30 kg/m^2^ and to explore their utility in diagnosis and nutritional guidance.

**Methods:** In a cross-sectional sample of 72 adults (36 per group), untargeted metabolomics by proton nuclear magnetic resonance (^1^H-NMR) was performed on plasma. Data were log-transformed and normalized. Group differences were assessed with Welch's t-tests and Benjamini-Hochberg false discovery rate (FDR) control. Significant metabolites were contextualized using pathway analysis and a metabolite-protein-disease network. Supervised models were trained with stratified splits and comprehensive validation including repeated cross-validation, permutation testing, and bootstrap confidence intervals to ensure robustness; performance metrics were computed on a hold-out test set.

**Results:** Thirty-four metabolites differed significantly between groups after FDR correction (p_adj < 0.05), including elevated levels of lysine, adenine, kynurenic acid, and alanine, and decreased levels of D-fucose, myo-inositol, and asparagine in individuals with BMI >= 30 kg/m^2^. Enrichment analysis pointed to amino-acid and carbohydrate metabolism (*e.g.*, phenylalanine/tyrosine/tryptophan biosynthesis; glycolysis/gluconeogenesis). Using repeated 10x5-fold cross-validation as the primary performance metric, support vector machine (SVM) achieved the best discrimination (CV accuracy: 72.4% +/- 9.1%, CV AUC: 0.735 +/- 0.136), significantly exceeding chance performance (permutation test p < 0.01). Lysine and quinolinic acid emerged as the most discriminatory metabolites.

**Conclusions:** Plasma ^1^H-NMR metabolomics reveals coordinated shifts in amino acid and carbohydrate pathways. Combined with machine learning, these findings may improve risk stratification and guide personalized nutritional strategies in obesity management.

**Keywords:** NMR-Metabolomics, body mass index, machine learning models, amino acid metabolism, metabolic biomarkers.

---

## Introduction

Obesity is a chronic, progressive, and multifactorial disease characterized by excessive fat deposits that impair health. It is estimated that 650 million adults are living with obesity, equivalent to 15.8% of the world population [1]. This condition is predominantly associated with lifestyle, including reduced physical activity and unhealthy eating habits. Furthermore, several genetic, hereditary, psychological, cultural, and ethnic factors also play significant roles in its manifestation and progression [2].

Inflammation, dyslipidemia, and insulin resistance resulting from obesity can promote atherogenesis through several mechanisms, such as the oxidation of lipoproteins [3]. Consequently, obesity exacerbates the risk for some metabolic and cardiovascular diseases (CVD) [4]. Thus, it is crucial to recognize the risk factors in the early stages of cardiovascular damage, as well as personalized risk stratification to implement prevention strategies focusing on weight management to reduce the impact of CVD [5].

Body mass index (BMI) is a practical tool for estimating adiposity in epidemiological and screening contexts; however, it has limitations in accurately reflecting body composition. Although a BMI >= 30 kg/m^2^ is commonly used to define obesity, clinical diagnosis requires confirmation through direct measures of body fat or additional anthropometric criteria [6]. Interestingly, the lowest mortality has been observed at a BMI around 25 kg/m^2^, despite the optimal range being 18.5-24.9 kg/m^2^ [7]. This paradox may be explained by the recognition that obesity phenotypes, rather than BMI alone, better predict health risks, highlighting the need for more comprehensive assessments beyond BMI [8].

Therefore, to improve understanding of obesity-associated health risks, the characterization of the metabolites can provide insights into the mechanisms that associate obesity with cardiometabolic consequences [9]. In this context, metabolomics has allowed the detection and measurement of changes in metabolite levels in response to a genetic variation, physiological, or pathological condition. This approach can lead to the identification of new biomarkers and clarify a diversity of disease mechanisms [10].

Metabolomics has proven useful in identifying metabolite alterations associated with obesity by comparing obese and healthy individuals [11]. These changes are linked to metabolic diseases and increased cardiovascular risk. Lipidomic analyses further demonstrate that obesity is associated with distinct lipid profile shifts, including correlations between BMI and specific phosphatidylcholine and lysophosphatidylcholine species [12]. Elevated lipid metabolite levels in individuals with BMI >= 30 kg/m^2^ support the existence of a metabolically distinct obesity phenotype [13].

Despite advances in metabolomics, there is still no clear consensus on which metabolites reliably indicate obesity-related metabolic dysfunction. To address this, the present study applies proton nuclear magnetic resonance spectroscopy (^1^H-NMR), a reproducible and non-destructive method, combined with machine learning to analyze plasma metabolite profiles. The goal is to identify metabolic signatures that differentiate individuals with BMI < 30 from those with BMI >= 30 kg/m^2^, aiming to support personalized risk assessment and guide targeted clinical and nutritional strategies to prevent cardiometabolic complications.

---

## Material and Methods

### Study Design and Population

This study is part of the multicenter BRAZUCA study (Brazilian Usual Consumption Assessment). It is a cross-sectional analysis of a subsample of adults and older individuals from the population-based research titled "Food Insecurity, Health, and Nutrition Conditions in the Adult and Elderly Population of a Northeastern Brazilian Capital: BRAZUCA Study - Natal." Both studies were approved by the Research Ethics Committees of the University of Sao Paulo and the Federal University of Rio Grande do Norte (CAAE: 96294718.4.1001.5421 - Ethical Approval Number 3.789.491, and CAAE: 96294718.4.2001.5292 - Ethical Approval Number 3.301.377, respectively).

The sampling plan considered a two-stage cluster probabilistic sample. A subsample of 344 individuals was planned, and data collection occurred starting in June 2019, and was discontinued in March 2020 due to the international public health emergency caused by COVID-19, resulting in a subsample of 112 participants. Exclusion criteria were pregnant and lactating women, illicit drug users, individuals who underwent chemotherapy and/or radiotherapy in the last 6 months, and those unable to respond to the research questions. The final sample size for this study was 72 individuals, matched by age and sex to reduce the likelihood of bias errors, divided into two groups: people with obesity (n = 36; BMI >= 30 kg/m^2^) and without obesity (n = 36; BMI < 30 kg/m^2^).

### Data Collection, Anthropometric, and Biochemical Parameters

Data collection was performed in households using a standardized, revised questionnaire. This instrument covered sociodemographic data (age, sex, race, income), lifestyle factors (physical activity, smoking, and alcohol use), and diagnosis of type 2 diabetes, arterial hypertension, dyslipidemias, and cardiovascular events (heart attack or stroke). The Global Risk Score was calculated and classified according to the Brazilian Dyslipidemia Guideline [14]. Blood pressure was measured using an Omron (HEM-7122, Kyoto, Japan) automatic arm monitor. Systolic and diastolic blood pressure were classified as elevated when they ranged >= 140 mmHg and >= 90 mmHg, respectively [15].

Body weight was measured using an electronic scale, with a capacity of 200 kg and 50 g precision (Lider P-200 M, Aracatuba, Sao Paulo, Brazil). Height measurement was performed with a portable stadiometer (Avanutri, Tres Rios, Rio de Janeiro, Brazil), 0.1 cm precision and a non-slip base. BMI cutoff was adopted according to the World Health Organization (WHO) [16]. Waist circumference was measured with a non-extendable ergonomic tape measure from Cescorf (Porto Alegre, Rio Grande do Sul, Brazil). We considered the cutoff points proposed by the WHO (2000) [16]. The visceral adiposity index (VAI) was calculated considering the cutoff points proposed previously [16].

Blood was collected by peripheral vein puncture in the morning after 10-12 hours of fasting. In tubes without anticoagulants, 4 mL of blood was allocated for biochemical analysis, and 4 mL was placed in EDTA tubes for metabolomics. Enzymatic methods were used to analyze fasting glycemia, total cholesterol, and triacylglycerols (TAG). A homogeneous enzymatic colorimetric assay measured HDL-c levels. An immunoassay was applied to determine insulin using the sandwich technique, and high-sensitivity C-reactive protein (hs-CRP) by immunoturbidimetric assay. All these analyses were performed using an automated device COBAS 6000 - Roche Professional Diagnostics (Risch-Rotkreuz, Switzerland). Low-density lipoprotein cholesterol (LDL-c) values were obtained using the Friedewald formula [14].

The cutoff point adopted for fasting glycemia is recommended by the Brazilian Diabetes Society [18]. HDL-c values were considered low when < 40 mg/dL. TAG and total cholesterol were considered high when > 150 mg/dL and > 190 mg/dL, respectively. Based on the cardiovascular risk category defined by the Global Risk Score calculation, the LDL-c targets were defined as proposed by the Brazilian Dyslipidemia Guideline [14]. LDL-c values were considered high when they were above the target for each individual. The cutoff points for hs-CRP for cardiovascular risk indication were as recommended by Ridker (2007) [19]. The Homeostatic Model Assessment of insulin resistance (HOMA-IR) was calculated and considered the reference point as suggested by Geloneze et al. (2009) [20].

### Metabolomics Profile Analysis

^1^H-NMR metabolomics was applied to 400 uL of plasma samples by adding 100 uL of deuterium oxide solvent (D~2~O, 99.9% with 0.03% trimethylsilyl propionic acid, TSP, from Sigma-Aldrich, Andover, Maryland, United States of America). Samples were then centrifuged at 12,000 x g for 2 min at 4C, and transferred to 5 mm tubes. The ^1^H-NMR high-resolution (noesypr1d), T~2~-edited (cpmgpr1d), and diffusion-edited (stebpgp1s191d) spectra were recorded using the Bruker AVANCE III 600 MHz instrument (Bruker Biospin, Karlsruhe, Germany) with the Triple Resonance Broadband Inverse (TBI) probe at 25C. Two-dimensional NMR total correlation spectroscopy (TOCSY) data were recorded using the mlevphpr pulse sequence, and heteronuclear quantum correlation spectral data were obtained using the hsqcedetgpsp pulse sequence (Supporting Information Fig S1). Scalar couplings were evaluated from J-RESolved spectroscopic experiments using the jresgpprqf pulse sequence. Metabolites were assigned based on chemical shifts, coupling constants, 2D NMR spectral features, and database searches, such as the Human Metabolome Database and BioMagResBank [21,22].

### Statistical Analysis and Pathway Exploration

Categorical variables were analyzed using Chi-square tests, and data integrity was verified before metabolomic analysis. Metabolite data were processed using the MetaboAnalyst platform, including log~10~ transformation, normalization, and constant sum scaling to correct for concentration differences [23]. Principal Component Analysis (PCA), Partial Least Squares Discriminant Analysis (PLS-DA), and Orthogonal PLS-DA (OPLS-DA) were applied to identify group separation and relevant metabolic features. Variable Importance in Projection (VIP) scores were calculated from the PLS-DA model [23], and statistical significance was defined as p < 0.05. Key metabolites are presented in the Supporting Information (Fig S2 and Table S1).

Following preprocessing and normalizing, log~2~-transformed metabolite intensities were compared between groups. A Welch's t-test (two-sided) was applied to account for potentially unequal variances and sample sizes. The resulting raw p-values were then adjusted using the Benjamini-Hochberg false discovery rate (FDR) procedure to mitigate the risk of multiple-testing errors. Metabolites were deemed significantly altered if they (i) satisfied p_adj < 0.05 and (ii) exhibited |log2 fold change (log2FC)| >= 0.1. The log2FC was calculated by subtracting the mean log2 intensity of each group. An additional manual inspection of borderline cases was performed to identify metabolites potentially overlooked by a single threshold. All significantly altered metabolites (p_adj < 0.05) were mapped to biochemical pathways Human Metabolome Database [21,22] for functional interpretation. Pathway impact scores were also computed, using centrality-based or betweenness-based metrics to pinpoint which pathways had the greatest potential biological effect. Pathways with an adjusted p-value < 0.05 were considered significantly perturbed using the Metaboanalyst 6.0 Pathway Analyst module [24].

Each significant metabolite was queried against the Search Tool for Interacting Chemicals (STITCH) database [25] to retrieve known or predicted protein interaction partners with a confidence score >= 0.70, and an organism filter (Homo sapiens) was applied to ensure higher specificity. The Human Metabolome Database [21,22] was used to find documented links between identified metabolites and reported diseases or phenotypes. Disease and Gene Network (DisGeNET) [26] provided curated associations between each protein and disease. Overlapping or convergent disease categories suggested common pathological mechanisms. Interactions were integrated in Cytoscape [27], generating a multi-level graph where nodes represented metabolites, proteins, and diseases, and edges indicated known functional or literature-based relationships.

### Machine Learning Classification

A supervised-learning pipeline was implemented to classify obesity status (BMI >= 30 vs < 30 kg/m^2^). **Critically, data were randomly split into training (70%, n=50) and test (30%, n=22) sets BEFORE any feature selection to prevent information leakage.** All statistical tests for identifying significant metabolites were performed exclusively on the training set, and the resulting feature set was then applied to both training and test data. Scaling was fit on the training set and applied to the test set.

Nine classification algorithms were evaluated: logistic regression, support vector machine (SVM), decision tree, k-nearest neighbors (KNN), naive Bayes, multi-layer perceptron (MLP), random forest, gradient boosting, and XGBoost. All models were implemented using scikit-learn (v1.0+) and XGBoost (v1.5+) with explicit hyperparameters documented in Supplementary Table S2 for reproducibility. Random seeds were fixed (random_state=42) throughout the analysis.

### Model Validation and Robustness Assessment

<mark>To address potential overfitting concerns inherent to small sample sizes, we implemented a comprehensive validation strategy. First, repeated stratified cross-validation served as the primary performance metric, derived from 10 repetitions of 5-fold stratified cross-validation on the full dataset (n=72), yielding 50 independent performance estimates per model. This approach provides robust mean estimates and standard deviations that account for sampling variability. Second, permutation testing with 1,000 iterations assessed whether model performance exceeded chance levels by randomly shuffling class labels, retraining models, and comparing the resulting null distribution against actual performance; a p-value < 0.05 indicates performance significantly better than random classification. Third, bootstrap confidence intervals for the held-out test set were computed using 1,000 bootstrap resamples with replacement, providing uncertainty estimates for all reported metrics. Fourth, learning curves plotting training and validation accuracy as a function of training set size (20%, 40%, 60%, 80%, 100%) were used to diagnose overfitting.</mark>

The final model was selected based on repeated cross-validation performance. We report accuracy, precision, recall, F1-score, Matthews correlation coefficient (MCC), and area under the receiver operating characteristic curve (AUC-ROC). All code and data are publicly available at https://github.com/omagebright/obesity.

SHapley Additive exPlanations (SHAP) were applied, enabling us to see how each metabolite shifts individual predictions. SHAP summary plots can offer sample-level and global insights, pinpointing whether particular metabolites (features) consistently drive the probability toward or away from the BMI >= 30 kg/m^2^ classification.

Metabolomic data preprocessing involved multiple quality control steps: (i) spectral alignment using the TSP reference peak, (ii) removal of water and solvent regions, (iii) binning to reduce dimensionality while preserving metabolite signals, (iv) normalization to total spectral area to account for concentration differences, (v) log~10~ transformation to achieve normal distribution, and (vi) Pareto scaling to give appropriate weight to both high and low-abundance metabolites. Quality control samples were interspersed throughout the analysis to monitor instrumental drift and ensure reproducibility (coefficient of variation < 15% for all identified metabolites).

---

## Results

### Population Characteristics

No significant differences were observed in the distribution of sex, age, smoking, or alcohol consumption between the two groups (all p > 0.05). Arterial hypertension, Type 2 diabetes, and dyslipidemias were more prevalent in the BMI >= 30 kg/m^2^ group, although the differences were not statistically significant. People with obesity had significantly higher waist circumference (p < 0.001). The BMI < 30 kg/m^2^ group had a higher proportion of people classified as overweight, and the BMI >= 30 kg/m^2^ group showed significantly high TAG levels (p = 0.009). For the Global Risk Score, people with obesity had significantly higher proportions of high or very high cardiovascular risk (p = 0.031) (Table 1).

**Table 1.** Characteristics of the individuals with a BMI >= 30 kg/m^2^ and BMI < 30 kg/m^2^ groups

| Variables | BMI < 30 kg/m^2^ (n = 36) | BMI >= 30 kg/m^2^ (n = 36) | p-value |
|-----------|---------------------------|----------------------------|---------|
| **Sex** | | | |
| Male | 12 (33.3) | 12 (33.3) | 1.000 |
| Female | 24 (66.7) | 24 (66.7) | |
| **Age (years)** | | | |
| > 60 | 21 (58.3) | 22 (61.1) | 0.157 |
| < 60 | 15 (41.7) | 14 (38.9) | |
| **Cigarette Smoking** | | | |
| Never smoker | 22 (61.1) | 23 (63.8) | 0.817 |
| Current smoker | 14 (38.9) | 13 (36.2) | |
| **Alcohol use** | | | |
| Never drink | 22 (61.1) | 22 (61.1) | 1.000 |
| Drink use | 14 (38.9) | 14 (38.9) | |
| **Waist circumference**^a^ | | | |
| Low risk of CVD | 8 (24.2) | 0 (0.0) | < 0.001 |
| Increased risk of CVD | 14 (42.4) | 0 (0.0) | |
| Substantial increased risk | 11 (33.3) | 33 (100.0) | |
| **High TAG** | 12 (33.4) | 23 (63.9) | 0.009 |
| **Global Risk Score** | | | |
| Low CVD risk | 10 (27.7) | 3 (8.3) | 0.031 |
| Intermediary CVD risk | 5 (13.8) | 2 (5.5) | |
| High CVD risk | 21 (58.3) | 31 (86.2) | |

*Abbreviation: BMI, body mass index; CVD, cardiovascular disease; TAG, triacylglycerol. Values are represented by frequency (percentage). ^a^Waist circumference data missing for 3 participants in each group (n=33 per group for this variable).*

### Metabolomic Profiling

Metabolomic profiling revealed significant differences in metabolite levels between groups. Following proper methodology where feature selection was performed on the training set only (n=50), thirty-four metabolites were identified as significantly different (p_adj < 0.05) after FDR adjustment (full list in Supplementary Table S3). The BMI >= 30 kg/m^2^ group exhibited higher levels of L-lysine (log2FC = 0.46, p_adj = 0.004), adenine (log2FC = 0.89, p_adj = 0.010), quinolinic acid (log2FC = 0.65, p_adj = 0.010), L-alanine (log2FC = 0.24, p_adj = 0.010), and kynurenic acid (log2FC = 0.64, p_adj = 0.010). Conversely, the BMI >= 30 kg/m^2^ group showed lower levels of D-fucose (log2FC = -1.97, p_adj = 0.010), myo-inositol (log2FC = -0.80, p_adj = 0.035), L-asparagine (log2FC = -0.66, p_adj = 0.029), and cadaverine (log2FC = -0.67, p_adj = 0.032) (Table 2).

**Table 2.** Top differential metabolites identified between BMI >= 30 kg/m^2^ and BMI < 30 kg/m^2^ groups (training set, n=50)

| Metabolites | BMI >= 30 (Mean) | BMI < 30 (Mean) | Trend | log2FC | Cohen's d | p_adj |
|-------------|------------------|-----------------|-------|--------|-----------|-------|
| L-Lysine | 0.00503 | 0.00365 | Up | 0.46 | 1.30 | 0.004 |
| Adenine | 0.00103 | 0.00056 | Up | 0.89 | 1.17 | 0.010 |
| Quinolinic Acid | 0.00260 | 0.00165 | Up | 0.65 | 1.12 | 0.010 |
| D-Fucose | 0.00168 | 0.00657 | Down | -1.97 | -1.11 | 0.010 |
| Kynurenic Acid | 0.00171 | 0.00110 | Up | 0.64 | 1.07 | 0.010 |
| L-Alanine | 0.00528 | 0.00448 | Up | 0.24 | 1.07 | 0.010 |
| Benzoic Acid | 0.00117 | 0.00070 | Up | 0.73 | 1.03 | 0.011 |
| Hypoxanthine | 0.00059 | 0.00031 | Up | 0.92 | 1.02 | 0.011 |
| Myo-Inositol | 0.00164 | 0.00286 | Down | -0.80 | -0.79 | 0.035 |
| L-Asparagine | 0.00309 | 0.00490 | Down | -0.66 | -0.84 | 0.029 |

*Abbreviation: FC, fold change. Up trend, relatively higher levels of metabolites present in the BMI >= 30 kg/m^2^ group; Down trend, relatively lower levels. p_adj, p-value adjusted for FDR using Benjamini-Hochberg method. Full list of 34 significant metabolites in Supplementary Table S3. Values are relative ASICS intensities.*

<mark>Multivariate and pathway-level analyses corroborate a distinct metabolic signature in individuals with BMI >= 30 kg/m^2^. Univariate testing is summarized in the volcano plot (Fig 1a), which highlights 34 significantly altered metabolites with L-lysine showing the largest effect size (Cohen's d = 1.30). A three-dimensional PCA score plot (Fig 1b) shows partial separation between the obese (BMI >= 30 kg/m^2^) and non-obese (BMI < 30 kg/m^2^) groups, with the first three principal components explaining 73.4%, 13.6% and 8.6% of the total variance, respectively. Pathway enrichment and topology analysis (Fig 1c) pinpointed lysine degradation and tryptophan metabolism as the most significantly affected pathways, followed by phenylalanine, tyrosine, and tryptophan biosynthesis, pyruvate metabolism, and inositol phosphate metabolism.</mark>

<mark>

**Figure 1. Differential Metabolite Analysis**

![Figure 1](figures/Figure1.png)

**Fig 1.** Differential metabolite analysis between BMI >= 30 kg/m^2^ and BMI < 30 kg/m^2^ groups. **(a)** Volcano plot showing 34 significantly altered metabolites (FDR < 0.05). Red points indicate metabolites elevated in obesity; blue points indicate metabolites reduced in obesity. L-Lysine shows the largest effect size (Cohen's d = 1.30). **(b)** Principal Component Analysis (3D) demonstrating partial separation between BMI groups. PC1, PC2, and PC3 explain 73.4%, 13.6%, and 8.6% of the total variance, respectively. **(c)** Pathway enrichment analysis highlighting lysine degradation and tryptophan metabolism as the most significantly affected pathways.

</mark>

### Protein-Metabolite-Disease Network

A protein-metabolite-disease interaction network was constructed to reveal the potential functional relationships among the significantly altered metabolites. Interactions between the metabolites and proteins were obtained from the STITCH database with a confidence threshold >= 0.70. Representative metabolites from the significant panel (lysine, leucine, alanine, glutamic acid, and isoleucine) were mapped to distinct proteins including: alanine aminotransferase (ALT), aspartate aminotransferase (AST), branched-chain aminotransferase 1 (BCAT1), branched-chain alpha-keto acid dehydrogenase complex (BCKDHA), tyrosine aminotransferase (TAT), leucyl-tRNA synthetase (LARS), lysyl-tRNA synthetase (KARS), glutamate dehydrogenase (GLUD1), and others involved in amino acid metabolism and energy production pathways (Supplementary Fig S2).

### Machine Learning Classification

We evaluated nine machine learning algorithms for distinguishing individuals with BMI >= 30 kg/m^2^ from those with BMI < 30 kg/m^2^ using the panel of 34 differentially abundant metabolites identified from the training set. <mark>Repeated 10×5-fold stratified cross-validation served as the primary performance metric, with permutation testing to assess statistical significance.</mark>

<mark>SVM achieved the highest cross-validation accuracy of 72.4% ± 9.1% (F1: 0.746 ± 0.078; AUC: 0.735 ± 0.136), followed by Random Forest at 70.6% ± 10.7% (Table 3, Fig 2a). The moderate standard deviations reflect expected variability given the sample size and represent realistic performance bounds rather than model instability. Permutation testing with 1,000 iterations confirmed that seven of nine models significantly exceeded chance-level classification (p < 0.05), with SVM, Random Forest, XGBoost, and Naive Bayes all achieving p ≤ 0.003 (Fig 2c). Decision Tree and Neural Network did not reach significance (p = 0.098 and p = 0.113, respectively), likely reflecting overfitting tendencies with limited training data.</mark>

**Table 3.** Machine learning model performance for obesity classification

| Model | CV Accuracy | CV F1-Score | CV AUC | Permutation p |
|-------|-------------|-------------|--------|---------------|
| SVM | 0.724 ± 0.091 | 0.746 ± 0.078 | 0.735 ± 0.136 | 0.002 |
| Random Forest | 0.706 ± 0.107 | 0.711 ± 0.102 | 0.750 ± 0.128 | 0.003 |
| Naive Bayes | 0.689 ± 0.106 | 0.656 ± 0.133 | 0.796 ± 0.112 | 0.002 |
| Gradient Boosting | 0.678 ± 0.124 | 0.688 ± 0.123 | 0.729 ± 0.126 | 0.015 |
| KNN | 0.663 ± 0.117 | 0.661 ± 0.126 | 0.720 ± 0.130 | 0.019 |
| XGBoost | 0.659 ± 0.118 | 0.656 ± 0.123 | 0.725 ± 0.119 | 0.002 |
| Logistic Regression | 0.646 ± 0.110 | 0.642 ± 0.124 | 0.716 ± 0.122 | 0.020 |
| Decision Tree | 0.602 ± 0.107 | 0.592 ± 0.122 | 0.600 ± 0.107 | 0.098^ns^ |
| Neural Network | 0.599 ± 0.131 | 0.599 ± 0.136 | 0.682 ± 0.143 | 0.113^ns^ |

*Values represent mean ± standard deviation from 10 repetitions of 5-fold stratified cross-validation. Permutation p-values derived from 1,000 label permutations. ^ns^Not significant (p ≥ 0.05).*

<mark>

**Figure 2. Machine Learning Model Performance**

![Figure 2](figures/Figure2.png)

**Fig 2.** Machine learning classification of obesity status using 34 differentially abundant metabolites. **(a)** Cross-validation performance comparison for the six top-performing classifiers, showing accuracy (blue), F1-score (green), and AUC (red). Error bars represent standard deviation across 50 cross-validation iterations. SVM achieved the highest accuracy (72.4% ± 9.1%). **(b)** Confusion matrix for SVM on the hold-out test set (n = 22), achieving 63.6% accuracy. **(c)** Permutation test results for seven models that exceeded chance performance (p < 0.05); Decision Tree and Neural Network were excluded as they did not reach statistical significance.

</mark>

<mark>On the independent hold-out test set (n = 22), SVM achieved 63.6% accuracy (Fig 2b). However, the small test set size warrants cautious interpretation; bootstrap 95% confidence intervals were wide (SVM: 45.5–81.8%), reflecting inherent uncertainty. Cross-validation metrics should therefore guide primary interpretation of model performance. Learning curves confirmed appropriate model complexity for SVM and Random Forest, with minimal gaps between training and validation accuracy at full sample size (Supplementary Fig S1a).</mark>

<mark>Feature importance analysis identified L-lysine and quinolinic acid as the most discriminatory metabolites, both elevated in the BMI ≥ 30 kg/m² group (Fig 3a). D-fucose showed the strongest negative association, consistent with its reduced levels in obesity. The six metabolites with the largest effect sizes are displayed with individual participant data points in Fig 3b.</mark>

<mark>

**Figure 3. Metabolite Effect Sizes**

![Figure 3](figures/Figure3.png)

**Fig 3.** Metabolite effect sizes and distributions. **(a)** Cohen's d effect sizes for the top 15 discriminatory metabolites. Dashed lines indicate large effect threshold (|d| = 0.8). **(b)** Box plots with individual data points for the six metabolites with largest effect sizes. Asterisks indicate significance: * p < 0.05, ** p < 0.01.

</mark>

---

## Discussion

This study using ^1^H-NMR metabolomics and machine-learning algorithms identified thirty-four significantly dysregulated metabolites in individuals with obesity and high cardiovascular risk, offering insights into how obesity disrupts key biochemical pathways. The marked increase in waist circumference in the obesity group indicates greater central adiposity, strongly associated with metabolic disturbances like insulin resistance and altered lipid metabolism [28]. This fat accumulation elevates triglyceride levels, a hallmark of atherogenic dyslipidemia driven by increased VLDL production and impaired triglyceride clearance [29]. These combined factors heighten cardiometabolic risk, as shown by significantly higher cardiovascular risk scores in the obesity group, emphasizing the need to address central adiposity and dyslipidemia in preventive healthcare and nutrition strategies [28].

L-lysine emerged as the most significantly elevated metabolite in the BMI >= 30 kg/m^2^ group (Cohen's d = 1.30), followed by adenine, quinolinic acid, and kynurenic acid. Elevated kynurenic acid and quinolinic acid are particularly noteworthy as components of the tryptophan-kynurenine pathway, which has been increasingly linked to obesity-related inflammation and metabolic dysfunction [30,31]. In our study, the BMI >= 30 kg/m^2^ group had a higher percentage of people with high HOMA-IR, indicating a greater risk of insulin resistance, possibly influenced by this metabolomic imbalance.

The elevated amino acid profiles observed (L-lysine, L-alanine, L-leucine) align with the well-documented rise of amino acids in insulin-resistant states [33,34]. This pattern suggests enhanced protein catabolism or altered amino acid uptake in tissues, consistent with metabolic dysfunction in obesity [35]. The significant elevation of benzoic acid and hippuric acid may reflect altered gut microbial metabolism, as these metabolites are largely derived from microbial processing of dietary polyphenols [36,37].

In contrast, D-fucose and myo-inositol were significantly diminished in the obesity group. Myo-inositol depletion is particularly interesting given its role as a second messenger in insulin signaling pathways [38]. Reduced myo-inositol has been associated with insulin resistance and may contribute to the metabolic dysfunction observed in obesity. The pronounced decrease in D-fucose may relate to altered glycoprotein metabolism or gut microbiome changes in obesity [39,40].

Pathway enrichment analysis highlighted phenylalanine, tyrosine, and tryptophan biosynthesis, pyruvate metabolism, and glycolysis/gluconeogenesis among the most significantly altered routes in individuals with BMI >= 30 kg/m^2^. These cascades are integral to energy production, redox balance, and amino-acid interconversion. Evidence suggests that obesity can reprogram these primary pathways, as stated earlier, contributing to insulin resistance, inflammation, and ectopic lipid deposition [6, 33].

<mark>Our protein-metabolite-disease interaction network (Supplementary Fig S2) showed that representative amino acid metabolites from our significant panel map onto obesity-related proteins, underscoring the functional interconnectedness of metabolic reconfigurations.</mark> We noted potential links between these disrupted metabolic profiles and various chronic diseases by incorporating annotations from STITCH, the Human Metabolome Database, and DisGeNET. This aligns with emerging models of obesity that position it not solely as a harbinger of conditions such as type 2 diabetes or cardiovascular disease but as an active contributor to organ dysfunction driven by inflammatory, endocrine, and hemodynamic stressors [40].

### Sample Size Considerations

Our sample size (n=72, 36 per group) warrants discussion. While larger cohorts would increase statistical power, our study demonstrates adequate power for detecting the observed effect sizes. L-lysine, showing the largest effect (Cohen's d = 1.30), achieved >95% power, while metabolites with smaller effects (d ~ 0.8) achieved approximately 80% power. The repeated cross-validation approach, which utilizes all samples for both training and validation across multiple iterations, provides more robust performance estimates than a single train-test split. Importantly, permutation testing confirmed that our classification results significantly exceed chance for most models (seven of nine models p < 0.05, with the best-performing SVM at p = 0.002), indicating genuine discriminatory signal despite the modest sample size. These findings should be interpreted as preliminary evidence requiring validation in larger, independent cohorts.

### Model Selection

Support vector machine (SVM) achieved the highest cross-validation accuracy (72.4%), followed by random forest (70.6%). The choice of SVM as the best-performing model reflects its ability to handle high-dimensional feature spaces effectively when the number of features (34 metabolites) approaches the sample size (n=72). SVM's performance was statistically validated through permutation testing (p < 0.01), confirming that the model captures genuine discriminatory patterns.

Notably, test set performance showed high variability across models (confidence intervals spanning 20-35 percentage points), reflecting the inherent uncertainty with small sample sizes. This underscores the importance of using cross-validation as the primary performance metric and interpreting test set results with appropriate caution.

### Population-Specific Validation

Our findings validate several metabolite associations previously reported in European and North American populations, including altered branched-chain amino acid profiles and dysregulated tryptophan-kynurenine pathway metabolites in obesity. Importantly, this study extends these observations to a Brazilian cohort with distinct dietary patterns and genetic backgrounds. The replication of core metabolic signatures across populations strengthens confidence in these biomarkers' generalizability, while population-specific validation remains essential for clinical translation in diverse settings.

### Clinical Implications

Beyond classification, our machine learning results have clinical translational value, hinting that refined biomarker panels, encompassing amino acids, short-chain organic acids, and possibly lipid intermediates, might surpass BMI alone in detecting early organ stress. The literature increasingly shows that "clinical obesity" should be differentiated from mere increases in body weight or BMI thresholds, focusing instead on organ-level dysfunction [37, 38]. Tracking signature metabolites (*e.g.*, lysine, quinolinic acid, kynurenic acid) could thus facilitate earlier interventions or more precise therapeutic efficacy monitoring.

These findings align with the evolving concept of obesity, moving beyond BMI to include excess adiposity and validated biomarkers of tissue and organ dysfunction for defining "clinical obesity" [40]. Our results demonstrate that specific metabolite changes and pathway disruptions may occur before or alongside clinical symptoms, highlighting the importance of targeted metabolite profiling-especially amino acid and short-chain metabolite panels-for enhanced risk stratification. Furthermore, this approach can guide personalized interventions, such as GLP-1 receptor agonists and lifestyle programs, as well as precision nutrition strategies tailored to individual metabolic phenotypes.

### Limitations

Despite the strengths of our integrative approach, several limitations should be noted. The cross-sectional design limits causal inference, as metabolic changes may be a cause or consequence of obesity. Prospective studies are needed to clarify temporal relationships and predictive metabolites. Semi-quantitative metabolomics provides broad coverage but has quantification limitations; targeted or absolute quantification of key metabolites is recommended for future validation. Additionally, the sample size and demographic characteristics may restrict the generalizability of findings, highlighting the need for larger, multicenter studies across diverse populations. Future studies could determine the thresholds of metabolite dysregulation that most accurately correlate with the onset of organ impairments. Mechanistic and interventional research, including Mendelian randomization, may further clarify whether modifying lysine, kynurenic acid, or quinolinic acid levels influences obesity-related phenotypes or the risk of comorbidities. Finally, integrating multi-omics data (epigenetics, proteomics, microbiomics) with established anthropometric and clinical assessments is likely to refine strategies for precision obesity management, enhancing both early detection of organ dysfunction and targeted therapeutic interventions.

---

## Conclusions

Our investigation of metabolites, pathways, and protein-metabolite-disease networks in individuals with BMI >= 30 kg/m^2^ offers evidence that obesity orchestrates metabolic remodeling affecting amino acids, tryptophan-kynurenine pathway metabolites, and sugar alcohols. L-lysine, quinolinic acid, and kynurenic acid elevations, alongside myo-inositol and D-fucose depletions, highlight disrupted nitrogen metabolism and insulin signaling pathways. The integration of machine learning with rigorous validation-using cross-validation as the primary metric (SVM: 72.4% accuracy) and statistical significance confirmed by permutation testing (p < 0.01)-demonstrates that these metabolite panels contain genuine discriminatory information. While these findings require validation in larger cohorts, they suggest potential for developing metabolite-based risk stratification tools. As obesity definitions become more function-centric, targeted medical and nutrition strategies guided by molecular biomarkers may enhance early detection and personalized interventions.

---

## Data Availability Statement

All data and code necessary to reproduce the analyses are publicly available at https://github.com/omagebright/obesity. This includes: (1) normalized metabolite intensities, (2) sample metadata (de-identified), (3) complete analysis notebooks, and (4) trained model parameters. Raw NMR spectra are available from the corresponding author upon reasonable request, subject to ethical approval constraints.

---

## Funding

<mark>This study was funded in part by the Coordination of Improvement of Higher Education Personnel (Coordenacao de Aperfeicoamento de Pessoal do Nivel Superior-CAPES, Grant number 001). This work was supported by the National Council for Scientific and Technological Development (Conselho Nacional de Desenvolvimento Cientifico e Tecnologico-CNPq, Grant numbers: 431053/2016-2, 405837/2016-0, and 308079/2021-3). We also thank INCTBio-Lauro Kubota and Sao Paulo Research Foundation (Fundacao de Amparo a Pesquisa do Estado de Sao Paulo-FAPESP), Grant numbers: #2025/23708-6, #2023/02691-2, #2022/11207-4, #2018/24069-3, #2016/20054-6, and #2014/50867-3.</mark>

## Conflicts of Interest

The authors declare no conflict of interest or any competing financial interests in relation to the work described.

## Author Contributions

S.C.V.C.L, C.O.L, D.M.L.M, K.C.M.S.E. investigation and performed clinical research, and collected data. P.E.N.R.B and A.C.S.J. organization and data analysis. E.S.B and L.G.M. performed NMR measurements and metabolomic analysis. F.B.O performed machine learning analysis and prepared the figures, data interpretation. K.C.M.S.E, A.C.S.J., L.T., and F.B.O. wrote the original draft of the manuscript. L.T. data interpretation, manuscript correction, and supervision of the study. S.C.V.C.L, C.O.L, D.M.L.M, L.F.C., F.B.J., and L.T. revised the manuscript and contributed to discussions. All authors read and approved the final version of the manuscript.

---

## Supplementary Materials

### Supplementary Table S2: Model Hyperparameters

| Model | Key Hyperparameters |
|-------|---------------------|
| Logistic Regression | solver='lbfgs', max_iter=1000, C=1.0 |
| Random Forest | n_estimators=100, max_depth=None, min_samples_split=2 |
| XGBoost | n_estimators=100, max_depth=6, learning_rate=0.1 |
| SVM | kernel='rbf', C=1.0, gamma='scale' |
| KNN | n_neighbors=5, weights='uniform' |
| Decision Tree | max_depth=None, min_samples_split=2 |
| Naive Bayes | Default (no hyperparameters) |
| Neural Network | hidden_layer_sizes=(100,), max_iter=1000 |
| Gradient Boosting | n_estimators=100, max_depth=3, learning_rate=0.1 |

*All models used random_state=42 for reproducibility.*

### Supplementary Table S3: Complete List of Significant Metabolites (Training Set, FDR < 0.05)

| Metabolite | Mean BMI>=30 | Mean BMI<30 | log2FC | Cohen's d | p_adj |
|------------|--------------|-------------|--------|-----------|-------|
| L-Lysine | 0.00503 | 0.00365 | 0.46 | 1.30 | 0.004 |
| Adenine | 0.00103 | 0.00056 | 0.89 | 1.17 | 0.010 |
| Quinolinic Acid | 0.00260 | 0.00165 | 0.65 | 1.12 | 0.010 |
| D-Fucose | 0.00168 | 0.00657 | -1.97 | -1.11 | 0.010 |
| Kynurenic Acid | 0.00171 | 0.00110 | 0.64 | 1.07 | 0.010 |
| 2-Picolinic Acid | 0.00219 | 0.00111 | 0.98 | 1.06 | 0.010 |
| L-Alanine | 0.00528 | 0.00448 | 0.24 | 1.07 | 0.010 |
| Benzoic Acid | 0.00117 | 0.00070 | 0.73 | 1.03 | 0.011 |
| Hypoxanthine | 0.00059 | 0.00031 | 0.92 | 1.02 | 0.011 |
| Phenylglyoxylic Acid | 0.00122 | 0.00077 | 0.66 | 1.04 | 0.011 |
| Nicotinic Acid | 0.00152 | 0.00085 | 0.83 | 0.99 | 0.013 |
| Indoxylsulfate | 0.00194 | 0.00125 | 0.63 | 1.00 | 0.013 |
| Oxypurinol | 0.00060 | 0.00018 | 1.76 | 0.97 | 0.015 |
| Mandelic Acid | 0.00049 | 0.00028 | 0.82 | 0.95 | 0.015 |
| 2-Hydroxybutyric Acid | 0.00501 | 0.00405 | 0.31 | 0.94 | 0.015 |
| Ethylmalonic Acid | 0.01221 | 0.01123 | 0.12 | 0.95 | 0.016 |
| Pyrocatechol | 0.00111 | 0.00049 | 1.18 | 0.92 | 0.017 |
| Sebacic Acid | 0.00456 | 0.00543 | -0.25 | -0.92 | 0.017 |
| 2-Methylglutaric Acid | 0.00291 | 0.00245 | 0.25 | 0.92 | 0.017 |
| L-Aspartate | 0.00322 | 0.00549 | -0.77 | -0.92 | 0.017 |
| L-Glutamic Acid | 0.00896 | 0.00690 | 0.38 | 0.90 | 0.018 |
| D-Glucose-6-Phosphate | 0.00076 | 0.00169 | -1.15 | -0.88 | 0.020 |
| Azelaic Acid | 0.00138 | 0.00076 | 0.87 | 0.88 | 0.020 |
| Hippuric Acid | 0.00118 | 0.00096 | 0.29 | 0.88 | 0.020 |
| Xylitol | 0.00125 | 0.00059 | 1.08 | 0.86 | 0.021 |
| L-Asparagine | 0.00309 | 0.00490 | -0.66 | -0.84 | 0.029 |
| 2-Oxoisovalerate | 0.00810 | 0.00724 | 0.16 | 0.81 | 0.032 |
| Cadaverine | 0.00081 | 0.00129 | -0.67 | -0.81 | 0.032 |
| 3-Hydroxybutyrate | 0.00699 | 0.00622 | 0.17 | 0.81 | 0.033 |
| L-Leucine | 0.00453 | 0.00398 | 0.19 | 0.80 | 0.033 |
| Myo-Inositol | 0.00164 | 0.00286 | -0.80 | -0.79 | 0.035 |
| L-Isoleucine | 0.00795 | 0.00717 | 0.15 | 0.79 | 0.037 |
| 4-AminoHippuric Acid | 0.00087 | 0.00063 | 0.46 | 0.77 | 0.039 |
| Propionate | 0.00058 | 0.00020 | 1.51 | 0.77 | 0.040 |

*Values represent relative ASICS intensities from training set only (n=50). p_adj, p-value adjusted for FDR using Benjamini-Hochberg method. Positive log2FC indicates higher levels in BMI >= 30 group.*

<mark>

### Supplementary Figure S1: Model Validation

![Supplementary Figure S1](figures/FigureS1_Validation.png)

**Fig S1.** Model validation analysis. **(a)** Learning curves for six classification models showing training and validation accuracy as a function of training set size. Shaded regions indicate +/-1 standard deviation across cross-validation folds. SVM shows the best generalization with minimal overfitting. **(b)** Permutation test distributions for top three models, demonstrating actual model performance (red line) significantly exceeds the null distribution.

</mark>

### Supplementary Figure S2: Protein-Metabolite-Disease Network

Protein-metabolite-disease interaction network showing functional relationships among significantly altered metabolites. Nodes represent metabolites, proteins, and diseases; edges indicate known functional or literature-based relationships from STITCH database (confidence >= 0.70).

---

## References

1. World Health Organization (WHO). Prevalence of obesity among adults, BMI >= 30 (age-standardized estimate) (%). Available at: https://www.who.int/data/gho/data/indicators/indicator-details/GHO/prevalence-of-obesity-among-adults-bmi--30-(age-standardized-estimate)-(-). Accessed May 4, 2024.

2. Endalifer ML, Diress G. Epidemiology, Predisposing Factors, Biomarkers, and Prevention Mechanism of Obesity: A Systematic Review. *Journal of Obesity* 2024;31:6134362.

3. Henning RJ. Obesity and Obesity-Induced Inflammatory Disease Contribute to Atherosclerosis: A Review of the Pathophysiology and Treatment of Obesity. *American Journal of Cardiovascular Disease* 2021;11:504-529.

4. Safaei M, Sundararajan EA, Driss M, Boulila W, Shapi'i A. A Systematic Literature Review on Obesity: Understanding the Causes and Consequences of Obesity and Reviewing Various Machine Learning Approaches Used to Predict Obesity. *Computers in Biology and Medicine* 2021;136:104754.

5. Ostrominski JW, Powell-Wiley TM. Risk Stratification and Treatment of Obesity for Primary and Secondary Prevention of Cardiovascular Disease. *Current Atherosclerosis Reports* 2024;26:11-23.

6. Rubino F, Cummings DE, Eckel RH, Cohen RV, Wilding JPH, Brown WA, et al. Definition and Diagnostic Criteria of Clinical Obesity. *The Lancet Diabetes & Endocrinology* 2025;13:221-262.

7. Nowak MM, Niemczyk M, Gołębiewski S, Pączek L. Impact of Body Mass Index on All-Cause Mortality in Adults: A Systematic Review and Meta-Analysis. *Journal of Clinical Medicine* 2024;13:2305.

8. Bosy-Westphal A, Müller MJ. Diagnosis of Obesity Based on Body Composition-Associated Health Risks: Time for a Change in Paradigm. *Obesity Reviews* 2021;22:e13190.

9. Regan JA, Shah SH. Obesity Genomics and Metabolomics: A Nexus of Cardiometabolic Risk. *Current Cardiology Reports* 2020;22:174.

10. Gonzalez-Covarrubias V, Martínez-Martínez E, Del Bosque-Plata L. The Potential of Metabolomics in Biomedical Applications. *Metabolites* 2022;12:194.

11. Rangel-Huerta OD, Pastor-Villaescusa B, Gil A. Are We Close to Defining a Metabolomic Signature of Human Obesity? A Systematic Review of Metabolomics Studies. *Metabolomics* 2019;15:93.

12. Bellot PENR, Moia MN, Reis BZ, et al. Are Phosphatidylcholine and Lysophosphatidylcholine Body Levels Potentially Reliable Biomarkers in Obesity? A Review of Human Studies. *Molecular Nutrition & Food Research* 2023;67:e2200568.

13. Bellot PENR, Braga ES, Omage FB, et al. Plasma Lipid Metabolites as Potential Biomarkers for Identifying Individuals at Risk of Obesity-Induced Metabolic Complications. *Scientific Reports* 2023;13:11729.

14. Sociedade Brasileira de Cardiologia. Atualização da Diretriz Brasileira de Dislipidemias e Prevenção da Aterosclerose – 2017. *Arquivos Brasileiros de Cardiologia* 2017;109:1-30.

15. Barroso WKS, Rodrigues CIS, Bortolotto LA, et al. Diretrizes Brasileiras de Hipertensão Arterial — 2020. *Arquivos Brasileiros de Cardiologia* 2021;116:516-658.

16. World Health Organization. Obesity: preventing and managing the global epidemic. Report of a WHO consultation. *World Health Organization Technical Report Series* 2000;894:253.

17. Amato MC, Giordano C. Visceral Adiposity Index: An Indicator of Adipose Tissue Dysfunction. *International Journal of Endocrinology* 2014;730827.

18. Rodacki M, Cobas RA, Zajdenverg L, Silva Júnior WS, Giacaglia L, Calliari LE, et al. Diagnóstico de Diabetes Mellitus. *Diretriz Oficial da Sociedade Brasileira de Diabetes* 2024. ISBN: 978-65-272-0704-7.

19. Ridker PM. C-Reactive Protein and the Prediction of Cardiovascular Events Among Those at Intermediate Risk: Moving an Inflammatory Hypothesis Toward Consensus. *Journal of the American College of Cardiology* 2007;49:2129-2138.

20. Geloneze B, Vasques AC, Stabe CF, Pareja JC, Rosado LE, Queiroz EC, et al. HOMA1-IR and HOMA2-IR Indexes in Identifying Insulin Resistance and Metabolic Syndrome: Brazilian Metabolic Syndrome Study (BRAMS). *Arquivos Brasileiros de Endocrinologia & Metabologia* 2009;53:281-287.

21. Ulrich EL, Akutsu H, Doreleijers JF, Harano Y, Ioannidis YE, Lin J, et al. BioMagResBank. *Nucleic Acids Research* 2008;36:D402-D408.

22. Wishart DS, Knox C, Guo AC, Eisner R, Young N, Gautam B, et al. HMDB: A Knowledgebase for the Human Metabolome. *Nucleic Acids Research* 2009;37:D603-D610.

23. Cloarec O, Dumas ME, Craig A, Barton RH, Trygg J, Hudson J, et al. Statistical Total Correlation Spectroscopy: An Exploratory Approach for Latent Biomarker Identification from Metabolic 1H NMR Data Sets. *Analytical Chemistry* 2005;77:1282-1289.

24. Pang Z, Lu Y, Zhou G, Hui F, Xu L, Viau C, Spigelman AF, MacDonald PE, Wishart DS, Li S, Xia J. MetaboAnalyst 6.0: towards a unified platform for metabolomics data processing, analysis and interpretation. *Nucleic Acids Research* 2024;52(W1):W398-W406.

25. Szklarczyk D, Santos A, von Mering C, Jensen LJ, Bork P, Kuhn M, et al. STITCH 5: Augmenting Protein-Chemical Interaction Networks with Tissue and Affinity Data. *Nucleic Acids Research* 2016;44:D380-D384.

26. Piñero J, Bravo À, Queralt-Rosinach N, Gutiérrez-Sacristán A, Deu-Pons J, Centeno E, et al. DisGeNET: A Comprehensive Platform Integrating Information on Human Disease-Associated Genes and Variants. *Nucleic Acids Research* 2017;45:D833-D839.

27. Shannon P, Markiel A, Ozier O, Baliga NS, Wang JT, Ramage D, et al. Cytoscape: A Software Environment for Integrated Models of Biomolecular Interaction Networks. *Genome Research* 2003;13:2498-2504.

28. Powell-Wiley TM, Poirier P, Burke LE, Després JP, Gordon-Larsen P, Lavie CJ, et al. Obesity and Cardiovascular Disease: A Scientific Statement From the American Heart Association. *Circulation* 2021;143:e984-e1010.

29. Vekic J, Stefanovic A, Zeljkovic A. Obesity and Dyslipidemia: A Review of Current Evidence. *Current Obesity Reports* 2023;12:207-222.

30. Bartoloni B, Mannelli M, Gamberi T, Fiaschi T. The Multiple Roles of Lactate in the Skeletal Muscle. *Cells* 2024;13:1177.

31. Lin Y, Bai M, Wang S, Chen L, Li Z, Li C, et al. Lactate Is a Key Mediator That Links Obesity to Insulin Resistance via Modulating Cytokine Production from Adipose Tissue. *Diabetes* 2022;71:637-652.

32. Tang Q, Tan P, Ma N, Ma X. Physiological Functions of Threonine in Animals: Beyond Nutrition Metabolism. *Nutrients* 2021;13:2592.

33. Newgard CB, An J, Bain JR, Muehlbauer MJ, Stevens RD, Lien LF, et al. A Branched-Chain Amino Acid–Related Metabolic Signature That Differentiates Obese and Lean Humans and Contributes to Insulin Resistance. *Cell Metabolism* 2009;9:311-326.

34. Vanweert F, Schrauwen P, Phielix E. Role of Branched-Chain Amino Acid Metabolism in the Pathogenesis of Obesity and Type 2 Diabetes-Related Metabolic Disturbances BCAA Metabolism in Type 2 Diabetes. *Nutrition & Diabetes* 2022;12:35.

35. Lynch CJ, Adams SH. Branched-Chain Amino Acids in Metabolic Signalling and Insulin Resistance. *Nature Reviews Endocrinology* 2014;10:723-736.

36. Neinast MD, Jang C, Hui S, Murashige DS, Chu Q, Morscher RJ, et al. Quantitative Analysis of the Whole-Body Metabolic Fate of Branched-Chain Amino Acids. *Cell Metabolism* 2019;30:1141-1156.e7.

37. Yamakado M, Tanaka T, Nagao K, Imaizumi A, Komatsu M, Daimon T, et al. Plasma Amino Acid Profile Associated with Fatty Liver Disease and Co-Occurrence of Metabolic Risk Factors. *Scientific Reports* 2017;7:14485.

38. Hoffer LJ, Forse RA. Protein Metabolic Effects of a Prolonged Fast and Hypocaloric Refeeding. *American Journal of Physiology* 1990;258:E832-E840.

39. Franckhauser S, Elias I, Rotter Sopasakis V, Ferré T, Nagaev I, Andersson CX, et al. Overexpression of IL-6 Leads to Hyperinsulinaemia, Liver Inflammation and Reduced Body Weight in Mice. *Diabetologia* 2008;51:1306-1316.

40. Dumas M-E, Barton RH, Toye A, Cloarec O, Blancher C, Rothwell A, et al. Metabolic Profiling Reveals a Contribution of Gut Microbiota to Fatty Liver Phenotype in Insulin-Resistant Mice. *Proceedings of the National Academy of Sciences of the United States of America* 2006;103:12511-12516.
