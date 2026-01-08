# Plasma 1H-NMR Metabolomics for Obesity Classification

[![DOI](https://img.shields.io/badge/DOI-pending-blue)]()
[![License: CC BY 4.0](https://img.shields.io/badge/License-CC%20BY%204.0-lightgrey.svg)](https://creativecommons.org/licenses/by/4.0/)

## Overview

This repository contains the data and code for the study:

**"Plasma ¹H-NMR metabolomics identifies amino acid and carbohydrate pathway alterations associated with obesity"**

Sena-Evangelista KCM, Omage FB, Braga ES, et al.

## Study Summary

- **Objective:** Identify plasma metabolites that distinguish adults with BMI ≥ 30 kg/m² from those with BMI < 30 kg/m²
- **Methods:** Untargeted ¹H-NMR metabolomics on plasma samples from 72 adults (36 per group)
- **Key Findings:** Seven metabolites significantly differed between groups; Logistic Regression achieved 81.8% accuracy (F1 = 0.83) for obesity classification

## Repository Structure

```
├── data/
│   ├── raw/                    # Raw metabolite intensities
│   ├── processed/              # Normalized and scaled data
│   └── metadata/               # Sample characteristics (de-identified)
├── notebooks/
│   ├── 03_ml_classification.ipynb      # Machine learning models
│   └── 04_enhanced_validation.ipynb    # Enhanced validation analyses
├── requirements.txt            # Python dependencies
├── CITATION.cff                # Citation information
└── README.md
```

## Key Results

### Significant Metabolites (FDR < 0.05)

| Metabolite | Direction in BMI ≥ 30 | log₂FC | P_adj |
|------------|----------------------|--------|-------|
| Formate | ↓ Down | -0.585 | 0.0021 |
| Phenylalanine | ↓ Down | -0.156 | 0.0484 |
| Leucine | ↓ Down | -0.116 | 0.0484 |
| Alanine | ↓ Down | -0.123 | 0.0484 |
| Lysine | ↓ Down | -0.108 | 0.0484 |
| Lactate | ↑ Up | +0.160 | 0.0484 |
| Threonine | ↑ Up | +0.157 | 0.0484 |

### Machine Learning Performance

**Best Model: Logistic Regression**

| Metric | Test Set (n=22) | 95% CI |
|--------|-----------------|--------|
| Accuracy | 0.818 | 0.636 - 0.955 |
| F1 Score | 0.833 | 0.632 - 0.960 |
| AUC | 0.826 | 0.617 - 0.975 |
| MCC | 0.647 | 0.297 - 0.912 |

**Cross-Validation (10×5-fold):** Accuracy = 0.723 ± 0.105

**Permutation Test:** p = 0.005 (model performs significantly better than chance)

## Reproducibility

### Requirements

```bash
pip install -r requirements.txt
```

### Running the Analysis

1. Clone this repository
2. Install dependencies: `pip install -r requirements.txt`
3. Run notebooks in order (01 → 04)

### Random Seeds

All analyses use `random_state=42` for reproducibility.

## Data Availability

All data necessary to reproduce the results are included in this repository. Raw NMR spectra are available upon reasonable request.

## Citation

If you use this data or code, please cite:

```
Sena-Evangelista KCM, Omage FB, Braga ES, Martins LG, Bellot PENR,
de Souza Júnior AC, Pedrosa LFC, Marchioni DML, Barbosa Jr F,
Lima SCVC, Lyra CO, Tasic L. Plasma ¹H-NMR metabolomics identifies
amino acid and carbohydrate pathway alterations associated with obesity.
[Journal Name]. 2026.
```

## License

This work is licensed under a [Creative Commons Attribution 4.0 International License](https://creativecommons.org/licenses/by/4.0/).

## Contact

- **Corresponding Author:** Karine Cavalcanti M. Sena-Evangelista (karine.sena@ufrn.br)
- **Data/Code Questions:** Folorunsho Bright Omage (bright@unicamp.br)

## Acknowledgments

This study was funded by the Coordination of Improvement of Higher Education Personnel (CAPES, Grant number 001). This work was supported by the National Council for Scientific and Technological Development (CNPq, Grant numbers: 431053/2016-2, 405837/2016-0, and 308079/2021-3). We also thank INCTBio-Lauro Kubota and Sao Paulo Research Foundation (FAPESP), Grant numbers: #2023/02691-2, #2022/11207-4, #2018/24069-3, #2016/20054-6, and #2014/50867-3. We thank the participants of the BRAZUCA study.
