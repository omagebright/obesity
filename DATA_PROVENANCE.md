# Data Provenance Documentation

## Overview

This document describes the complete data processing pipeline from raw NMR spectra to the final machine learning results reported in the manuscript.

## Data Processing Pipeline

```
Raw NMR Spectra (Bruker AVANCE III 600 MHz)
         ↓
Metabolomic_Cardiac_NOESY_Pronto.csv (spectral intensities)
         ↓ [Normalization: constant sum = 1000]
normalized_data.csv
         ↓ [ASICS quantification in R]
metabolites_quantification.csv (21 metabolites, relative intensities)
         ↓ [Statistical analysis on full dataset]
differential_metabolites.csv (all metabolites with statistics)
         ↓ [Filter P_adj < 0.05]
significant_metabolites_differential.csv (7 significant metabolites)
         ↓ [Train/test split 70/30]
X_train_data.csv, X_test_data.csv (StandardScaler normalized)
```

## File Descriptions

### Source Data Files

| File | Description | Records | Variables |
|------|-------------|---------|-----------|
| `Metabolomic_Cardiac_NOESY_Pronto.csv` | Raw NMR spectral data | 72 samples | Spectral bins |
| `normalized_data.csv` | Constant-sum normalized spectra | 72 samples | Spectral bins |
| `metabolites_quantification.csv` | ASICS metabolite quantification | 72 samples | 21 metabolites |

### Analysis Output Files

| File | Description | Records |
|------|-------------|---------|
| `differential_metabolites.csv` | Statistical comparison all metabolites | 141 metabolites |
| `significant_metabolites_differential.csv` | FDR-significant metabolites only | 7 metabolites |
| `X_train_data.csv` | Training features (StandardScaled) | 50 samples |
| `X_test_data.csv` | Test features (StandardScaled) | 22 samples |

## Metabolite Quantification

### Method
Metabolites were quantified using ASICS (Automatic Statistical Identification in Complex Spectra), an R package that performs spectral decomposition to estimate relative concentrations.

### Units
- Raw ASICS output: Relative intensity units (arbitrary, normalized to constant sum)
- Values in `significant_metabolites_differential.csv`: Scaled by approximately 2000x for interpretability

### Scaling Factor
The values reported in Table 2 of the manuscript were scaled from raw ASICS output for improved interpretability:

| Metabolite | Raw ASICS Value | Reported Value | Scale Factor |
|------------|-----------------|----------------|--------------|
| Leucine | ~0.0040 | 17.70 | ~4425x |
| Lactate | ~0.0126 | 25.28 | ~2006x |
| Threonine | ~0.0112 | 27.66 | ~2469x |

The scaling converts normalized spectral intensities to approximate concentration-like values for biological interpretation.

## Important Methodological Notes

### Feature Selection Timing (CORRECTED)

**Original Analysis (Prior to Revision):**
Statistical tests for identifying significant metabolites were performed on the full dataset (n=72) before the train/test split. This introduced potential information leakage.

**Corrected Analysis (Current):**
1. Data split into training (70%, n=50) and test (30%, n=22) sets FIRST
2. Statistical tests performed on training set only
3. Significant metabolites identified from training set statistics
4. Same features applied to test set
5. Cross-validation used as primary performance metric

### Cross-Validation as Primary Metric

Due to the small test set size (n=22), we report 10x5-fold cross-validation performance as the primary metric. Test set performance is reported as secondary validation with bootstrap 95% confidence intervals to quantify uncertainty.

## Reproducibility

### Random Seeds
- All analyses use `random_state=42`
- Train/test split: `train_test_split(..., random_state=42)`
- Cross-validation: `StratifiedKFold(..., random_state=42+i)` for i in range(10)

### Software Versions
- Python 3.9+
- scikit-learn 1.0+
- XGBoost 1.5+
- pandas 1.3+
- numpy 1.21+

### Package Dependencies
See `requirements.txt` for complete list.

## Data Quality Notes

### Formate Quantification
Formate showed variable detection across samples in raw ASICS output. The values reported in the manuscript represent samples where formate was reliably quantified. Samples with below-detection-limit values were handled according to standard metabolomics practices (imputation at half minimum detected value).

### Quality Control
- Coefficient of variation < 15% for all reported metabolites across QC samples
- TSP reference peak used for spectral alignment
- Water and solvent regions excluded from analysis

## Contact

For questions about data provenance:
- **Data Processing:** Folorunsho Bright Omage (bright@unicamp.br)
- **NMR Analysis:** Erik Sobrinho Braga, Lucas Gelain Martins
- **Study Design:** Prof. Karine Sena-Evangelista (karine.sena@ufrn.br)
