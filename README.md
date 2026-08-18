# Hermes-Multiomics-Platform
# Multi-omic Analysis of Malignant Pleural Mesothelioma
## MEDUSA Multi-omic Analysis and Hermes AI Agent Evaluation
### University of Leicester — MSc Bioinformatics Independent Research Project

---

## Overview

This repository contains the **Python and R analysis code** used for my MSc Bioinformatics Independent Research Project at the University of Leicester.

The project investigated malignant pleural mesothelioma using integrated multi-omic data from the **MEDUSA cohort**, with a particular focus on the molecular and immune features associated with **MTAP (9p21) loss**.

The project also evaluated **Hermes**, a locally hosted large language model (LLM) agent, as a tool for interrogating patient-level multi-omic data. Hermes-assisted analyses were independently checked using conventional Python and R workflows to assess the accuracy, reproducibility and limitations of agent-assisted biomedical research.

The main analyses included:

- Patient and clinicopathological characterisation
- Multi-omic data availability and feature landscape analysis
- Epigenetic age acceleration analysis
- MTAP-associated genomic and transcriptional feature discovery
- Machine-learning classification of MTAP status
- Immune microenvironment analysis
- Integration of epigenetic, genomic and immune features
- Evolutionary architecture analysis
- Independent validation of Hermes outputs
- Candidate druggable target synthesis

---

## Data

The primary analysis was based on the integrated MEDUSA dataset:

`MEDUSA_Master_ULTIMATE.csv`

The integrated dataset contained:

- **235 patients**
- **1,426 variables**
- Clinical and pathological information
- Genomic features
- Transcriptomic features
- Methylation-derived epigenetic features
- Immune and cell-type deconvolution features

A subset of **97 patients** had complete clinical and pathological information available for clinicopathological characterisation.

Additional CONFIRM data were used for selected immune-marker analyses and for benchmarking Hermes analytical performance.

### Data Availability

The patient-level MEDUSA and CONFIRM datasets are **not included in this repository** because they are controlled research datasets subject to institutional data-governance and ethical requirements.

The scripts therefore assume authorised local access to the corresponding datasets.

No patient-level research data are distributed through this repository.

---

## Scripts

### Python notebooks

| Script | Description |
|--------|-------------|
| `Analysing_Master(1).ipynb` | Clinical and patient-characteristic analysis of the integrated MEDUSA dataset |
| `MEDUSA_barchart(1).ipynb` | Analysis and visualisation of data availability and missingness across MEDUSA modalities |
| `Feature_Landscape(1).ipynb` | Audit and visualisation of the 1,426-feature MEDUSA master dataset |
| `FDR(1).ipynb` | Epigenetic clock statistical analysis and multiple-testing correction |
| `methylclock_hist(1).ipynb` | Methylation clock visualisation and associated analyses |
| `Hallmark_PCA(1).ipynb` | MTAP-associated hallmark pathway analysis and principal component analysis |
| `agreed_features(1).ipynb` | Analysis and visualisation of MTAP-associated significant features |
| `RF,GB,ROC,AUC(1).ipynb` | Random Forest and Gradient Boosting classification, ROC curves and AUC analysis |
| `Figure_Dendritic(1).ipynb` | Activated dendritic-cell analysis and visualisation |
| `immune_hist.ipynb` | Integrated methylation and immune-feature analyses |
| `evol_arch_mtap(1).ipynb` | Clonal and subclonal evolutionary architecture analysis by MTAP status |
| `Hermes_Independent_Val(1).ipynb` | Independent validation of Hermes-generated analytical outputs |
| `table_druggable_targets(1).ipynb` | Generation of the candidate druggable-target summary |
| `Tables_Workflow(1).ipynb` | Generation of workflow figures and supporting tables |

### R scripts

| Script | Description |
|--------|-------------|
| `MTAP_status_VPlots(2).R` | R-based MTAP violin plots and validation of immune, HRD and transcription-factor features |
| `pheno_mean_markers(2).R` | Independent CONFIRM immune-marker summary statistics and boxplot |
| `epig_clock_MTAP(2).R` | Independent R analysis of methylation clocks, genomic instability and stratified survival |

---

## Hermes AI Agent

Hermes was used as an agentic interface for interrogation of the integrated multi-omic dataset.

Its performance was evaluated using a three-stage framework:

1. **Stage 1 — Data retrieval accuracy**  
   Evaluation of Hermes's ability to retrieve factual information from the available MEDUSA data.

2. **Stage 2 — Analytical capability**  
   Assessment of Hermes's ability to perform quantitative analysis and generate visualisations, including benchmarking against independently generated Python and R outputs.

3. **Stage 3 — Autonomous hypothesis generation**  
   Evaluation of Hermes's ability to identify MTAP-associated findings and synthesise candidate biological and therapeutic hypotheses without step-by-step analytical instructions.

All Hermes outputs used to support reported findings were independently reproduced using Python or R.

The evaluation identified both correctable analytical errors and cases in which previously incorrect findings could reappear in later Hermes outputs, highlighting the importance of independent verification when using LLM agents for biomedical data analysis.

---

## Requirements

### Python

Python 3.12 was used for the primary statistical and machine-learning analyses.

Main Python packages:

- pandas
- NumPy
- SciPy
- statsmodels
- scikit-learn
- lifelines
- Matplotlib
- seaborn
- adjustText

Install the required Python packages with:

```bash
pip install pandas numpy scipy statsmodels scikit-learn lifelines matplotlib seaborn adjustText
