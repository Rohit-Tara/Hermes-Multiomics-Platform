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

Hermes was used as an agentic interface for interrogation of the integrated multi-omic dataset. Unlike the independently written Python and R analyses contained in this repository, Hermes analyses were initiated through natural-language prompts, with the agent autonomously generating the code and outputs required to respond.

Its performance was evaluated using a three-stage framework:

1. **Stage 1 — Data retrieval accuracy**  
   Evaluation of Hermes's ability to retrieve factual information from the available MEDUSA data.

2. **Stage 2 — Analytical capability**  
   Assessment of Hermes's ability to perform quantitative analysis and generate visualisations, including benchmarking against independently generated Python and R outputs.

3. **Stage 3 — Autonomous hypothesis generation**  
   Evaluation of Hermes's ability to identify MTAP-associated findings and synthesise candidate biological and therapeutic hypotheses without step-by-step analytical instructions.

Hermes-generated findings used to support the reported results were independently reproduced or verified using Python or R. The evaluation identified both correctable analytical errors and instances in which previously incorrect findings could reappear in later Hermes outputs, demonstrating the importance of independent verification when using LLM agents for biomedical data analysis.

---

## Requirements

The analyses in this repository were performed using **Python 3.12.12** and **R 4.6.1**.

### Python Requirements

Python was used for the primary data processing, statistical analysis, machine-learning analysis and figure generation.

The Python environment used for this project contained:

| Package | Version |
|---------|---------|
| pandas | 3.0.0 |
| NumPy | 2.4.1 |
| SciPy | 1.17.0 |
| statsmodels | 0.14.6 |
| scikit-learn | 1.8.0 |
| lifelines | 0.30.3 |
| Matplotlib | 3.10.8 |
| seaborn | 0.13.2 |
| adjustText | 1.4.0 |

The required Python packages can be installed using:

```bash
pip install pandas==3.0.0 numpy==2.4.1 scipy==1.17.0 statsmodels==0.14.6 scikit-learn==1.8.0 lifelines==0.30.3 matplotlib==3.10.8 seaborn==0.13.2 adjustText==1.4.0
```

### R Requirements

R 4.6.1 was used for selected independent statistical analyses, validation and figure generation.

The R packages used were:

- dplyr
- readr
- tidyr
- stringr
- ggplot2
- survival

The required R packages can be installed from within R using:

```r
install.packages(c(
    "dplyr",
    "readr",
    "tidyr",
    "stringr",
    "ggplot2",
    "survival"
))
```

---

## Software Versions

| Software / Package | Version |
|--------------------|---------|
| Python | 3.12.12 |
| pandas | 3.0.0 |
| NumPy | 2.4.1 |
| SciPy | 1.17.0 |
| statsmodels | 0.14.6 |
| scikit-learn | 1.8.0 |
| lifelines | 0.30.3 |
| Matplotlib | 3.10.8 |
| seaborn | 0.13.2 |
| adjustText | 1.4.0 |
| R | 4.6.1 |
| dplyr | 1.2.1 |
| readr | 2.2.0 |
| tidyr | 1.3.2 |
| ggplot2 | Version used in the R environment |
| survival | 3.8-9 |
| Hermes | Local deployment |

---

## Usage

1. Clone or download this repository.

2. Install the required Python packages:

```bash
pip install pandas==3.0.0 numpy==2.4.1 scipy==1.17.0 statsmodels==0.14.6 scikit-learn==1.8.0 lifelines==0.30.3 matplotlib==3.10.8 seaborn==0.13.2 adjustText==1.4.0
```

3. Install the required R packages:

```r
install.packages(c(
    "dplyr",
    "readr",
    "tidyr",
    "stringr",
    "ggplot2",
    "survival"
))
```

4. Obtain authorised access to the required MEDUSA and/or CONFIRM datasets.

5. Place the required data files in an appropriate local directory.

6. Update the local file paths within the Python notebooks and R scripts where required.

7. Run the relevant Python notebook according to the analysis being reproduced.

8. Run the R scripts where independent analysis or figure reproduction is required.

The scripts generate the statistical analyses, tables and figures used to support the findings reported in the dissertation.

Because the underlying patient-level datasets cannot be distributed publicly, the complete analyses cannot be reproduced from this repository without authorised access to the corresponding research data.

---

## Repository Structure

```text
Hermes-Multiomics-Platform/
│
├── README.md
│
├── Python/
│   ├── agreed_features(1).ipynb
│   ├── Analysing_Master(1).ipynb
│   ├── evol_arch_mtap(1).ipynb
│   ├── FDR(1).ipynb
│   ├── Feature_Landscape(1).ipynb
│   ├── Figure_Dendritic(1).ipynb
│   ├── Hallmark_PCA(1).ipynb
│   ├── Hermes_Independent_Val(1).ipynb
│   ├── immune_hist.ipynb
│   ├── MEDUSA_barchart(1).ipynb
│   ├── methylclock_hist(1).ipynb
│   ├── RF,GB,ROC,AUC(1).ipynb
│   ├── table_druggable_targets(1).ipynb
│   └── Tables_Workflow(1).ipynb
│
└── R/
    ├── MTAP_status_VPlots(2).R
    ├── pheno_mean_markers(2).R
    └── epig_clock_MTAP(2).R
```
### Hermes-Generated Analyses and Code Availability

Hermes was operated through natural-language prompts and autonomously determined analytical steps where required, including file selection, data parsing, statistical analysis, cross-file integration and generation of outputs and visualisations.

Code and analytical outputs generated autonomously by Hermes are not included in this repository as researcher-authored source code. The Python notebooks and R scripts provided here instead represent the independently implemented analyses used throughout the project, including conventional analyses used to reproduce, verify and evaluate Hermes-generated findings.

This distinction reflects the methodology used in the dissertation. Hermes-generated outputs were not treated as ground truth and were independently checked against the underlying MEDUSA data using Python or R before being retained for downstream interpretation. Outputs that could not be independently reproduced or whose provenance could not be established were excluded from downstream use.

The absence of Hermes-generated scripts from this repository therefore reflects their provenance as autonomous AI-generated analytical code rather than researcher-authored code. The dissertation documents the natural-language prompting framework, validation procedure, performance assessment and identified failure modes of the Hermes agent.
