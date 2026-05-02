
# Replication Package: A Multiway Cluster-Robust Test for Stochastic Dominance

**Authors:** Yaqian Wu, Qingliang Fan, and Yu-Chin Hsu

This repository contains all the necessary R code to reproduce the simulation results and empirical findings presented in the manuscript **"A Multiway Cluster-Robust Test for Stochastic Dominance"**.

---

## 1. System Requirements

The code is implemented in **R**. Most scripts utilize parallel computing and will automatically detect available CPU cores to optimize runtime.

### Required R Packages
You can install all dependencies by running the following command in your R console:

```r
install.packages(c("data.table", "foreach", "parallel", "doParallel", 
                   "doSNOW", "tidyverse", "stringr"))
```

*The code was developed and tested using R version 4.0.0 or higher.*

---

## 2. Directory Structure

The `replication.zip` package is organized as follows:

```text
/R
  └── myFUN.R               # Centralized utility functions (CDFs, Bootstrap functions)

/Simulations
  ├── sim_datagen.R         # Data-generating processes (DGP)
  └── Table4.1_Example1.R   # Monte Carlo scripts for Tables 4.1 through 4.8
      ... [Table4.8_Example8.R]

/Empirical
  ├── empirical_dataclean.R                 # Data-cleaning script for IPUMS extract
  ├── Table5.1_descriptive_statistics.R     # Generates Summary Statistics
  ├── Table5.2&D.1_sex_unconditional_2way.R # 2-Way analysis scripts
  ├── ...
  └── TableD.8_race_by_education_3way.R     # 3-Way robustness scripts
```

---

## 3. Empirical Data Preparation

Due to copyright restrictions, the raw empirical dataset is **not included**. Users must obtain the data from **IPUMS USA** separately.

### (1) Obtain Raw Data
1. Visit [IPUMS USA](https://usa.ipums.org/usa/).
2. Download an extract containing the following samples: **2002 ACS, 2012 ACS, and 2022 ACS**.
3. Your extract **must include** these 15 variables:
   `YEAR`, `SAMPLE`, `CPI99`, `NCHILD`, `SEX`, `AGE`, `RACE`, `RACED`, `EDUC`, `EDUCD`, `OCC`, `IND`, `INCTOT`, `INCWAGE`, `PWSTATE2`.
4. Save the raw file (e.g., `usa_00026.csv.gz`) into the project root or the `/Empirical` folder.

### (2) Data Cleaning
Run `Empirical/empirical_dataclean.R` to process the raw extract. This script performs:
- **Sample Selection:** Restricts age to 18–65 and keeps only White (`RACE=1`) and Black (`RACE=2`) individuals.
- **Filtering:** Removes observations with missing values or invalid occupation/state codes.
- **Variable Construction:**
  - `real_Inwage`: Calculated as `(INCWAGE * CPI99 / 10000)`, representing pre-tax wage income in 1999 dollars.
  - `child`: Binary indicator for parental status (`NCHILD > 0`).
  - `age_group` & `edu_group`: Categorical variables for conditional testing.

---

## 4. Detailed Instructions for Reproduction

> [!IMPORTANT]
> All scripts rely on **relative paths**. You **must** set your R working directory to the **PROJECT ROOT DIRECTORY** (the folder containing this README) before execution.

### Execution Examples
**In R Console:**
```r
source("Simulations/Table4.1_Example1.R")
```

**In Command Line / Shell:**
```bash
Rscript Simulations/Table4.1_Example1.R
Rscript Empirical/Table5.2\&D.1_sex_unconditional_2way.R
```

---

### Suggested Running Order

1. **Step 1: Monte Carlo Simulations (Tables 4.1 – 4.8)**
   Run the scripts in `/Simulations`. Results will be saved as `.csv` files in a timestamped output folder.
   
2. **Step 2: Prepare Empirical Data**
   Run `Empirical/empirical_dataclean.R`. This script generates the cleaned dataset used for all subsequent empirical analysis.

3. **Step 3: Descriptive Statistics (Table 5.1)**
   Run `Empirical/Table5.1_descriptive_statistics.R`.

4. **Step 4: 2-Way Empirical Analysis (Tables 5.2 – 5.5, and D.1 – D.3)**
   Run the 2-way analysis scripts in the `/Empirical` folder.

5. **Step 5: 3-Way Empirical Analysis (Tables D.4 – D.8)**
   Run the 3-way analysis scripts in the `/Empirical` folder to reproduce the appendix robustness checks.

---

## Contact
For any technical inquiries regarding this replication package, please contact:
**Yaqian Wu** 
Email: [wuyq2024@hust.edu.cn](mailto:wuyq2024@hust.edu.cn)

---
