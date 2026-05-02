
# REPLICATION PACKAGE FOR: A Multiway Cluster-Robust Test for Stochastic Dominance

**AUTHORS:** Yaqian Wu, Qingliang Fan, and Yu-Chin Hsu

This replication package contains all the necessary R code to reproduce the simulation and empirical results presented in the manuscript "A Multiway Cluster-Robust Test for Stochastic Dominance".

---

## 1. SYSTEM REQUIREMENTS

The code is written in R. Most scripts utilize parallel computing and will automatically detect available CPU cores to optimize runtime.

**Required R packages include:**
- `data.table`
- `foreach`
- `parallel`
- `doParallel`
- `doSNOW`
- `tidyverse`
- `stringr`

You can install all necessary packages by running the following command in R:
```r
install.packages(c("data.table", "foreach", "doParallel", "doSNOW", "tidyverse", "stringr"))
```

---

## 2. DIRECTORY STRUCTURE

The `replication.zip` package is organized as follows:

- **/R**
  - `myFUN.R`: Centralized file containing common utility functions (e.g., empirical CDFs, integrated CDFs, and bootstrap resampling functions).

- **/Simulations**
  - `sim_datagen.R`: Data-generating processes used across simulation scripts.
  - `Table4.1_Example1.R`: Monte Carlo script for Table 4.1.
  - `Table4.2_Example2.R`: Monte Carlo script for Table 4.2.
  - `Table4.3_Example3.R`: Monte Carlo script for Table 4.3.
  - `Table4.4_Example4.R`: Monte Carlo script for Table 4.4.
  - `Table4.5_Example5.R`: Monte Carlo script for Table 4.5.
  - `Table4.6_Example6.R`: Monte Carlo script for Table 4.6.
  - `Table4.7_Example7.R`: Monte Carlo script for Table 4.7.
  - `Table4.8_Example8.R`: Monte Carlo script for Table 4.8.

- **/Empirical**
  - `empirical_dataclean.R`: Data-cleaning script that prepares the empirical analysis dataset from the original IPUMS USA extract.
  - `Table5.1_descriptive_statistics.R`: Generates descriptive statistics presented in Table 5.1.
  - `Table5.2&D.1_sex_unconditional_2way.R`: Analysis for Table 5.2 and D.1 (the first unconditional panel).
  - `Table5.2&D.1_sex_by_child_2way.R`: Analysis for Table 5.2 and D.1 (the second and third panels).
  - `Table5.3&D.2_child_by_sex_2way.R`: Analysis for Table 5.3 and D.2.
  - `Table5.4&D.3_child_by_age_2way.R`: Analysis for Table 5.4 and D.3.
  - `Table5.5_race_by_education_2way.R`: Analysis for Table 5.5.
  - `TableD.4_sex_unconditional_3way.R`: Analysis for Table D.4 (the first unconditional panel).
  - `TableD.4_sex_by_child_3way.R`: Analysis for Table D.4 (the second and third panels).
  - `TableD.5_child_by_sex_3way.R`: Analysis for Table D.5.
  - `TableD.6&D.7_child_by_age_3way.R`: Analysis for Table D.6 and D.7.
  - `TableD.8_race_by_education_3way.R`: Analysis for Table D.8.

---

## 3. EMPIRICAL DATA PREPARATION

The empirical dataset is not included in this replication package for copyright reasons. Users need to obtain the IPUMS USA data (Ruggles et al., 2023) separately and prepare the empirical dataset before running the empirical scripts.

### (1) Download the raw data (e.g., a file named "usa_00026.csv.gz").
- The data are available from IPUMS USA at: https://usa.ipums.org/usa/
- Users should download an IPUMS USA extract based on the following samples:
  - 2002 ACS, 2012 ACS, 2022 ACS
- The extract must includes the following variables:
  - `YEAR`, `SAMPLE`, `CPI99`, `NCHILD`, `SEX`, `AGE`, `RACE`, `EDUC`, `OCC`, `IND`, `INCWAGE`, `PWSTATE2`

### (2) Run the script `Empirical/empirical_dataclean.R` to perform the following cleaning and construction steps:
- Restrict sample to individuals aged 18 to 65.
- Keep only White (`RACE = 1`) and Black (`RACE = 2`) individuals.
- Restrict to valid occupation, industry, and state-of-work codes.
- Drop missing values.
- **Variable Construction:**
  - `real_Inwage`: Constructed as `(INCWAGE * CPI99 / 10000)`. This represents pre-tax wage income adjusted to 1999 dollars.
  - `child`: Dummy variable for parental status (`NCHILD > 0`).
  - `age_group` & `edu_group`: Discretized categories for conditional tests

---

## 4. DETAILED INSTRUCTIONS FOR REPRODUCTION

**IMPORTANT:** All scripts rely on relative paths. You must set your R working directory to the **PROJECT ROOT DIRECTORY** (the folder containing this README) before executing any script.

**Example to run a script in the R console:**
```r
source("Simulations/Table4.1_Example1.R")
```

**Example to run a script from the command line / shell:**
```bash
$ Rscript Simulations/Table4.1_Example1.R
$ Rscript Empirical/Table5.2\&D.1_sex_unconditional_2way.R
```

### --- SUGGESTED RUNNING ORDER ---
To systematically reproduce all tables in the main text and the appendix, please follow this sequence:

- **Step 1: Monte Carlo Simulations (Tables 4.1 - 4.8)**
  Run the simulation scripts located in the `/Simulations` folder sequentially. Output files (`.csv`) will be automatically written to a timestamped output folder.

- **Step 2: Prepare Empirical Data**
  Download the required IPUMS USA data and run `Empirical/empirical_dataclean.R`. This script processes the raw IPUMS extract and saves the cleaned dataset required for all subsequent empirical tables.

- **Step 3: Descriptive Statistics (Tables 5.1)**
  Run `Empirical/Table5.1_descriptive_statistics.R` to generate the summary statistics.

- **Step 4: 2-Way Empirical Analysis (Tables 5.2 to 5.5, and D.1 to D.3)**
  Run the corresponding 2-way empirical scripts in the `/Empirical` folder. Each script reads the provided empirical dataset and writes a `.csv` output file.

- **Step 5: 3-Way Empirical Analysis (Tables D.4 to D.8)**
  Run the corresponding 3-way empirical scripts in the `/Empirical` folder to reproduce the appendix tables

---

**For any technical inquiries regarding this replication package, please contact:**
Yaqian Wu at wuyq2024@hust.edu.cn

---
