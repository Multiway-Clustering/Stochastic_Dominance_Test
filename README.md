-----------------------------------------------------------------------------
REPLICATION PACKAGE FOR: A Multiway Cluster-Robust Test for Stochastic Dominance
-----------------------------------------------------------------------------
AUTHORS: Yaqian Wu, Qingliang Fan, and Yu-Chin Hsu



This replication package contains all the necessary R code, custom functions, and empirical datasets to fully reproduce the simulation and empirical results presented in the manuscript "A Multiway Cluster-Robust Test for Stochastic Dominance".

-----------------------------------------------------------------------------
1. SYSTEM REQUIREMENTS
-----------------------------------------------------------------------------
The code is written in R. Most scripts utilize parallel computing and will automatically detect available CPU cores to optimize runtime.

Required R packages include:
- data.table
- foreach
- parallel
- doParallel
- doSNOW
- tidyverse
- stringr

You can install all necessary packages by running the following command in R:
install.packages(c("data.table", "foreach", "doParallel", "doSNOW", "tidyverse", "stringr"))

-----------------------------------------------------------------------------
2. DIRECTORY STRUCTURE
-----------------------------------------------------------------------------
The "replication.zip" package is organized as follows:

/R
  - myFUN.R : Centralized file containing common utility functions (e.g., empirical CDFs, integrated CDFs, and bootstrap resampling functions).

/Simulations
  - sim_datagen.R             : Data-generating processes used across simulation scripts.
  - Table4.1_Example1.R       : Monte Carlo script for Table 4.1.
  - Table4.2_Example2.R       : Monte Carlo script for Table 4.2.
  - Table4.3_Example3.R       : Monte Carlo script for Table 4.3.
  - Table4.4_Example4.R       : Monte Carlo script for Table 4.4.
  - Table4.5_Example5.R       : Monte Carlo script for Table 4.5.
  - Table4.6_Example6.R       : Monte Carlo script for Table 4.6.
  - Table4.7_Example7.R       : Monte Carlo script for Table 4.7.
  - Table4.8_Example8.R       : Monte Carlo script for Table 4.8.

/Empirical
  - data_occ_ind_state_sex_child_edu_age.csv : The raw dataset used in the empirical analysis.
  - Table5.1_descriptive_statistics.R        
  - Table5.2&D.1_sex_unconditional_2way.R
  - Table5.2&D.1_sex_by_child_2way.R
  - Table5.3&D.2_child_by_sex_2way.R
  - Table5.4&D.3_child_by_age_2way.R
  - Table5.5_race_by_education_2way.R
  - TableD.4_sex_unconditional_3way.R
  - TableD.4_sex_by_child_3way.R
  - TableD.5_child_by_sex_3way.R
  - TableD.6&D.7_child_by_age_3way.R
  - TableD.8_race_by_education_3way.R

-----------------------------------------------------------------------------
3. DETAILED INSTRUCTIONS FOR REPRODUCTION
-----------------------------------------------------------------------------
IMPORTANT: All scripts rely on relative paths. You must set your R working directory to the PROJECT ROOT DIRECTORY (the folder containing this README) before executing any script.

Example to run a script in the R console:
> source("Simulations/Table4.1_Example1.R")

Example to run a script from the command line / shell:
$ Rscript Simulations/Table4.1_Example1.R
$ Rscript Empirical/Table5.2\&D.1_sex_unconditional_2way.R

--- SUGGESTED RUNNING ORDER ---
To systematically reproduce all tables in the main text and the appendix, please follow this sequence:

Step 1: Descriptive Statistics
Run `Empirical/Table5.1_descriptive_statistics.R` to generate the summary statistics.

Step 2: Monte Carlo Simulations (Tables 4.1 - 4.8)
Run the simulation scripts located in the `/Simulations` folder sequentially. 
Note: Simulation scripts may take substantial time to finish depending on the number of repetitions, bootstrap draws, and your hardware capabilities. Output files (.csv) will be automatically written to a timestamped output folder.

Step 3: 2-Way Empirical Analysis (Tables 5.2 to 5.5, and D.1 to D.3)
Run the corresponding 2-way empirical scripts in the `/Empirical` folder. Each script reads the provided empirical dataset and writes a `.csv` output file.

Step 4: 3-Way Empirical Analysis (Tables D.4 to D.8)
Run the corresponding 3-way empirical scripts in the `/Empirical` folder to reproduce the appendix tables.

-----------------------------------------------------------------------------

For any technical inquiries regarding this replication package, please contact:
Yaqian Wu at wuyq2024@hust.edu.cn

