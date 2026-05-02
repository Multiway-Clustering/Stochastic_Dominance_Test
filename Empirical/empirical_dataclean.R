
####============================clean data=================================####
data_raw<-read.csv("/usa_00026.csv.gz")
size1<-dim(data_raw)[1];size1

# In IPUMS coding, 0 or special top-coded values indicate missing or invalid observations for selected variables.
# Drop observations with missing wage income coded as 999999.
# Restrict the sample to working-age adults between 18 and 65.
# Keep only White (1) and Black (2) individuals.
# Sex is coded as 1 for male and 2 for female.
library(tidyverse)

data <- data_raw %>%
  filter(OCC > 0 & PWSTATE2 >= 1 & PWSTATE2 <= 56 & INCWAGE!=999999 & AGE >= 18 
         & AGE <= 65 & EDUC != 99 & RACE <=2) %>%
  mutate(
    # 1. Create age-group categories.
    age_group = case_when(
    AGE <= 24 ~ "1_18-24",
    AGE <= 34 ~ "2_25-34",
    AGE <= 49 ~ "3_35-49",
    AGE <= 65 ~ "4_50-65"),
    # 2. Create education-group categories.
    edu_group = case_when(
      EDUC >= 0 & EDUC <= 5 ~ "1_Less_than_HS",
      EDUC == 6             ~ "2_HS_Graduate",
      EDUC >= 7 & EDUC <= 8 ~ "3_Some_College",
      EDUC == 10            ~ "4_Bachelors",
      EDUC == 11            ~ "5_Advanced_Degree"),
   # 3. Deflate wage income and convert the unit to ten-thousands.
   real_Inwage = INCWAGE * CPI99 / 10000,
  
   # 4. Create an indicator for having at least one child.
   child = if_else(NCHILD > 0, 1, 0) 
  )


final_data <- data %>%
  select(
    YEAR, NCHILD, SEX, AGE, RACE, EDUC, OCC, IND, PWSTATE2, 
    age_group, edu_group, real_Inwage, child
  )

str(final_data)

write.csv(final_data,"Empirical/data_occ_ind_state_sex_child_edu_age.csv")
