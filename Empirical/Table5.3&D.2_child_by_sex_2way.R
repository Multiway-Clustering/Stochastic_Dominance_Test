####==================================Data==================================####
data<-read.csv("Empirical/data_occ_ind_state_sex_child_edu_age.csv")
data<-data[,-1]

####=============================Report Result===============================####
library(tidyverse)

meta_data <- expand.grid(
  Year = c(2002, 2012, 2022),
  Sex  = c(1, 2), # 1: Male, 2: Female
  Race = c(1, 2), # 1: White, 2: Black
  Direction = c("With_vs_Without", "Without_vs_With")
) %>% 
  arrange(Direction, Race, Sex, Year) 

col_names <- c("SD1", "SD2", "SD3", "DH1", "DH2", "DH3")
colnames(Test_stat_M) <- col_names
colnames(pvalues_M) <- col_names

df_stat <- as.data.frame(Test_stat_M) %>%
  bind_cols(meta_data) %>%
  mutate(Metric = "Test Stat")

df_p <- as.data.frame(pvalues_M) %>%
  bind_cols(meta_data) %>%
  mutate(Metric = "p values")

final_df <- bind_rows(df_stat, df_p) %>%
  pivot_longer(cols = all_of(col_names), names_to = "Order", values_to = "Val") %>%
  pivot_wider(names_from = c(Direction, Order), values_from = Val) %>%
  mutate(
    Race_Lab = factor(ifelse(Race == 2, "Black", "White"), levels = c("Black", "White")),
    Sex_Lab  = factor(ifelse(Sex == 1, "Male", "Female"), levels = c("Male", "Female")),
    Metric   = factor(Metric, levels = c("Test Stat", "p values"))
  ) %>%
  arrange(Race_Lab, Sex_Lab, Year, Metric)

table_output <- final_df %>%
  mutate(across(where(is.numeric) & !c(Year, Sex, Race), ~ sprintf("%.3f", .))) %>%
  group_by(Race_Lab) %>%
  mutate(Race_Display = ifelse(row_number() == 1, as.character(Race_Lab), "")) %>%
  group_by(Race_Lab, Sex_Lab) %>%
  mutate(Sex_Display = ifelse(row_number() == 1, as.character(Sex_Lab), "")) %>%
  group_by(Race_Lab, Sex_Lab, Year) %>%
  mutate(Year_Display = ifelse(row_number() == 1, as.character(Year), "")) %>%
  ungroup() %>%
  select(
    Race = Race_Display, Sex = Sex_Display, Year = Year_Display, ` ` = Metric,
    starts_with("With_vs_Without"),
    starts_with("Without_vs_With")
  )

print(table_output)

write.csv(table_output, "2_2way_Empirical_Results_child_by_sex.csv", row.names = FALSE)