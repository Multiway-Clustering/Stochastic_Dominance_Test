library(data.table)

data<-read.csv("Empirical/data_occ_ind_state_sex_child_edu_age.csv")

data<-data[,-1]
setDT(data)


dt <- data[RACE %in% c(1, 2), ]

dt[, RACE_label := ifelse(RACE == 1, "White", "Black")]
dt[, SEX_label := ifelse(SEX == 1, "Male", "Female")]

calc_summary <- function(d) {
  ind_counts <- d[, .N, by = IND]$N
  occ_counts <- d[, .N, by = OCC]$N
  
  return(list(
    Sample_Size   = nrow(d),   
    Income_Mean   = mean(d$real_Inwage, na.rm = TRUE),
    Income_Median = median(d$real_Inwage, na.rm = TRUE),
    Income_SD     = sd(d$real_Inwage, na.rm = TRUE),
    
    Child_Mean    = mean(d$NCHILD, na.rm = TRUE),
    Child_Ratio   = mean(d$NCHILD > 0, na.rm = TRUE),
    
    Age_Mean      = mean(d$AGE, na.rm = TRUE),
    Age_SD        = sd(d$AGE, na.rm = TRUE),
    
    Edu_Mean      = mean(d$EDUC, na.rm = TRUE),
    Edu_SD        = sd(d$EDUC, na.rm = TRUE),
    
    Ind_Num       = uniqueN(d$IND),
    #Ind_SD        = if(length(ind_counts) > 1) sd(ind_counts) else 0,
    Occ_Num       = uniqueN(d$OCC),
    #Occ_SD        = if(length(occ_counts) > 1) sd(occ_counts) else 0
    State_Num     = uniqueN(d$PWSTATE2)
  ))
}


stat_sex <- dt[, calc_summary(.SD), by = .(YEAR, RACE_label, SEX_label)]

stat_total <- dt[, calc_summary(.SD), by = .(YEAR, RACE_label)]
stat_total[, SEX_label := "Total"] 

final_stat <- rbind(stat_sex, stat_total, use.names = TRUE)


final_stat[, SEX_label := factor(SEX_label, levels = c("Female", "Male", "Total"))]
setorder(final_stat, YEAR, RACE_label, SEX_label)
print(head(final_stat, 10))

write.csv(final_stat, 'desc_stat_by_race_sex.csv', row.names = FALSE)
