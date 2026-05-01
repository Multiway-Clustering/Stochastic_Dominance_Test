####==================================Data==================================####
data<-read.csv("Empirical/data_occ_ind_state_sex_child_edu_age.csv")
data<-data[,-1]

####===============================Parallel Setup===========================####

library(foreach)
library(doSNOW)
library(doParallel)
library(parallel)
library(data.table)
source('R/myFUN.R')
no_cores <- detectCores()-1
cl = makeCluster(no_cores, type = "SOCK")
registerDoSNOW(cl)
clusterSetRNGStream(cl,)


clusterEvalQ(cl, {
  library(data.table)
  source('R/myFUN.R')
  
})

####=================================Test==================================####

reps = 499
year_set<-c(2002,2012,2022)
race_set<-c(1,2)
sex_set<-rbind(c(1,2),c(2,1))
alpha<-0.05
Test_stat_M<-NULL
pvalues_M<-NULL

for (sex in 1:2) {
  cat(sprintf("\n--- Running for SEX = %d vs. SEX = %d ---\n", sex_set[sex,1],sex_set[sex,2]))
  for (race in race_set){
    cat(sprintf("  --> Running for RACE = %d ---\n", race))
    for (year in year_set) {
      cat(sprintf("    ->> Running  for YEAR = %d\n", year))
      
      data_now<-data[data$YEAR==year,]
      
      data1<-data_now[data_now$SEX==sex_set[sex,1]&data_now$RACE==race,]
      data2<-data_now[data_now$SEX==sex_set[sex,2]&data_now$RACE==race,]
      
      N1<-length(unique(data1$IND));M1<-length(unique(data2$IND));
      N2<-length(unique(data1$OCC));M2<-length(unique(data2$OCC))
      Nbar<-min(N1,N2);Mbar<-min(M1,M2)
      
      ####============================Test Statistic==============================####
      
      y1<-data1$real_Inwage
      y2<-data2$real_Inwage
      
      z<-sort(unique(c(y1,y2)))
      Fy1<-myECDF2(y1,z)$Fy
      Fy2<-myECDF2(y2,z)$Fy
      
      plot(z,Fy1,type='l',ylab="ECDF",xlab='Income ($10,000)',main=c(year))
      lines(z,Fy2,lty=4,col='blue')
      legend("bottomright",legend = c('Male','Female'),lty=c(1,4),col=c('black','blue')) 
      
      FFy1 <- myFFy(y1,z)$FFy
      FFy2 <- myFFy(y2,z)$FFy
      
      FFFy1 <- myFFFy(y1,z)$FFFy
      FFFy2 <- myFFFy(y2,z)$FFFy
      
      TN1<-sqrt(Nbar*Mbar/(Nbar+Mbar))*max(Fy1-Fy2)
      TN2<-sqrt(Nbar*Mbar/(Nbar+Mbar))*max(FFy1-FFy2)
      TN3<-sqrt(Nbar*Mbar/(Nbar+Mbar))*max(FFFy1-FFFy2)
      
      N<-as.numeric(length(y1));M<-as.numeric(length(y2))
      
      TN1_DH<-sqrt(N*M/(N+M))*max(Fy1-Fy2)
      TN2_DH<-sqrt(N*M/(N+M))*max(FFy1-FFy2)
      TN3_DH<-sqrt(N*M/(N+M))*max(FFFy1-FFFy2)
      
      Test_statistics<-c(TN1,TN2,TN3,TN1_DH,TN2_DH,TN3_DH)
      Test_stat_M<-rbind(Test_stat_M,Test_statistics)
      #round(Test_statistics,3)
      
      ####==========================Recentering Function==========================####
      
      delta1<--0.1*sqrt(log(log(Nbar+Mbar)))
      #delta2<--0.01*sqrt(log(log(Nbar+Mbar)))
      #delta3<--0.001*sqrt(log(log(Nbar+Mbar)))
      #delta4<-0
      
      delta<-c(delta1)#,delta2,delta3,delta4)
      
      mu1<-vapply(delta, function(i) 
        (Fy1-Fy2)*((sqrt(Nbar*Mbar/(Nbar+Mbar))*(Fy1-Fy2))<i),numeric(length(z)))
      mu2<-vapply(delta, function(j) 
        (FFy1-FFy2)*((sqrt(Nbar*Mbar/(Nbar+Mbar))*(FFy1-FFy2))<j),numeric(length(z)))
      mu3<-vapply(delta, function(k) 
        (FFFy1-FFFy2)*((sqrt(Nbar*Mbar/(Nbar+Mbar))*(FFFy1-FFFy2))<k),numeric(length(z)))
      
      delta1_DH<--0.1*sqrt(log(log(N+M)))

      delta_DH<-c(delta1_DH)#,delta2_DH,delta3_DH,delta4_DH)
      
      mu1_DH<-vapply(delta_DH, function(i) 
        (Fy1-Fy2)*((sqrt(N*M/(N+M))*(Fy1-Fy2))<i),numeric(length(z)))
      mu2_DH=vapply(delta_DH, function(j) 
        (FFy1-FFy2)*((sqrt(N*M/(N+M))*(FFy1-FFy2))<j),numeric(length(z)))
      mu3_DH=vapply(delta_DH, function(k) 
        (FFFy1-FFFy2)*((sqrt(N*M/(N+M))*(FFFy1-FFFy2))<k),numeric(length(z)))
      
      
      ####============================Bootstrap===================================####
      
      system.time({
        result =  foreach(k = 1:reps, .combine = 'rbind')%dopar%{
          
          #####================================WFH=================================#####
          Y1_bs<-bs_sample_empirical('IND','OCC','real_Inwage',data1)
          Fy1_bs<-myECDF2(Y1_bs,z)$Fy
          
          Y2_bs<-bs_sample_empirical('IND','OCC','real_Inwage',data2)
          Fy2_bs<-myECDF2(Y2_bs,z)$Fy
          
          FFy1_bs <- myFFy(Y1_bs,z)$FFy
          FFFy1_bs <- myFFFy(Y1_bs,z)$FFFy
          
          FFy2_bs <- myFFy(Y2_bs,z)$FFy
          FFFy2_bs <- myFFFy(Y2_bs,z)$FFFy
          
          
          TN_bs1<-apply(mu1,2,function(mu) sqrt(Nbar*Mbar/(Nbar+Mbar))*max((Fy1_bs-Fy2_bs)-(Fy1-Fy2)+mu))
          TN_bs2<-apply(mu2,2,function(mu) sqrt(Nbar*Mbar/(Nbar+Mbar))*max((FFy1_bs-FFy2_bs)-(FFy1-FFy2)+mu))
          TN_bs3<-apply(mu3,2,function(mu) sqrt(Nbar*Mbar/(Nbar+Mbar))*max((FFFy1_bs-FFFy2_bs)-(FFFy1-FFFy2)+mu))
          
          
          #####==================================DH================================#####
          y1_bs_DH<-y1[sample(1:N,N,replace = TRUE)]
          Fy1_bs_DH<-myECDF2(y1_bs_DH,z)$Fy
          
          y2_bs_DH<-y2[sample(1:M,M,replace = TRUE)]
          Fy2_bs_DH<-myECDF2(y2_bs_DH,z)$Fy
          
          FFy1_bs_DH <- myFFy(y1_bs_DH,z)$FFy
          FFFy1_bs_DH <- myFFFy(y1_bs_DH,z)$FFFy
          
          FFy2_bs_DH <- myFFy(y2_bs_DH,z)$FFy
          FFFy2_bs_DH <- myFFFy(y2_bs_DH,z)$FFFy
          
          TN_bs_DH1<-apply(mu1_DH,2,function(mu) sqrt(N*M/(N+M))*max((Fy1_bs_DH-Fy2_bs_DH)-(Fy1-Fy2)+mu))
          TN_bs_DH2<-apply(mu2_DH,2,function(mu) sqrt(N*M/(N+M))*max((FFy1_bs_DH-FFy2_bs_DH)-(FFy1-FFy2)+mu))
          TN_bs_DH3<-apply(mu3_DH,2,function(mu) sqrt(N*M/(N+M))*max((FFFy1_bs_DH-FFFy2_bs_DH)-(FFFy1-FFFy2)+mu))
          
          RES<-c(TN_bs1,TN_bs2,TN_bs3,TN_bs_DH1,TN_bs_DH2,TN_bs_DH3)
          
          return(RES)
        }
      })  

      
      ####==============================p value===================================####
      Test_statistics2<-rep(Test_statistics,each=length(delta))
      
      pvalues<-(rowSums(apply(result,1,function(x) x>Test_statistics2))+1)/(reps+1)
      pvalues_M<-rbind(pvalues_M,pvalues)
      
    }
  }
}


####=============================Report Result===============================####
library(tidyverse)
meta_data <- expand.grid(
  Year      = year_set, # c(2002, 2012, 2022)
  Race      = race_set, # c(1, 2)
  Direction = c("Male_vs_Female", "Female_vs_Male")
)

col_names <- c("SD1", "SD2", "SD3", "DH1", "DH2", "DH3")
colnames(Test_stat_M) <- col_names
colnames(pvalues_M)   <- col_names

df_stat <- as.data.frame(Test_stat_M) %>%
  bind_cols(meta_data) %>%
  mutate(Metric = "Test Stat")

df_p <- as.data.frame(pvalues_M) %>%
  bind_cols(meta_data) %>%
  mutate(Metric = "p values")

combined_df <- bind_rows(df_stat, df_p) %>%
  pivot_longer(cols = all_of(col_names), names_to = "Order", values_to = "Val") %>%
  pivot_wider(names_from = c(Direction, Order), values_from = Val)

final_table <- combined_df %>%
  mutate(
    Race_Lab = factor(ifelse(Race == 2, "Black", "White"), levels = c("Black", "White")),
    Metric   = factor(Metric, levels = c("Test Stat", "p values"))
  ) %>%
  arrange(Race_Lab, Year, Metric) %>%
  mutate(across(where(is.numeric) & !Year, ~ sprintf("%.3f", .)))

table_output <- final_table %>%
  group_by(Race_Lab) %>%
  mutate(Race_Display = ifelse(row_number() == 1, as.character(Race_Lab), "")) %>%
  group_by(Race_Lab, Year) %>%
  mutate(Year_Display = ifelse(row_number() == 1, as.character(Year), "")) %>%
  ungroup() %>%
  select(
    Race = Race_Display, 
    Year = Year_Display, 
    ` ` = Metric,
    starts_with("Male_vs_Female"),
    starts_with("Female_vs_Male")
  )

print(table_output)
 write.csv(table_output, "0_2way_Emprical_Gender_Dominance_By_Race.csv", row.names = FALSE)
