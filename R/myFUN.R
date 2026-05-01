myECDF<-function(Y,z){

  Iy <- outer(Y, z, `<=`)
  Fy<-colMeans(Iy)
  
  return(list(Fy=Fy))
}


myECDF2<-function(Y,z){

  Fy <- vapply(z, function(zj) mean(Y <= zj), numeric(1))
  
  return(list(Fy=Fy))
}

myFFy<-function(Y,z){
  
  FFy <- vapply(z, function(zj) {
    mean((Y <= zj) * (zj - Y))
  }, numeric(1))
  
  return(list(FFy=FFy))
}

myFFFy<-function(Y,z){
  
  FFFy <- vapply(z, function(zj) {
    mean((Y <= zj) * (zj - Y)^2/2)
  }, numeric(1))
  
  return(list(FFFy=FFFy))
}

#For 2-way, only one unit in each cell
bs_sample<-function(N1,N2,Y){
  Ym<-matrix(Y,N1,N2,byrow = TRUE)
  
  N1bs<-sample(1:N1,N1,replace=TRUE)
  N2bs<-sample(1:N2,N2,replace=TRUE)
  
  Ymbs<-Ym[N1bs,N2bs];
  Y_bs<-c(Ymbs)
  return(Y_bs=Y_bs)
}

#For 2-way, more than one unit in each cell
bs_sample2<-function(N1,N2,data_y){
  #Ym<-matrix(Y,N1,N2,byrow = TRUE)
  
  N1bs<-sample(1:N1,N1,replace=TRUE)
  N2bs<-sample(1:N2,N2,replace=TRUE)
  
  library(data.table)
  bootstrap_clusters <- CJ(cluster1_ind = N1bs, cluster2_ind = N2bs, unique = FALSE)
  
  dt<-as.data.table(data_y)
  names(dt)<-c('Y','cluster1_ind','cluster2_ind')
  y_bs <- dt[bootstrap_clusters, on = .(cluster1_ind, cluster2_ind), Y,allow.cartesian = TRUE]
  
  y_bs_clean<-na.omit(y_bs)
  return(Y_bs=y_bs_clean)
}

bs_sample_cor<-function(N1,N2,X,Y){
  Ym<-matrix(Y,N1,N2,byrow = TRUE)
  Xm<-matrix(X,N1,N2,byrow = TRUE)
  
  N1bs<-sample(1:N1,N1,replace=TRUE)
  N2bs<-sample(1:N2,N2,replace=TRUE)
  
  Ymbs<-Ym[N1bs,N2bs];
  Y_bs<-c(Ymbs)
  
  Xmbs<-Xm[N1bs,N2bs];
  X_bs<-c(Xmbs)
  return(list(Y_bs=Y_bs,X_bs=X_bs))
}

#IND<-data2022$IND;STATE<-data2022$PWSTATE2;DATA<-data2022

bs_sample_empirical<-function(cluster1_name,cluster2_name,y_name,DATA){
  
  dt<-as.data.table(DATA)
  c1_values <- unique(dt[[cluster1_name]])
  c2_values <- unique(dt[[cluster2_name]])
  
  ID1 = sample(c1_values, length(c1_values), replace = TRUE)
  bs1 <- data.table(ID1)
  bs1_counts <- bs1[, .N, by = ID1]
  setnames(bs1_counts, c(cluster1_name, "n1"))
  
  ID2 = sample(c2_values, length(c2_values), replace = TRUE)
  bs2 <- data.table(ID2)
  bs2_counts <- bs2[, .N, by = ID2]
  setnames(bs2_counts, c(cluster2_name, "n2"))
  
  dt_bs <- merge(dt, bs1_counts, by = cluster1_name, all = FALSE)
  dt_bs <- merge(dt_bs, bs2_counts, by = cluster2_name, all = FALSE)
  dt_bs[, bs_weight := n1 * n2]
  
  y_bs_res <- dt_bs[rep(1:.N, bs_weight), get(y_name)]
  
  return(y_bs_res)
}


bs_sample_empirical_3way <- function(cluster1_name, cluster2_name, cluster3_name, y_name, DATA) {
  dt <- as.data.table(DATA)
  
  c1_values <- unique(dt[[cluster1_name]])
  bs1 <- data.table(ID1 = sample(c1_values, length(c1_values), replace = TRUE))
  bs1_counts <- bs1[, .N, by = ID1]
  setnames(bs1_counts, c(cluster1_name, "n1"))
  
  c2_values <- unique(dt[[cluster2_name]])
  bs2 <- data.table(ID2 = sample(c2_values, length(c2_values), replace = TRUE))
  bs2_counts <- bs2[, .N, by = ID2]
  setnames(bs2_counts, c(cluster2_name, "n2"))
  
  c3_values <- unique(dt[[cluster3_name]])
  bs3 <- data.table(ID3 = sample(c3_values, length(c3_values), replace = TRUE))
  bs3_counts <- bs3[, .N, by = ID3]
  setnames(bs3_counts, c(cluster3_name, "n3"))
  

  dt_bs <- merge(dt, bs1_counts, by = cluster1_name, all = FALSE)
  dt_bs <- merge(dt_bs, bs2_counts, by = cluster2_name, all = FALSE)
  dt_bs <- merge(dt_bs, bs3_counts, by = cluster3_name, all = FALSE)
  
  dt_bs[, bs_weight := as.numeric(n1) * n2 * n3]
  
  y_bs_res <- dt_bs[rep(1:.N, bs_weight), get(y_name)]
  
  return(y_bs_res)
}
