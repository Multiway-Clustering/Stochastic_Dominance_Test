

#N1=N2=50;M1=M2=50;mu_y=0;sigma_y=1

#Example 1 & 2
data_gen1_2<-function(N1,N2,M1,M2,mu_y,sigma_y){
  
  N=N1*N2;M=M1*M2
  
  y1<-rep(rnorm(N1,mu_y,sigma_y),each=N2)  #row common shock
  y2<-rep(rnorm(N2,mu_y,sigma_y),times=N1) #column common shock
  y12<-rnorm(N,mu_y,sigma_y)               #individual shock         
  Y<-(y1+y2+sqrt(3)*y12)/sqrt(5)

  x1<-rep(rnorm(M1),each=M2)
  x2<-rep(rnorm(M2),times=M1)
  x12<-rnorm(M)
  X<-(x1+x2+sqrt(3)*x12)/sqrt(5)
  
  Y<-pnorm(Y);X<-pnorm(X)
  return(list(X=X,Y=Y))
}

#Example 3&4
data_gen3_4<-function(N1,N2,M1,M2,mu_y,sigma_y){
  
  N=N1*N2;M=M1*M2
  
  Ncell_Y<-1+rpois(N,4)
  Ncell_X<-1+rpois(M,4) 
  
  simple_size_Y<-sum(Ncell_Y)   
  simple_size_X<-sum(Ncell_X)  
  
  Ncell_Y_mat<-matrix(Ncell_Y,N1,N2,byrow = TRUE)
  Ncell_X_mat<-matrix(Ncell_X,M1,M2,byrow = TRUE)
  
  y1<-rep(rep(rnorm(N1,mu_y,sigma_y),each=N2),times=Ncell_Y)    #row common shock
  y2<-rep(rep(rnorm(N2,mu_y,sigma_y),times=N1),times=Ncell_Y)    #column common shock
  y12<-rep(rnorm(N,mu_y,sigma_y),times=Ncell_Y)                  #cell shock
  y13<-rnorm(sum(Ncell_Y),mu_y,sigma_y)                          #individual shock         
  Y<-(y1+y2+sqrt(2)*y12+y13)/sqrt(5)
  
  x1<-rep(rep(rnorm(M1,0,1),each=M2),times=Ncell_X)    #row common shock
  x2<-rep(rep(rnorm(M2),times=M1),times=Ncell_X)    #column common shock
  x12<-rep(rnorm(M),times=Ncell_X)                  #cell shock
  x13<-rnorm(sum(Ncell_X))                          #individual shock         
  X<-(x1+x2+sqrt(2)*x12+x13)/sqrt(5)
  
  Y<-pnorm(Y);X<-pnorm(X)
  
  cluster_Y_ind1<-rep(1:N1,times=rowSums(Ncell_Y_mat))
  cluster_Y_ind2<-rep(rep(1:N2,times=N1),times=Ncell_Y)
  
  cluster_X_ind1<-rep(1:M1,times=rowSums(Ncell_X_mat))
  cluster_X_ind2<-rep(rep(1:M2,times=M1),times=Ncell_X)
  
  
  data_x<-data.frame(X,cluster_X_ind1,cluster_X_ind2)
  data_y<-data.frame(Y,cluster_Y_ind1,cluster_Y_ind2)
  return(list(data_x=data_x,data_y=data_y))
}

#Example 5 & 6
data_gen5_6<-function(N1,N2,M1,M2,mu_y,sigma_y,rho){
  
  if(N1!=M1 & N2!=M2){
    M1=N1; M2=N2
    cat('N1(N2) need to be equal to M1(M2)')
  }
  
  N=N1*N2;M=M1*M2
  
  x1<-rep(rnorm(M1),each=M2)
  x2<-rep(rnorm(M2),times=M1)
  x12<-rnorm(M)
  X<-(x1+x2+sqrt(3)*x12)/sqrt(5)
  
  u1<-rep(rnorm(M1),each=M2)
  u2<-rep(rnorm(M2),times=M1)
  u12<-rnorm(M)
  
  y1<-rho*x1+sqrt(1-rho^2)*u1
  y2<-rho*x2+sqrt(1-rho^2)*u2
  y12<-rho*x12+sqrt(1-rho^2)*u12
  
  
  Y<-sigma_y*(y1+y2+sqrt(3)*y12)/sqrt(5)+mu_y
  
  Y<-pnorm(Y);X<-pnorm(X)
  return(list(X=X,Y=Y))
}
