ISTEP <- function(case,d){
  ipw <- IPW(case[,2],case[,1],case[,-c(1,2)])
  finite_rows <- apply(ipw, 1, function(x) all(is.finite(x)))
  # 选择不包含Inf值的行
  ipw <- ipw[finite_rows, ]
  
  SIS <- rank.sis.lrm(d,ipw);
  MicSIS <- Mic.SIS(d,ipw);
  MrMicSIS <- MrMic.SIS(d,ipw);
  DPDSIS <- dpd.sis(d,ipw);
  MrDcSIS <- MrDc.SIS(d,ipw);
  I <- cbind(SIS[,1],MicSIS[,1],MrMicSIS[,1],DPDSIS[,1],MrDcSIS[,1]);
  TU <- cbind(SIS[1,3],MicSIS[1,3],MrMicSIS[1,3],DPDSIS[1,3],MrDcSIS[1,3])
  count <- matrix(rep(0,5*5),ncol =5 )
  for (i in c(1:4)) {
    for (j in c(1:5)) {
      if(i %in% I[,j]){
        count[i,j] = 1
      }
    }
  }
  for (j in c(1:5)) {
    count[5,j] = sum(count[,j])/4
  }
  result <- list(count=count,TU = TU)
  return(result)
  
}

diedai <-function(d,mu,SigmaC1,U1,betaC1){
  M = matrix(rep(0,5*5),ncol =5 )
  
  library("doParallel")
  library("foreach")
  nprocs=detectCores() - 1;
  cl=makeCluster(nprocs);
  registerDoParallel(cl);
  TU = rep(0,4);
  est<-foreach(w=1:500,
               .combine=rbind,
               .packages =c("sfsmisc","MASS","mvtnorm","sn","stats4","minerva","energy","gsl"),
               .export = c("Moni1","IPW","ISTEP","equals","Sc","MrDc.SIS","MrMic.SIS","dpd.sis",
                           "Censor1","Mic.SIS","rank.sis.lrm","ab")
  ) %do% {
    Case1 <- Moni1(betaC1,300,sig = SigmaC1,mu,U1);
    Case1 <- na.omit(Case1)
    CO <- ISTEP(Case1,d);
    M = M + CO$count
    TU = rbind(TU,CO$TU)
    if(w%%5 == 0){
      print(w)
    }
  }  #close for cycle w
  stopCluster(cl)
  M = M/500;
  result = list(M=M[1:4,],TU=TU[-1,])
  return(result)
}

Sigma = diag(c(2,4,rep(1,998)));
for (i in c(1:1000)) {
  for (j in c(1:1000)) {
    if(i!=j){
      Sigma[i,j] = 0.5^abs(i-j)
    }
  }
}

betaC1 = c(3,1,3,-4,rep(0,996))
betaC2 = c(4,-1,1,3,rep(0.996))
U1 = 0.3
U2 = 0.6

case1res = diedai(300,c(3,2,1,1,rep(0,996)),Sigma,0.3,c(3,1,3,-4,rep(0,996)))
