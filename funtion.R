library(minerva)  #MIC
#Parallel algorithms
library(foreach)
library(doParallel)
library(iterators)
# calculate DC
library(energy)
library(gsl)

library(MASS)  
library(sn)
library(stats4)
library(mvtnorm)

library(glmnet)
library(survival)
library(bujar)

library(dplyr)
library(splines)

library(ncvreg)


# first step

#1,simulated data censoring
Censor1 <- function(U1,Y,n){
  U2 <- runif(n,min = 0,max = 1);
  for (i in 1:n) {
    if(U2[i] <= U1){
      Y[i,2] = 0
      Y[i,1] = log(runif(1,min = 0 + .Machine$double.eps,max = exp(Y[i,1]) ))
    }else{
      Y[i,2] = 1
    }
  }
  return(Y)
}

#2,mapping
equals <- function(x){
  n=length(x);
  xr <- rep(0,n);
  xr[order(x)] <- c(1:length(x));
  equals = seq(0,1,length.out = n);
  iassign <- equals[xr];
  return(iassign)
}

#3,  predators
xxt <- function(e){
  a = sum(abs(e)>=abs(e[1]))
  b = sum(abs(e)>=abs(e[2]))
  c = sum(abs(e)>=abs(e[3]))
  d = sum(abs(e)>=abs(e[4]))
  return(max(a,b,c,d))
}

#4,KM estimator
Sc <-function(delta,X){
  #n <- dim(X)[1]
  n <- length(X)
  Xsort = sort(X);
  index <- order(X);
  delta_ = 1 - delta;
  Sc <- rep(1,n);
  Sc[1] <- 1 - delta_[index[1]]/sum(X>=Xsort[1]);
  for (i in 2:n) {
    Sc[i] = Sc[i-1]*(1-delta_[index[i]]/sum(X>=Xsort[i]))
  }
  
  Sc[n] = Sc[n-1];
  Sc <- as.data.frame(Sc);
  Sc$index <-index;
  return(Sc)
}

#5,

IPW <- function(delta,y,x){
  Sc <- Sc(1- delta,y);
  Yipw <- data.frame(y,delta,c(1:length(y)));
  colnames(Yipw) <- c("y","delta","number");
  KM1 <- Sc$Sc[order(Sc$index)];
  Yipw$y <- Yipw$delta*log(exp(Yipw$y)/KM1);
  Yipw <- Yipw[Yipw$delta == 1,]
  Yipwc <- Yipw[,1];
  Xipw <- x[Yipw$number,];
  result <- data.frame(Yipwc,Xipw)
  return(result)
}

#6,
MrDc.SIS <- function(d,ipw){
  ipwzr <- apply(ipw, 2, equals)
  m = rep(0,dim(ipw)[2]-1)
  for (i in 2:dim(ipw)[2]) {
    m[i-1] = dcor(ipwzr[,1],ipwzr[,i])
  }
  M = rep(0,d)
  M[1] = xxt(m)
  I=sort(abs(m), decreasing = TRUE, index.return=TRUE)
  I1=I$ix[1:d]
  beta=m[I1]
  result=cbind(I1, beta,M)
  return(result)
}







# step two

#7,soft threshold function
SoftThreshold  <- function(c,t){
  xx <- pmax(abs(c)-t,0);
  return(sign(c)*xx)
}



#8,Define the p-norm
norm_p <- function(x,p) {
  norm_p <- sum(abs(x)^p)^(1/p)
  return(norm_p)
}



#9,Define the Adaptive step size 
alpha <- function(x,p){
  #x <- drop(x);
  alpha_p <- norm_p(x,1)/norm_p(x,p);
  return(alpha_p)
}



#10,Define the unit vector
e <- function(pos,n){
  e <- as.matrix(rep(0,n));
  e[pos,] <- 1;
  return(e)
}



#11,Define the numerator of the gradient under the p-norm
fenzi <- function(x,p){
  n <- length(x);
  fenzi <- 0;
  for (i in 1:n) {
    fenzi = fenzi + abs(x[i])**(p-1) * sign(x[i])*e(i,n)
  }
  return(fenzi)
}

#12,
#Define the gradient under the p-norm counterpart
gradient <- function(la,x,p){
  alpha_l = alpha(x,p);
  fenzi_1 <- fenzi(x,p);
  v <- la*alpha_l*fenzi_1/(norm_p(x,p))**(p-1)
  return(v)
}

# 13
BJ <- function(X,Y,delta,beta){
  status <- delta;
  yy <- Y; 
  N <- length(yy); 
  timeorig <- yy; 
  order.orig <- 1:N;
  dummystrat <- factor(rep(1, N));  
  ypred <- X%*%beta;
  ehat <- timeorig - ypred; 
  state <- status;
  state[ehat == max(ehat)] <- 1;
  S <- structure(cbind(ehat, state), class = "Surv", type = "right");
  KM.ehat <- survfitKM(dummystrat, S, conf.type = "none", se.fit = FALSE);
  n.risk <- KM.ehat$n.risk;
  surv <- KM.ehat$surv;
  repeats <- c(diff( - n.risk), n.risk[length(n.risk)]);
  surv <- rep(surv, repeats);
  w <-  - diff(c(1, surv)); 
  m <- order(ehat,  - status);
  bla <- cumsum((w * ehat[m]));
  bla <- (bla[length(bla)] - bla)/(surv + state[m]);	## Put bla back into original order
  bl <- bla;
  bl[(1:N)[m]] <- bla;
  yhat <-  bl + ypred; #(X%*%beta)
  yy[state == 0] <- yhat[state == 0];
  y.imputed <- yy;
  return(y.imputed)
}


#14,BJSPP

DCA_BJ_lp <- function(X,y,delta,beta_c,la,p){
  X <- as.matrix(X);
  m <- dim(X)[1];
  n <- dim(X)[2];
  beta <- beta_c;
  w <- beta_c;
  z <- rep(0,m);
  I <- diag(rep(1,n));
  #la = 5
  rho1 = 1
  rho2 = 1
  miu <- rep(0,n);
  eta <- rep(0,m);
  l=0;
  #Tol <- ce(la)
  for (i in 1:5) {
    repeat {
      beta1 <- beta;
      v <- gradient(la,beta,p);
      beta <- SoftThreshold(w-(miu-v)/rho1,la/rho1);
      w <- solve(rho2 * t(X)%*%X + rho1 * I) %*% (rho2 * t(X) %*% (y + z) + rho1 * beta + miu - t(X) %*% eta );
      z <- (eta + rho2 * X %*% w - rho2 * y)/(1 + rho2);
      miu <- miu + rho1 * (beta - w);
      eta <- eta + rho2 * (X %*% w - y - z);
      
      
      l <- l + 1;
      beta2 <- beta;
      #print(norm_p(beta2-beta1,2)/norm_p(beta1, 2))
      if ((norm_p(beta2-beta1,2)/norm_p(beta1, 2) < 0.5) && (l>=100) ) {
        #print(l)  
        break  # 跳出循环
      }
    }
    y <- BJ(X,y,delta,beta);
    #print(i)
  }
  daic <- DAIC(delta,y,X%*%beta,beta)
  result1 = list(beta = beta,daic = daic);
  return(result1)
}

# B-spline
Bzk <- function(X){
  xx = bs(X[,1], df = 3)
  for (i in 2:dim(X)[2]) {
    m = bs(X[,i],df = 3)
    xx <- cbind(xx,m)
  }
  #xx <- as.matrix(x)
  return(xx)
}


# delta-AIC
co_k <- function(beta){
  co <- c(rep(0,5))
  s <- 3
  k <- c(0)
  for (i in seq(2,length(beta),s)) {
    three_sum <- sum(abs(beta[i:(i+2)]))
    k <- c(k, three_sum)
  }
  co[2:5] <- as.numeric(k[2:5]>0)
  co[1] <- sum(k>0)
  return(co)
}


DAIC <- function(delta,y,y_hat,beta){
  re <- c(0,0,0,rep(0,5))
  m <- sum(delta)
  re[1] <- sum(delta*(y-y_hat)^2)/m;
  re[4:8] <- co_k(beta)
  re[2] <- 2*re[4] + m* log(re[1])
  re[3] <- log(m)*re[4] + m * log(re[1])
  return(re)
}


Moni1 <- function(beta,n,sig,mu,U1){
  epsilon = rnorm(n, mean = 0, sd = 1);
  X = as.data.frame(mvrnorm(n, mu, sig, tol = 1e-6, empirical = FALSE, EISPACK = FALSE));
  x = as.matrix(X);
  Y = matrix(c(rep(0,n*2)),n,2);
  Y[,1] = cbind(sin(x[,1]),x[,2],x[,3]^2,x[,4],x[,-c(1,2,3,4)])%*%beta+epsilon;
  Ym <- Censor1(U1,Y,n);
  Ym <- as.data.frame(Ym);
  colnames(Ym) <- c("y","delta");
  result = data.frame(Ym,X);
  return(result)
}

Moni2 <- function(beta,n,sig,mu,U1){
  epsilon = rnorm(n, mean = 0, sd = 1);
  #epsilon2 <- rt(n,2)
  X = as.data.frame(mvrnorm(n, mu, sig, tol = 1e-6, empirical = FALSE, EISPACK = FALSE));
  x = as.matrix(X);
  Y = matrix(c(rep(0,n*2)),n,2);
  Y[,1] = cbind(sin(x[,1]),exp(x[,2])/(1+x[2]^2),x[,3]^2,x[,4],x[,-c(1,2,3,4)])%*%beta+epsilon;
  Ym <- Censor1(U1,Y,n);
  Ym <- as.data.frame(Ym);
  colnames(Ym) <- c("y","delta");
  result = data.frame(Ym,X);
  return(result)
}

Moni3 <- function(beta,n,sig,mu,U1){
  epsilon = rnorm(n, mean = 0, sd = 1);
  #epsilon2 <- rt(n,2)
  X = as.data.frame(rmvt(n, sig, df = 20));
  x = as.matrix(X);
  Y = matrix(c(rep(0,n*2)),n,2);
  Y[,1] = cbind(sin(x[,1]),exp(x[,2])/(1+x[2]^2),x[,3]^2,x[,4],x[,-c(1,2,3,4)])%*%beta+epsilon;
  Ym <- Censor1(U1,Y,n);
  Ym <- as.data.frame(Ym);
  colnames(Ym) <- c("y","delta");
  result = data.frame(Ym,X);
  return(result)
}

function(beta,n,sig,mu,U1){
  epsilon = rnorm(n, mean = 0, sd = 1);
  #epsilon2 <- rt(n,2)
  X = as.data.frame(rmvt(n, sig, df = 20));
  x = as.matrix(X);
  Y = matrix(c(rep(0,n*2)),n,2);
  Y[,1] = cbind(sin(x[,1]),x[,2],x[,3]^2,x[,4],x[,-c(1,2,3,4)])%*%beta+epsilon;
  Ym <- Censor1(U1,Y,n);
  Ym <- as.data.frame(Ym);
  colnames(Ym) <- c("y","delta");
  result = data.frame(Ym,X);
  return(result)
}


Bjsc <- function(X,y,delta){
  m <- cbind(rep(1,dim(X)[1]),X)
  for (i in 1:5) {
    
    fit <- cv.ncvreg(X, y, penalty = "SCAD")
    t = fit$lambda.min
    a = coef(fit, s = t)
    #3print(a)
    #y <- BJ(cbind(rep(1,dim(X)[1]),X),y,delta,a)
    y <- BJ(X,y,delta,a[2:length(a)])
  }
  
  DAICl <- DAIC(delta,y,m%*%a,a)
  #print(cbind(y-m%*%a,delta))
  result <- list(beta = a,DAICl = DAICl)
  print(t)
  return(result)
}

rank.sis.lrm <- function(d,ipw) {
  
  n <- length(ipw[,1])
  p=ncol(ipw[,-1])
  # T<-solve(t(X)%*%X)
  est<-matrix(0,p,1);
  
  library("doParallel")
  library("foreach")
  nprocs=detectCores() - 1;
  cl=makeCluster(nprocs);
  registerDoParallel(cl);
  est<-foreach(w=1:p,.combine=rbind,.packages = c("sfsmisc","MASS")) %dopar% {
    y=as.numeric(ipw[,1])
    X=ipw[,-1] 
    X1=X[,w]
    
    est[w]=cor(y,X1,method="spearman")
  }  #close for cycle w
  stopCluster(cl)
  M = rep(0,d)
  M[1] = xxt(est)
  I=sort(abs(est), decreasing = TRUE, index.return=TRUE)
  I1=I$ix[1:d]
  beta=est[I1]
  result=cbind(I1, beta,M)
  return(result)
}


Mic.SIS <-function(d,ipw){
  #ipwzr <- apply(ipw, 2, equals)
  est = mine(ipw[,-1], ipw[,1], alpha = 0.05)$MIC
  M = rep(0,d)
  M[1] = xxt(est)
  I=sort(abs(est), decreasing = TRUE, index.return=TRUE)
  I1=I$ix[1:d]
  beta=est[I1]
  result=cbind(I1, beta,M)
  return(result)
  
}

MrMic.SIS <-function(d,ipw){
  ipwzr <- apply(ipw, 2, equals)
  est = mine(ipwzr[,-1], ipwzr[,1], alpha = 0.05)$MIC
  M = rep(0,d)
  M[1] = xxt(est)
  I=sort(abs(est), decreasing = TRUE, index.return=TRUE)
  I1=I$ix[1:d]
  beta=est[I1]
  result=cbind(I1, beta, M)
  return(result)
  
}


lmdpd <- function(y,X,alpha,method = "L-BFGS-B",p=ncol(X), Initial=matrix(c(median(y),rep(1,p-1),log(mad(y))))){
  
  n <- length(y)
  X = as.matrix(X)
  T<-solve(t(X)%*%X)
  
  #Objective function for the computation of MDPDE (Minimum DPD estimators)
  objfunction_Linear_reg <- function(t) {
    beta<-t[1:p]
    sigma<-exp(t[p+1])
    r<- (-((y-X%*%beta)^2)/(2*sigma*sigma))
    c<-1/(sigma^alpha)
    
    if ( alpha == 0)
    { f = -mean(r)+t[p+1]}
    else
    { f =c*((1/sqrt(1+alpha))-((1+alpha)*mean(exp(r))/alpha))}
    return(f)
  }
  
  
  ## Obtain the MDPDE in 'est'
  ## Default method used in optim is "L-BFGS-B", in case of non-convergence, this method can be changed
  result<-optim(Initial, objfunction_Linear_reg , gr = NULL,
                method =  method,
                lower = -Inf, upper = Inf, control = list(), hessian = FALSE)
  est=result$par
  s<-exp(est[p+1])  ###sigma estimate
  est[p+1]<-s
  
  
  ## Record convergence indicator in 'conv'
  conv=result$convergence
  
  
  ## Output
  output<-list(est, conv)
  return(output)
}

dpd.sis <- function(d,ipw,alpha=0.5,Method="L-BFGS-B",p=ncol(X),Initial=matrix(rep(1,p))){
  
  n <- length(ipw[,1])
  p=ncol(ipw[,-1])
  est<-matrix(0,p,1);
  
  #n <- length(y)
  #p=ncol(X)
  
  #est<-matrix(0,p-1,1);
  
  library("doParallel")
  nprocs=detectCores() - 2
  cl=makeCluster(nprocs)
  registerDoParallel(cl)
  est<-foreach(w=2:p,.combine=rbind,.packages = c("sfsmisc","MASS"),.export = c("lmdpd")) %do% {
    y=as.numeric(ipw[,1])
    X = ipw[,-1]
    X1=X[,c(1,w)]
    init=matrix(c(median(y), Initial[w],log(mad(y))))
    
    ldpd<-lmdpd(y,X1,alpha,Method,2,init)
    paste(ldpd[[2]])
    betaw=ldpd[[1]]
    est[w]=betaw[2]
    #output for parallel processing Out<-est
  }#close for cycle w
  
  stopCluster(cl)
  
  M = rep(0,d)
  M[1] = xxt(est)
  I=sort(abs(est), decreasing = TRUE, index.return=TRUE)
  I1=I$ix[1:d]
  beta=est[I1]
  
  result=cbind(I1, beta,M)
  return(result)
}


dpd.sis <- function(d,ipw,alpha=0.5,Method="L-BFGS-B",p=ncol(X),Initial=matrix(rep(1,p))){
  
  n <- length(ipw[,1])
  p=ncol(ipw[,-1])
  est<-matrix(0,p,1);
  
  #n <- length(y)
  #p=ncol(X)
  
  #est<-matrix(0,p-1,1);
  
  library("doParallel")
  nprocs=detectCores() - 2
  cl=makeCluster(nprocs)
  registerDoParallel(cl)
  est<-foreach(w=2:p,.combine=rbind,.packages = c("sfsmisc","MASS"),.export = c("lmdpd")) %do% {
    y=as.numeric(ipw[,1])
    X = ipw[,-1]
    X1=X[,c(1,w)]
    init=matrix(c(median(y), Initial[w],log(mad(y))))
    
    ldpd<-lmdpd(y,X1,alpha,Method,2,init)
    paste(ldpd[[2]])
    betaw=ldpd[[1]]
    est[w]=betaw[2]
    #output for parallel processing Out<-est
  }#close for cycle w
  
  stopCluster(cl)
  
  M = rep(0,d)
  M[1] = xxt(est)
  I=sort(abs(est), decreasing = TRUE, index.return=TRUE)
  I1=I$ix[1:d]
  beta=est[I1]
  
  result=cbind(I1, beta,M)
  return(result)
}

#BJASS
BJASS <- function(X,Y,delta,beta_ini,k,S,beta_S){
  iter <- 1;
  maxiter <- 500;
  tol <- 1;
  X <- as.matrix(X)
  eig <- eigen(crossprod(X,X));
  u = max(eig$values);#找出最大特征值
  beta_old = beta_ini;
  n <- dim(X)[1];
  p <- dim(X)[2];
  gamma = c(rep(0,p));
  
  while ((iter <= maxiter) && (tol > 1E-3)) {  #2个条件，都满足就可以停止迭代
    e = Y - X[,S] %*% beta_S;
    
    Y <- BJ(X[,S],Y,delta,beta_old)
    s = Y - X %*% beta_old;#beta demension
    s = as.matrix(s)
    lold = - t(s)%*%s;
    gamma = beta_old + 1/u*(t(X)%*%s);
    tmp = sort(abs(gamma),decreasing = T);
    beta_new = gamma * (abs(gamma)>=tmp[k]);
    t <-as.matrix(Y-X%*%beta_new)
    
    
    lnew = - t(t) %*% t;
    if(lnew < lold){
      u = 2*u
    }
    tol = sqrt(sum((beta_new-beta_old)^2))
    #tol = norm(beta_new-beta_old);
    beta_old = beta_new;
    iter = iter + 1;
  }
  
  S = which(abs(beta_old) > 0);
  beta_S = beta_old[S];
  daic <- DAIC(delta,Y,X%*%beta_old,beta_old)
  result = list(S=S,beta_S=beta_S,daic=daic)
  return(result)
}

Bjsc <- function(X,y,delta){
  m <- cbind(rep(1,dim(X)[1]),X)
  for (i in 1:5) {
    
    fit <- cv.ncvreg(X, y, penalty = "SCAD")
    t = fit$lambda.min
    a = coef(fit, s = t)
    #3print(a)
    #y <- BJ(cbind(rep(1,dim(X)[1]),X),y,delta,a)
    y <- BJ(X,y,delta,a[2:length(a)])
  }
  
  DAICl <- DAIC(delta,y,m%*%a,a)
  #print(cbind(y-m%*%a,delta))
  result <- list(beta = a,DAICl = DAICl)
  print(t)
  return(result)
}





