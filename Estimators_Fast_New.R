################################################
##### Implementation of K-class estimators #####
################################################
K_class <- function(k,Z,Y,n,P_A){
  Wk <- t(Z) %*% ((1-k)*diag(n)+k* P_A )  %*%Z
    a <- solve(Wk) %*% t(Z) %*% ( (1-k)*diag(n) + k*P_A ) %*% Y 
  return(a)
}


K_class_lambda <- function(lambda,Z,Y,n,P_A){
  Wl <- t(Z) %*% (diag(n)+as.vector(lambda)* P_A )  %*%Z
    a <- solve(Wl) %*% t(Z) %*% ( diag(n) + as.vector(lambda)*P_A ) %*% Y 
  return(a)
}



LIML_k <- function(Y,X,A_1,n,P_A){
  if(A_1 == "none"){
    M_A_1 <- diag(n)
  }
  else{
    M_A_1 <- diag(n)- A_1%*%solve(t(A_1)%*%(A_1))%*%t(A_1)
  }
  M_A = diag(n) - P_A
  
  W = t(cbind(Y,X))%*%M_A%*%cbind(Y,X)
  W_1 = t(cbind(Y,X))%*%M_A_1%*%cbind(Y,X)
  
  
  return(min(abs(eigen(W_1%*%solve(W))$values)))
}


FULLER_k <- function(alpha,Y,X,A_1,n,dA,P_A){
  LIML_kappa <- LIML_k(Y,X,A_1,n,P_A)
  return(LIML_kappa-alpha/(n-dA))
}

########################################################################
###### Implementation of P-value approach with BinarySearch ############
########################################################################


Test_Statistic <- function(a,YtP_AY,ZtP_AZ,YtP_AZ,YtY,ZtZ,YtZ,n){
  return(n*(YtP_AY+t(a)%*%ZtP_AZ%*%a-2%*%YtP_AZ%*%a)/(YtY+t(a)%*%ZtZ%*%a-2%*%YtZ%*%a))
}


BinarySearch <- function(Z,Y,dZ,dA,p,N,n,YtP_AY,ZtP_AZ,YtP_AZ,YtY,ZtZ,YtZ,P_A,A_1,X)
{
  q <- qchisq(1-p,df=dA, ncp=0,lower.tail = TRUE,log.p = FALSE)
  
  if(Test_Statistic(K_class(1,Z,Y,n,P_A),YtP_AY,ZtP_AZ,YtP_AZ,YtY,ZtZ,YtZ,n)>= q){
    Fuller4Kappa <- FULLER_k(4,Y,X,A_1=A_1,n,dA,P_A)
    alpha <- K_class(Fuller4Kappa,Z=Z,Y,n,P_A)
  }
  else if(Test_Statistic(K_class(0,Z,Y,n,P_A),YtP_AY,ZtP_AZ,YtP_AZ,YtY,ZtZ,YtZ,n)>= q){
    alpha <- K_class(0,Z=Z,Y,n,P_A)
    }
  else{
    
  lmax <- 2
  lmin <- 0
  
  while(Test_Statistic(K_class_lambda(lmax,Z,Y,n,P_A),YtP_AY,ZtP_AZ,YtP_AZ,YtY,ZtZ,YtZ,n)>q){
    lmin <- lmax
    lmax <- lmax^2
  }
  
  Delta <- lmax-lmin
  while(Delta > 1/N || Accept == FALSE){
    l <- (lmin+lmax)/2
    alpha <- K_class_lambda(l,Z,Y,n,P_A)
    TestStatistic <- Test_Statistic(alpha,YtP_AY,ZtP_AZ,YtP_AZ,YtY,ZtZ,YtZ,n)
    if(TestStatistic>q){
      Accept <- FALSE
      lmin <- l
      
    }
    else {
      Accept <- TRUE
      lmax <- l
    }
    Delta <- lmax-lmin
  }
  alpha <- K_class_lambda(lmax,Z,Y,n,P_A)
  }
  return(alpha)
}


