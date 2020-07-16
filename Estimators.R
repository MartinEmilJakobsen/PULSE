################################################
##### Implementation of K-class estimators #####
################################################
K_class <- function(k,A,Z,Y,n){
  P_A <- A%*%solve(t(A)%*%A)%*%t(A)
  Wk <- t(Z) %*% ((1-k)*diag(n)+k* P_A )  %*%Z
    a <- solve(Wk) %*% t(Z) %*% ( (1-k)*diag(n) + k*P_A ) %*% Y 
  return(a)
}


K_class_lambda <- function(lambda,A,Z,Y,n){
  P_A <- A%*%solve(t(A)%*%A)%*%t(A)
  Wl <- t(Z) %*% (diag(n)+as.vector(lambda)* P_A )  %*%Z
    a <- solve(Wl) %*% t(Z) %*% ( diag(n) + as.vector(lambda)*P_A ) %*% Y 
  return(a)
}



LIML_k <- function(A,A_1,X,Y,n){
  P_A <- A%*%solve(t(A)%*%A)%*%t(A)
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


FULLER_k <- function(alpha,A,A_1,X,Y,n,dA){
  P_A <- A%*%solve(t(A)%*%A)%*%t(A)
  LIML_kappa <- LIML_k(A,A_1,X,Y,n)
  return(LIML_kappa-alpha/(n-dA))
}

########################################################################
###### Implementation of P-value approach with BinarySearch ############
########################################################################


Test_Statistic <- function(a,A,Z,Y,n){
  P_A <- A%*%solve(t(A)%*%A)%*%t(A)
  YtP_AY <- t(Y)%*%P_A%*%Y
  ZtP_AZ <- t(Z)%*%P_A%*%Z
  YtP_AZ <- t(Y)%*%P_A%*%Z
  YtY <- t(Y)%*%Y
  ZtZ <- t(Z)%*%Z
  YtZ <- t(Y)%*%Z
  return(n*(YtP_AY+t(a)%*%ZtP_AZ%*%a-2%*%YtP_AZ%*%a)/(YtY+t(a)%*%ZtZ%*%a-2%*%YtZ%*%a))
}


PULSE <- function(A,A_1,X,Y,p,N,n)
{
  if(A_1 == "none")
  {
    Z <- X
  } else {
    Z <- bind_cols(A_1,X)
  }
  dZ <- ncol(Z)
  dA <- ncol(A)
  n <- nrow(A)
  
  P_A <- A%*%solve(t(A)%*%A)%*%t(A)
  YtP_AY <- t(Y)%*%P_A%*%Y
  ZtP_AZ <- t(Z)%*%P_A%*%Z
  YtP_AZ <- t(Y)%*%P_A%*%Z
  YtY <- t(Y)%*%Y
  ZtZ <- t(Z)%*%Z
  YtZ <- t(Y)%*%Z
  
  q <- qchisq(1-p,df=dA, ncp=0,lower.tail = TRUE,log.p = FALSE)
  
  if(Test_Statistic(K_class(1,A,Z,Y,n),A,Z,Y,n)>= q){
    message("TSLS rejected")
    Fuller4Kappa <- FULLER_k(4,A,A_1,X,Y,n,dA)
    alpha <- K_class(Fuller4Kappa,A,Z,Y,n)
  } else if(Test_Statistic(K_class(0,A,Z,Y,n),A,Z,Y,n)<= q){
    message("OLS Accepted")
    alpha <- K_class(0,A,Z,Y,n)
  }
  else{
    
  lmax <- 2
  lmin <- 0
  
  while(Test_Statistic(K_class_lambda(lmax,A,Z,Y,n),A,Z,Y,n)>q){
    lmin <- lmax
    lmax <- lmax^2
  }
  
  Delta <- lmax-lmin
  
  while(Delta > 1/N || Accept == FALSE){
    l <- (lmin+lmax)/2
    alpha <- K_class_lambda(l,A,Z,Y,n)
    TestStatistic <- Test_Statistic(alpha,A,Z,Y,n)
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
  alpha <- K_class_lambda(lmax,A,Z,Y,n)
  }
  return(alpha)
}



