library(tidyverse)
library(magrittr)
library(matlib)
n <- 80000
alpha <- 3
delta1 <- 10
delta2 <- 5
gamma <- 3
A <- rnorm(n,mean=0,sd=1) %>% as.matrix
H <- rnorm(n,mean=0,sd=1) %>% as.matrix
X1 <- 2*A+ delta1*H + rnorm(n,mean=0,sd=1) %>% as.matrix
Y <- alpha*X1+delta2*H+rnorm(n,mean=0,sd=1) %>% as.matrix
X2 <- gamma*Y + rnorm(n,mean=0,sd=1) %>% as.matrix
X <- cbind(X1,X2)
Z <- X

p<- 0.05

#Estimators_Slow functions needed.

################################################
##### Implementation of K-class estimators #####
################################################
inv <- function(A){
  if(ncol(A)==1){
    r <- 1/A
  } else(
    r<- matlib::inv(A)
  )
  return(r)
}

K_class <- function(k,A,Z,Y,n){
  QA<- qr.Q(qr(A))
  RA<- qr.R(qr(A))
  temp <- solve((1-k)*t(Z)%*%Z + k* t(Z)%*%A%*%inv(RA)%*%t(QA)%*%Z) 
  a<-(1-k)*temp%*%t(Z) %*% Y +  k*temp%*%t(Z) %*%A%*%inv(RA)%*%t(QA)%*% Y
  return(a)
}

K_class_lambda <- function(lambda,A,Z,Y,n){
  QA<- qr.Q(qr(A))
  RA<- qr.R(qr(A))
  temp <- solve(t(Z)%*%Z + as.vector(lambda)* t(Z)%*%A%*%inv(RA)%*%t(QA)%*%Z) 
  a<-temp%*%t(Z) %*% Y +  as.vector(lambda)*temp%*%t(Z) %*%A%*%inv(RA)%*%t(QA)%*% Y
  return(a)
}

LIML_k <- function(A,A_1,X,Y,n){
  
  QA<- qr.Q(qr(A))
  RA<- qr.R(qr(A))
  if(norm(A_1) == 0){
    W_1 = t(cbind(Y,X))%*%cbind(Y,X)
  }
  else{
    #    M_A_1 <- diag(n)- A_1%*%solve(t(A_1)%*%(A_1))%*%t(A_1)
    W_1 = t(cbind(Y,X))%*%cbind(Y,X)  - t(cbind(Y,X))%*%A_1%*%solve(t(A_1)%*%(A_1))%*%t(A_1)%*%cbind(Y,X)
  }
  W = t(cbind(Y,X))%*%cbind(Y,X) - t(cbind(Y,X))%*%A%*%inv(RA)%*%t(QA)%*%cbind(Y,X)
  return(min(abs(eigen(W_1%*%solve(W))$values)))
}

FULLER_k <- function(alpha,A,A_1,X,Y,n,dA){
  LIML_kappa <- LIML_k(A,A_1,X,Y,n)
  return(LIML_kappa-alpha/(n-dA))
}

########################################################################
###### Implementation of P-value approach with BinarySearch ############
########################################################################

Test_Statistic <- function(a,A,Z,Y,n,q){
  QA<- qr.Q(qr(A))
  RA<- qr.R(qr(A))
  # YtP_AY <- t(Y)%*%A%*%inv(RA)%*%t(QA)%*%Y
  # ZtP_AZ <- t(Z)%*%A%*%inv(RA)%*%t(QA)%*%Z
  # YtP_AZ <- t(Y)%*%A%*%inv(RA)%*%t(QA)%*%Z
  # YtY <- t(Y)%*%Y
  # ZtZ <- t(Z)%*%Z
  # YtZ <- t(Y)%*%Z
  #(n-ncol(A)+q)*(YtP_AY+t(a)%*%ZtP_AZ%*%a-2%*%YtP_AZ%*%a)/(YtY+t(a)%*%ZtZ%*%a-2%*%YtZ%*%a)
  a<-(n-ncol(A)+q)*norm(A%*%(inv(RA)%*%t(QA)%*%Y)-A%*%(inv(RA)%*%t(QA)%*%Z%*%a),type="2")^2/norm(Y-Z%*%a,type="2")^2
  return(a)
}

PULSE <- function(A,A_1,X,Y,p,N,n)
{
  if(sum(A_1) == 0)
  {
    Z <- X
  } else {
    Z <- cbind(A_1,X)
  }
  dZ <- ncol(Z)
  dA <- ncol(A)
  n <- nrow(A)
  
  q <- qchisq(1-p,df=dA, ncp=0,lower.tail = TRUE,log.p = FALSE) 
  
   if(Test_Statistic(K_class(0,A,Z,Y,n),A,Z,Y,n,q)<= q){
    message(paste("OLS Accepted: T=",Test_Statistic(K_class(0,A,Z,Y,n),A,Z,Y,n,q),", q=",q))
    alpha <- K_class(0,A,Z,Y,n)
    m <- "OLS Accepted"
    t <- Test_Statistic(K_class(0,A,Z,Y,n),A,Z,Y,n,q)
    l <- 0 
    k <- 0
  }
  else{
    
    lmax <- 2
    lmin <- 0
    
    while(Test_Statistic(K_class_lambda(lmax,A,Z,Y,n),A,Z,Y,n,q)>q){
      lmin <- lmax
      lmax <- lmax^2
    }
    
    Delta <- lmax-lmin
    
    while(Delta > 1/N || Accept == FALSE){
      l <- (lmin+lmax)/2
      alpha <- K_class_lambda(l,A,Z,Y,n)
      TestStatistic <- Test_Statistic(alpha,A,Z,Y,n,q)
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
    m <- ""
    t<- Test_Statistic(alpha,A,Z,Y,n,q)
    l <- lmax
    k <- l/(1+l)
  }
  return(data.frame(alpha=alpha,m=m,t=t,q=q,l=l,k=k))
}


PULSE(A,A_1=0,X,Y,p,N=10000,n)



K_class_lambda(9999999,A,Z,Y,n)

(1-0.33196653*3)*3

