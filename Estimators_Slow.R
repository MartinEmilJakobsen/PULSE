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
  if(nrow(A_1) == 1){
    W_1 = t(solve(Y,X))%*%solve(Y,X)
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
  
  if(Test_Statistic(K_class(1,A,Z,Y,n),A,Z,Y,n,q)>= q){
    message(paste("TSLS rejected: T=",Test_Statistic(K_class(1,A,Z,Y,n),A,Z,Y,n,q),", q=",q))
    Fuller4Kappa <- FULLER_k(4,A,A_1,X,Y,n,dA)
    alpha <- K_class(Fuller4Kappa,A,Z,Y,n)
    m <- "TSLS rejected"
    t <- Test_Statistic(K_class(1,A,Z,Y,n),A,Z,Y,n,q)
  } else if(Test_Statistic(K_class(0,A,Z,Y,n),A,Z,Y,n,q)<= q){
    message(paste("OLS Accepted: T=",Test_Statistic(K_class(0,A,Z,Y,n),A,Z,Y,n,q),", q=",q))
    alpha <- K_class(0,A,Z,Y,n)
    m <- "OLS Accepted"
    t <- Test_Statistic(K_class(0,A,Z,Y,n),A,Z,Y,n,q)
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
  }
  return(data.frame(alpha=alpha,m=m,t=t,q=q))
}



