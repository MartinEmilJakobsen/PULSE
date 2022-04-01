library(tidyverse)
library(magrittr)
library(matlib)
library(quadprog)
library(gridExtra)

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


#############################################################################################
###### Implementation of P-value approach with BinarySearch (for underid. setup) ############
#############################################################################################

Test_Statistic <- function(a,A,Z,Y,n,q){
  QA<- qr.Q(qr(A))
  RA<- qr.R(qr(A))
  a<-(n-ncol(A)+q)*norm(A%*%(inv(RA)%*%t(QA)%*%Y)-A%*%(inv(RA)%*%t(QA)%*%Z%*%a),type="2")^2/norm(Y-Z%*%a,type="2")^2
  return(a)
}


lIV <- function(a,A,Z,Y){
  QA<- qr.Q(qr(A))
  RA<- qr.R(qr(A))
  norm(A%*%(inv(RA)%*%t(QA)%*%Y)-A%*%(inv(RA)%*%t(QA)%*%Z%*%a),type="2")^2
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
    #message(paste("OLS Accepted: T=",Test_Statistic(K_class(0,A,Z,Y,n),A,Z,Y,n,q),", q=",q))
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



##########################################################################################

set.seed(1000)

N <- 7
alphas <- runif(N,1,2)
delta1s <- runif(N,1,2)
delta2s <- runif(N,1,2)
gammas  <- runif(N,1,2)
etas <- runif(N,0.1,1)
samplesize <- c(seq(50,400,50),1000,5000)

data.frame(alpha=alphas,delta1=delta1s,delta2=delta2s,gamma=gammas,eta=etas)

Data <- data.frame(alpha=alphas,delta1=delta1s,delta2=delta2s,gamma=gammas,eta=etas) %>% 
  mutate(model=1:n()) %>% 
  expand_grid(samplesize,
              rep=1:100) 

genest <- function(alpha,delta1,delta2,gamma,eta,samplesize){
  n <- samplesize
  A <- rnorm(n,mean=0,sd=1) %>% as.matrix
  H <- rnorm(n,mean=0,sd=1) %>% as.matrix
  X1 <- eta*A+ delta1*H + rnorm(n,mean=0,sd=1) %>% as.matrix
  Y <- alpha*X1+delta2*H+rnorm(n,mean=0,sd=1) %>% as.matrix
  X2 <- gamma*Y + rnorm(n,mean=0,sd=1) %>% as.matrix
  X <- cbind(X1,X2)
  Z <- X
  
  pulse <- PULSE(A,A_1=0,X,Y,p=0.05,N=100000000,n)
  bPULSE <- pulse$alpha
  PULSEqOLS <- ifelse(sum(pulse$l)==0,1,0)
  b2true <- (gamma*(delta2^2+1))/(1+gamma^2*(delta2^2+1))
  b1true <- (1-b2true*gamma)*alpha
  btrue <- c(b1true,b2true)
   
  Dmat <- 2*t(Z)%*%Z
  dvec <- 2*t(Y)%*%Z
  Amat <- t(A%*%solve(t(A)%*%A)%*%t(A)%*%Z)
  bvec <- A%*%solve(t(A)%*%A)%*%t(A)%*%Y
  Amat <- cbind(Amat,-Amat)
  bvec <- c(bvec-rep(0.0000000001,n),-bvec-rep(0.0000000001,n))
  bIV <- solve.QP(Dmat,dvec,Amat,bvec,meq=0,factorized=FALSE)$solution
  
  
  data.frame(TracePULSE  = sum((bPULSE - btrue)^2),
             TraceIV = sum((bIV - btrue)^2),
             PULSEqOLS = PULSEqOLS)
}


 
Dat <- Data %>% 
  mutate(res = pmap(.l=list(alpha,delta1,delta2,gamma,eta,samplesize),.f=genest)) %>% 
  unnest(cols=c('res'))


Results <- Dat %>% group_by(model,samplesize) %>% summarize(proc = mean(PULSEqOLS), TracePULSE=mean(TracePULSE),TraceIV=mean(TraceIV))
Results

Results %>% mutate(dec=(TraceIV-TracePULSE)/TraceIV) %>% group_by(samplesize) %>% summarize(meandec=mean(dec))

PlotData <- Results %>%  gather(Method,Trace,c(TracePULSE,TraceIV)) %>% 
  mutate(model=as.character(model))

p1 <- ggplot(data=PlotData,aes(x=samplesize,y=Trace,color=model,linetype=Method,group=interaction(model,Method)))+
  geom_line()+
  geom_point()+
  scale_x_continuous(trans="log")+
  scale_y_continuous(limits=c(0,0.1))

p2 <- ggplot(data=PlotData,aes(x=samplesize,y=proc,color=model,linetype=Method,group=interaction(model,Method)))+
  geom_line()+
  geom_point()

grid.arrange(p1, p2, nrow = 1)



p1 <- ggplot(data=PlotData %>% filter(Method=="TracePULSE"),aes(x=samplesize,y=Trace,group=interaction(model,Method)))+
  geom_line()+
  ylab("Trace MSE")+
  xlab("Sample size")+
  scale_x_continuous(trans="log10",breaks=c(50,100,300,1000,5000))

p1  
ggsave(paste0("Plots/UnderidentifiedConvergence.eps"), plot = last_plot(), device = "eps", path = NULL,
         scale = 1, width = 12, height = 4, units = c("in"),
         dpi = 200, limitsize = FALSE)
  



############### SMALL SAMPLESIZES MANY MODELS #################
set.seed(20000)

N <- 1000
alphas <- runif(N,1,2)
delta1s <- runif(N,1,2)
delta2s <- runif(N,1,2)
gammas  <- runif(N,1,2)
etas <- runif(N,0.1,1)
samplesize <- c(seq(50,300,50))

data.frame(alpha=alphas,delta1=delta1s,delta2=delta2s,gamma=gammas,eta=etas)

Data <- data.frame(alpha=alphas,delta1=delta1s,delta2=delta2s,gamma=gammas,eta=etas) %>% 
  mutate(model=1:n()) %>% 
  expand_grid(samplesize,
              rep=1:100) 

Dat <- Data %>% 
  mutate(res = pmap(.l=list(alpha,delta1,delta2,gamma,eta,samplesize),.f=genest)) %>% 
  unnest(cols=c('res'))


Results <- Dat %>% group_by(model,samplesize) %>% summarize(proc = mean(PULSEqOLS), TracePULSE=mean(TracePULSE),TraceIV=mean(TraceIV))
Results

Results %>% mutate(dec=(TraceIV-TracePULSE)/TraceIV) %>% group_by(samplesize) %>% summarize(meandec=mean(dec))





