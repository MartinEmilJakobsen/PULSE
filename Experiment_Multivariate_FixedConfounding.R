library(MASS) # mvrnorm
library(matrixcalc)
library(stats)
library(lqmm) #makepositive.
library(expm) #sqrtm A=SS
library(furrr)
library(tidyverse)
library(metR)
library(gridExtra)
library(rlang)
library(stringr)

setwd("/home/lnd974/PULSE")
source("Estimators_Fast.R")

timestamp <- as.character(Sys.time()) %>% {str_replace_all(.,"[: -]","")}

set.seed(1)

dA <- 2
dX <- 2
dZ <- 2

truealpha <- c(0,0)
dimalpha <- ifelse(dZ == length(truealpha), dZ, NA)

CovStructure <- data.frame(
  phi1 = c(0.75,  0.353553, 0.424264,  0.9),
  phi2 = c(0.75,  0.353553, 0.424264,  0.9),
  eta =  c(0.25,  0.25,     0.8,       0.8)
) %>%  mutate(RhoNorm = sqrt((phi1^2+phi2^2-2*eta*phi1*phi2)/(1-eta^2)))


Simulate <- function(nSim,n,Cov){

  phi1 <- CovStructure[Cov,1]
  phi2 <- CovStructure[Cov,2]
  eta <-  CovStructure[Cov,3]
  
  
  xi11 <-runif(1,min=-2,max=2)
  xi12 <-runif(1,min=-2,max=2)
  xi21 <-runif(1,min=-2,max=2)
  xi22 <-runif(1,min=-2,max=2)
  
  xi = matrix(c(xi11,xi12,xi21,xi22),ncol=2,byrow=TRUE)
  
  Coefs <- list(xi11=xi11,xi12=xi12,xi21=xi21,xi22=xi22,phi1=phi1,phi2=phi2,eta=eta)
  
  tempdat <- list()
  
  for(i in 1:nSim){
    U_XY = mvrnorm(n,mu=c(0,0,0),Sigma=matrix(c(1,   eta,   phi1,
                                                eta,   1,   phi2,
                                                phi1, phi2, 1   ),ncol=3,byrow=TRUE))  
    
    A <- mvrnorm(n,mu=c(0,0),Sigma=matrix(c(1,0,0,1),ncol=2,byrow=TRUE))
    X <- A%*%xi+ U_XY[,1:2]
    Y <- X%*%matrix(truealpha,ncol=1,byrow=TRUE) + U_XY[,3]
    
    A_1 <- "none"
    Z <- X
    
    P_A <- A%*%solve(t(A)%*%A)%*%t(A)
    YtP_AY <- t(Y)%*%P_A%*%Y
    ZtP_AZ <- t(Z)%*%P_A%*%Z
    YtP_AZ <- t(Y)%*%P_A%*%Z
    YtY <- t(Y)%*%Y
    ZtZ <- t(Z)%*%Z
    YtZ <- t(Y)%*%Z
    
    Gn <- ((n-dA)/dA) * t(X)%*%P_A%*%X%*%solve(t(X)%*%(diag(n)-P_A)%*%X)
    
    tryCatch(  
      tempdat[[i]] <- tibble( OLS = list(solve(t(Z)%*%Z)%*%t(Z)%*%Y),
                              Ful1 = list(FULLER_k(1,Y,X,A_1=A_1,n,dA,P_A) %>% K_class(.,Z=Z,Y,n,P_A)),
                              Ful4 = list(FULLER_k(4,Y,X,A_1=A_1,n,dA,P_A) %>% K_class(.,Z=Z,Y,n,P_A)),
                              PULSE05 = list(BinarySearch(Z,Y,dZ,dA,p=0.05,N=10000,n=n,YtP_AY,ZtP_AZ,YtP_AZ,YtY,ZtZ,YtZ,P_A,A_1,X))
      )  %>%
        gather(Type,alpha)  %>%
        rowwise()  %>%
        mutate(MSESample = list((alpha-truealpha)%*%t(alpha-truealpha)),
               BiasSample = list(alpha- truealpha),
               Gn = list(Gn))  %>%
        ungroup()  , 
      error = function(e) {
        return(as.numeric(NA))
      } ) 
  }

  
  Data <- bind_rows(tempdat) %>% 
    ungroup() %>% 
    group_by(Type)  %>%
    summarise(MeanGn = list(rowMeans(array(unlist(Gn),dim=c(dX,dX,nrow(.))),dims=2)),
              MSE = list(rowMeans(array(unlist(MSESample),dim=c(dimalpha,dimalpha,nrow(.))),dims=2)),
              MeanEstimate = list(rowMeans(array(unlist(alpha),dim=c(dimalpha,nrow(.))),dims=1)),
              Bias = list(rowMeans(array(unlist(BiasSample),dim=c(dimalpha,nrow(.))),dims=1))) %>%
    mutate(Determinant = map_dbl(.x= MSE,.f= function(x) {det(x)}),
           Trace = map_dbl(.x= MSE,.f= function(x) {sum(diag(x))}),
           Bias = map_dbl(.x= Bias,.f= function(x) {sqrt(sum(x^2))}),
           ModelCoefs = list(Coefs)) %>% 
    ungroup
  
  return(Data)
}

# %% Executing simulation

nSim <- 5000
nObsPerSim <-  c(50)
Cov <- c(1,2,3,4)
nModel <-  seq(1,5000,1)

message(paste("Initializing Parallization code @",as.character(Sys.time())))

plan(multiprocess, workers = 62)

dat <- expand_grid(nModel=nModel, nSim = nSim, n=nObsPerSim, Cov = Cov) %>%
  mutate(Res = future_pmap(list(nSim=nSim, n=n, Cov=Cov),.f=Simulate,.progress = TRUE)) %>%
  unnest(cols=c(Res))

message("Parallization code has been executed without error")

DataNameWOtimestamp <- paste0("Experiment_Multivariate_FixedConfounding_nSim_",paste0(nSim,collapse="_"),"_nObsPerSim_",paste0(nObsPerSim,collapse="_"),"_nModel_",length(nModel),"_")
saveRDS(dat,file=paste0("Data/",DataNameWOtimestamp,timestamp,".RDS"))

message(paste0("Data-save has been executed without error: ",DataNameWOtimestamp,timestamp,".RDS"))

file.copy(from="Experiment_Multivariate_FixedConfounding.R", to=paste0("Data/Log/",DataNameWOtimestamp,timestamp,".R"), overwrite = TRUE, copy.mode = TRUE, copy.date = FALSE)

message("Log-save has been executed without error")
