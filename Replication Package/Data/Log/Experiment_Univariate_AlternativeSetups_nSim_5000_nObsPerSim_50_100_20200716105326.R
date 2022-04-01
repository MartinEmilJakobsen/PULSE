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


Simulate <- function(nSim,n,nInst,rho,Rsq,truealpha,xis=xis){
  dA <- nInst
  dX <- 1
  dZ <- 1
  xiEntry <- sqrt(dA*(1-Rsq)*Rsq)/(dA*(1-Rsq))

  if(xis== "Positive"){
    xi <- matrix(c(rep(xiEntry,dA)),ncol=1)
  }
  else if(xis == "Negative"){
    xi <- -matrix(c(rep(xiEntry,dA)),ncol=1)
  
  } else{
    xi <- diag(sample(c(-1,1),dA,replace=TRUE))%*%matrix(c(rep(xiEntry,dA)),ncol=1)
    
  }
    
  
  
  Coefs <- list(xiEntry=xiEntry,
                nInst=nInst,
                rho=rho,
                Rsq=Rsq,
                truealpha= truealpha) 
  
  tempdat <- list()
  
  for(i in 1:nSim)
  {
    
    U <- mvrnorm(n,mu=rep(0,2),Sigma=matrix(c(1,rho,rho,1),ncol=2,byrow=TRUE))
    A <- mvrnorm(n,mu=rep(0,dA),Sigma=diag(dA))
    X <- A%*%xi+  U[,1]
    Y <- X*truealpha+ U[,2]
    
    A_1= "none"
    Z <- X

    
    
    P_A <- A%*%solve(t(A)%*%A)%*%t(A)
    YtP_AY <- t(Y)%*%P_A%*%Y
    ZtP_AZ <- t(Z)%*%P_A%*%Z
    YtP_AZ <- t(Y)%*%P_A%*%Z
    YtY <- t(Y)%*%Y
    ZtZ <- t(Z)%*%Z
    YtZ <- t(Y)%*%Z
    
    
    Gn <- ((n-dA)/dA) * t(X)%*%P_A%*%X%*%solve(t(X)%*%(diag(n)-P_A)%*%X)
    
    tempdat[[i]] <- tibble( OLS = solve(t(Z)%*%Z)%*%t(Z)%*%Y,
                            TSLS = K_class(1,Z=Z,Y,n,P_A),
                            LIML = LIML_k(Y,X,A_1=A_1,n,P_A) %>% K_class(.,Z=Z,Y,n,P_A),
                            Ful1 = FULLER_k(1,Y,X,A_1=A_1,n,dA,P_A) %>% K_class(.,Z=Z,Y,n,P_A),
                            Ful4 = FULLER_k(4,Y,X,A_1=A_1,n,dA,P_A) %>% K_class(.,Z=Z,Y,n,P_A),
                            PULSE05 = BinarySearch(Z,Y,dZ,dA,p=0.05,N=10000,n=n,YtP_AY,ZtP_AZ,YtP_AZ,YtY,ZtZ,YtZ,P_A,A_1,X)
    )  %>%
      gather(Type,alpha)  %>%
      rowwise()  %>%
      mutate(MSESample = (alpha-truealpha)%*%t(alpha-truealpha),
             BiasSample = alpha- truealpha,
             Gn = Gn)  %>%
      ungroup()  
  }
  
  Data <- bind_rows(tempdat) %>% 
    ungroup() %>% 
    group_by(Type)  %>%
    summarise(MeanGn = mean(Gn),
              MSE = mean(MSESample),
              RMSE = sqrt(MSE),
              MeanEstimate = mean(alpha),
              VarEstimate = var(alpha),
              MeanBias = mean(BiasSample),
              IQRestimate = IQR(alpha),
              MedianBias = median(BiasSample),
              AbsMedianBias = abs(MedianBias)) %>%
    ungroup
  
  return(Data)
}

# %% Executing simulation

nSim <- 5000
nObsPerSim <-  c(50,100)
nInst <- c(5,10,30)
rho <- c(-0.9,-0.5,-0.3,-0.1,0.1,0.3,0.5,0.9)
Rsq <- c(0.0001,0.001,0.01,0.1,0.2,0.3)
truealpha <- c(-1,0,1)
xis <- c("Positive","Negative","Mixed")

message(paste("Initializing Parallization code @",as.character(Sys.time())))

plan(multiprocess, workers = 62)

dat <- expand_grid(nInst=nInst,rho=rho,Rsq=Rsq, nSim = nSim, n=nObsPerSim,truealpha= truealpha,xis=xis) %>%
  arrange(-nInst,-n) %>% 
  mutate(Res = future_pmap(list(nSim=nSim, n=n,nInst=nInst,rho=rho,Rsq=Rsq,truealpha= truealpha,xis=xis),.f=Simulate,.progress = TRUE)) %>%
  unnest(cols=c(Res))

message("Parallization code has been executed without error")

DataNameWOtimestamp <- paste0("Experiment_Univariate_AlternativeSetups_nSim_",paste0(nSim,collapse="_"),"_nObsPerSim_",paste0(nObsPerSim,collapse="_"),"_")
saveRDS(dat,file=paste0("Data/",DataNameWOtimestamp,timestamp,".RDS"))

message(paste0("Data-save has been executed without error: ",DataNameWOtimestamp,timestamp,".RDS"))

file.copy(from="Experiment_Univariate_AlternativeSetups.R", to=paste0("Data/Log/",DataNameWOtimestamp,timestamp,".R"), overwrite = TRUE, copy.mode = TRUE, copy.date = FALSE)

message("Log-save has been executed without error")
