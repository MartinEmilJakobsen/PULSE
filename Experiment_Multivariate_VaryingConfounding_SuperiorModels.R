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


Models <- readRDS(file="Data/VaryingConfounding_MSESuperiorModelData_20200716080129.RDS")

Simulate <- function(nSim,n,delta11,delta12,delta21,delta22,mu11,mu22,RhoSq,VepX1,VepX2,xi11,xi12,xi21,xi22){
  
  xi = matrix(c(xi11,xi12,xi21,xi22),ncol=2,byrow=TRUE)
  delta = matrix(c(delta11,delta12,delta21,delta22),ncol=2,byrow=TRUE)
  mu = matrix(c(mu11,mu22),ncol=1,byrow=TRUE) 
  
  Coefs <- list(xi11=xi11,
                xi12=xi12,
                xi21=xi21,
                xi22=xi22,
                delta11=delta11,
                delta12=delta12,
                delta21=delta21,
                delta22=delta22,
                mu11 = mu11,
                mu22 = mu22,
                VepX1 = VepX1,
                VepX2 = VepX2,
                RhoSq = RhoSq)
  
  tempdat <- list()

  for(i in 1:nSim){
    Ep = mvrnorm(n,mu=rep(0,dA+2+dX+1),Sigma=diag(dA+2+dX+1))  
    
    A <- Ep[,1:2]
    H <- Ep[,3:4]
    X <- A%*%xi+ H%*%delta + Ep[,5:6]%*%matrix(c(VepX1,0,0,VepX2),ncol=2,byrow=TRUE)
    Y <- X%*%matrix(truealpha,ncol=1,byrow=TRUE)+ H%*% mu + Ep[,7]
    
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

nSim <- 25000
nObsPerSim <-  c(50)

message(paste("Initializing Parallization code @",as.character(Sys.time())))

plan(multiprocess, workers = 62)

dat <- Models %>% expand_grid(nSim = nSim, n=nObsPerSim) %>%
  mutate(Res = future_pmap(list(nSim=nSim, n=n,delta11=delta11,delta12=delta12,
                                delta21=delta21,delta22=delta22,mu11=mu11,
                                mu22=mu22,RhoSq=RhoSq,VepX1=VepX1,VepX2=VepX2,
                                xi11=xi11,xi12=xi12,xi21=xi21,xi22=xi22),.f=Simulate,.progress = TRUE)) %>%
  unnest(cols=c(Res))

message("Parallization code has been executed without error")

DataNameWOtimestamp <- paste0("Experiment_Multivariate_VaryingConfounding_SuperiorModels_nSim_",paste0(nSim,collapse="_"),"_nObsPerSim_",paste0(nObsPerSim,collapse="_"),"_")
saveRDS(dat,file=paste0("Data/",DataNameWOtimestamp,timestamp,".RDS"))

message(paste0("Data-save has been executed without error: ",DataNameWOtimestamp,timestamp,".RDS"))

file.copy(from="Experiment_Multivariate_VaryingConfounding_SuperiorModels.R", to=paste0("Data/Log/",DataNameWOtimestamp,timestamp,".R"), overwrite = TRUE, copy.mode = TRUE, copy.date = FALSE)

message("Log-save has been executed without error")
