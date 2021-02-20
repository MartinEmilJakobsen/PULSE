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

setwd("/home/lnd974/causing/drafts/AR_and_Kclass/Simulations")
source("Estimators_Fast.R")


#Generating Coefficients
EXX <- 2
EXA <- 1
EAX <- EXA
EAA <- 1
EXY <- 2.5
EAY <- 1

PopTSLS <- (EXA*EAA^(-1)*EAX)^(-1)*EXA*(EAA)^(-1)*EAY
PopOLS <- (EXX)^(-1)*EXY 
PopK1o2 <- ((1-1/2)*EXX+(1/2)*EXA*EAA^(-1)*EAX)^(-1)*((1-1/2)*EXY+(1/2)*EXA*(EAA)^(-1)*EAY)
PopK3o4 <- ((1-3/4)*EXX+(3/4)*EXA*EAA^(-1)*EAX)^(-1)*((1-3/4)*EXY+(3/4)*EXA*(EAA)^(-1)*EAY)
PopK8o9 <- ((1-8/9)*EXX+(8/9)*EXA*EAA^(-1)*EAX)^(-1)*((1-8/9)*EXY+(8/9)*EXA*(EAA)^(-1)*EAY)



MSE <- function(coef,InstStrength){
  InstStrength^2 + 3+ InstStrength^2*coef^2+coef^2-2*InstStrength^2*coef-2*coef -2*0.5*coef
}

######### n=2000 #########


set.seed(1)

temp <- list()
for(i in 1: 50){
  n<- 2000
  A <- rnorm(n,mean=0,sd=1)
  U <- mvrnorm(n,mu=c(0,0),Sigma=matrix(c(1,0.5,0.5,1),ncol=2))
  X <-A + U[,1]
  Y <-X+ U[,2]
  Z <- X
  
  P_A <- A%*%solve(t(A)%*%A)%*%t(A)
  temp[i] <- list(data.frame(
    TSLS = K_class(1,Z=Z,Y,n,P_A) %>% as.numeric(),
    OLS = solve(t(Z)%*%Z)%*%t(Z)%*%Y %>% as.numeric(),
    K3o4 = K_class(3/4,Z=Z,Y,n,P_A) %>% as.numeric())
  )
}


IDs <- temp %>% bind_rows() %>% 
  mutate(.,ID=seq(1,nrow(.),1)) %>% 
  bind_rows(.,data.frame(TSLS = 1, OLS = 1.25,K3o4= 1.1, ID=9999)) %>% 
  gather(Type,Coef, c(-ID)) %>%  
  expand_grid(InstStrength=as.integer(seq(0,10000,1))) %>% 
  rowwise() %>% 
  mutate(MSE = MSE(Coef,InstStrength/100) ,
         InstStrength = abs(InstStrength)) %>%
  arrange(Type,Coef,ID,InstStrength) %>%
  ungroup %>% 
  group_by(Type,Coef,ID) %>% 
  mutate(SupMSE=cummax(MSE)) %>% 
  unique() %>% 
  ungroup %>% 
  select(-MSE,-Coef) %>%
  arrange(ID,InstStrength,SupMSE) %>% 
  ungroup %>% 
  group_by(ID,InstStrength) %>% 
  filter(SupMSE == min(SupMSE)) %>% 
  ungroup %>% 
  group_by(ID,Type) %>% 
  filter(InstStrength == min(InstStrength)|
           InstStrength == max(InstStrength)) %>% 
  mutate(Length = max(InstStrength)-min(InstStrength))

IDs %>% print(.,n=200)
IDs %>% filter(ID %in% c(13,3,9999)) %>% print



temp %>% 
  bind_rows() %>% 
  mutate(.,ID=seq(1,nrow(.),1)) %>% 
  bind_rows(.,data.frame(TSLS = 1, OLS = 1.25,K3o4= 1.1, ID=9999)) %>% 
  gather(Type,Coef,c(-ID)) %>% 
  expand_grid(InstStrength=seq(-600,600,1)) %>% 
  rowwise() %>% 
  mutate(MSE = MSE(Coef,InstStrength/100) ,
         InstStrength = abs(InstStrength)) %>%
  arrange(Type,Coef,InstStrength) %>%
  group_by(Type,Coef) %>% 
  mutate(SupMSE=cummax(MSE)) %>%
  group_by(Type,Coef,InstStrength) %>% 
  unique() %>% 
  ungroup %>%
  mutate(Type=ifelse(Type=="K3o4","kappa=3/4",Type)) %>% 
  {
  ggplot(data=.) +
  geom_line(aes(x=InstStrength/100,y=SupMSE,color=Type,group = Coef),alpha=0.2)+
  geom_line(data=. %>% filter(ID==13),aes(x=InstStrength/100,y=SupMSE,color=Type,group = Coef),size=0.8,alpha=1,linetype="dotted")+
  geom_line(data=. %>% filter(ID==3),aes(x=InstStrength/100,y=SupMSE,color=Type,group = Coef),size=0.8,alpha=1,linetype="dashed")+
  geom_line(data=. %>% filter(ID==9999),aes(x=InstStrength/100,y=SupMSE,group = Coef),color="black",size=0.5,alpha=0.9)+
  ylab("Worst case mean squared prediction error")+
  xlab("Maximum Intervention Strength sup|v|")+ 
  coord_cartesian(ylim=c(0.8, 1.2))+
  theme(plot.margin = unit(c(0,0.2,0,0), "cm"))+ 
  scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0))+
  theme(legend.position = "bottom")
}


ggsave("Plots/Dist_Robustness.png", plot = last_plot(), device = NULL, path = NULL,
       scale = 1, width = 12, height = 6, units = c("in"),
       dpi = 200, limitsize = FALSE)
ggsave("Plots/Dist_Robustness.eps", plot = last_plot(), device = cairo_ps, path = NULL,
       scale = 1, width = 12, height = 6, units = c("in"),
       dpi = 200, limitsize = FALSE)


