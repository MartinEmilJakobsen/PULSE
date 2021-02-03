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
library(AER)
library(matlib)
library(ivpack)




getwd()
source("../Estimators_Slow.R")

'Proximity_Data' <- read.delim(file="Data/nls.dat",sep="",header=FALSE)  %>% rename(
  id                =  V1,
  nearc2            =  V2,
    nearc4          =  V3,
    nearc4a         =  V4,
    nearc4b         =  V5,
    ed76            =  V6,
  ed66              =  V7,
  age76             =  V8,
    daded           =  V9,
  nodaded           =  V10,
    momed           =  V11,
  nomomed           =  V12,
    weight          =  V13,
    momdad14        =  V14,
    sinmom14        =  V15,
    step14          =  V16,
    reg661          =  V17,
    reg662          =  V18,
    reg663          =  V19,
    reg664          =  V20,
  reg665            =  V21,
  reg666            =  V22,
  reg667            =  V23,
  reg668            =  V24,
  reg669            =  V25,
    south66         =  V26,
    work76          =  V27,
    work78          =  V28,
    lwage76         =  V29,
    lwage78         =  V30,
    famed           =  V31,
  black             =  V32,
    smsa76r         =  V33,
    smsa78r         =  V34,
    reg76r          =  V35,
  reg78r            =  V36,
  reg80r            =  V37,
    smsa66r         =  V38,
    wage76          =  V39,
  wage78            =  V40,
  wage80            =  V41,
  noint78           =  V42,
    noint80         =  V43,
  enroll76          =  V44,
    enroll78        =  V45,
  enroll80          =  V46,
  kww               =  V47,
  iq                =  V48,
  marsta76          =  V49,
    marsta78        =  V50,
  marsta80          =  V51,
  libcrd14          =  V52
) %>% as_tibble() %>% 
  mutate(Intercept = 1)%>% 
  mutate(lwage76= ifelse(lwage76==".","0",lwage76) )%>%
  mutate(lwage76 = as.numeric(lwage76)) %>% 
  mutate(iq = ifelse(iq==".","0",iq)) %>% 
  mutate(iq = as.numeric(iq)) %>% 
  mutate(exp = age76-ed76-6) %>% 
  mutate(expsq = exp^2) %>% 
  mutate(age=age76,agesq=age76^2) %>% 
  mutate(f1 = ifelse(famed==1,1,0),
         f2 = ifelse(famed==2,1,0),
         f3 = ifelse(famed==3,1,0),
         f4 = ifelse(famed==4,1,0),
         f5 = ifelse(famed==5,1,0),
         f6 = ifelse(famed==6,1,0),
         f7 = ifelse(famed==7,1,0),
         f8 = ifelse(famed==8,1,0))
#################################################################
#################################################################
#################################################################
######## TABLE IV: MEN BORN 1930-1939 - 1980 CENSUS DATA ######## 
#################################################################
#################################################################
#################################################################

Selected_Data <- Proximity_Data 

############################################################################
# With family dummies
############################################################################
n <- nrow(Selected_Data)
#Target
Y <- Selected_Data %>% select(lwage76)  %>% as.matrix
#Included Exogenous:
A_1 <- Selected_Data %>% select(black,smsa76r,reg76r,smsa66r,reg662,reg663,reg664,reg665,reg666,reg667,reg668,reg669,Intercept,daded,momed,nodaded,nomomed ,f1,f2,f3,f4,f5,f6,f7,f8,momdad14,sinmom14) %>%  as.matrix
#All Exogenous
A <- Selected_Data %>% select(nearc4,age,agesq,black,smsa76r,reg76r,smsa66r,reg662,reg663,reg664,reg665,reg666,reg667,reg668,reg669,Intercept,daded,momed,nodaded,nomomed ,f1,f2,f3,f4,f5,f6,f7,f8,momdad14,sinmom14) %>% as.matrix
#Included Endogenous
X <- Selected_Data %>% select(ed76,exp,expsq) %>%  as.matrix
#Included (all)
Z <- cbind(X,A_1)

#OLS
K_class(0,A,Z,Y,n)

#IV
K_class(1,A,Z,Y,n)

#PULSE
PULSE(A,A_1,X,Y,p=0.05,N=10000,n)

#LM and IVREG check
lm(lwage76~ed76+exp+expsq+black+smsa76r+reg76r+smsa66r+reg662+reg663+reg664+reg665+reg666+reg667+reg668+reg669+daded+momed+nodaded+nomomed +f1+f2+f3+f4+f5+f6+f7+f8+momdad14+sinmom14,data=Selected_Data)
ivfit<-ivreg(logpgp95~avexpr|logem4,data=Selected_Data,x=TRUE)
summary(ivfit)

anderson.rubin.ci(ivfit)

############################################################################
# w/o family dummies
############################################################################
n <- nrow(Selected_Data)
#Target
Y <- Selected_Data %>% select(lwage76)  %>% as.matrix
#Included Exogenous:
A_1 <- Selected_Data %>% select(black,exp,expsq,smsa76r,reg76r,smsa66r,reg662,reg663,reg664,reg665,reg666,reg667,reg668,reg669,Intercept) %>%  as.matrix
#All Exogenous
A <- Selected_Data %>% select(nearc4,exp,expsq,black,smsa76r,reg76r,smsa66r,reg662,reg663,reg664,reg665,reg666,reg667,reg668,reg669,Intercept) %>% as.matrix
#Included Endogenous
X <- Selected_Data %>% select(ed76) %>%  as.matrix
#Included (all)
Z <- cbind(X,A_1)



#OLS
K_class(0,A,Z,Y,n)

#IV
K_class(1,A,Z,Y,n)

#PULSE
PULSE(A,A_1,X,Y,p=0.05,N=10000,n)

#LM and IVREG check
lm(lwage76~ed76+exp+expsq+black+smsa76r+reg76r+smsa66r+reg662+reg663+reg664+reg665+reg666+reg667+reg668+reg669,data=Selected_Data)
ivfit<-ivreg(lwage76~ed76+exp+expsq+black+smsa76r+reg76r+smsa66r+reg662+reg663+reg664+reg665+reg666+reg667+reg668+reg669|nearc4+age+agesq+black+smsa76r+reg76r+smsa66r+reg662+reg663+reg664+reg665+reg666+reg667+reg668+reg669+Intercept,data=Selected_Data,x=TRUE)
summary(ivfit)

anderson.rubin.ci(ivfit)