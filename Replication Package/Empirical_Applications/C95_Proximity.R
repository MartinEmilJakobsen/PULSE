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
library(ivpack)
library(kableExtra)



getwd()
source("../Estimators_Slow.R")

Proximity_Data <- read.delim(file="Data/C95_PROXIMITY.dat",sep="",header=FALSE)  %>% rename(
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
  filter(lwage76 != ".") %>% 
  mutate(lwage76 = as.numeric(lwage76)) %>% 
  mutate(iq = ifelse(iq==".","0",iq)) %>% 
  mutate(iq = as.numeric(iq)) %>% 
  mutate(kww = ifelse(kww==".","0",kww)) %>% 
  mutate(kww = as.integer(kww)) %>% 
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
####### Card (1995) Proximity to college - TABLE 3 and 4 ######## #################################################################
#################################################################
#################################################################

Selected_Data <- Proximity_Data 

Proximity_Data %>% 

############################################################################
# w/o family dummies
############################################################################
n <- nrow(Selected_Data)
#Target
Y <- Selected_Data %>% select(lwage76)  %>% as.matrix
#Included Exogenous:
A_1 <- Selected_Data %>% select(black,smsa76r,reg76r,smsa66r,reg662,reg663,reg664,reg665,reg666,reg667,reg668,reg669,Intercept) %>%  as.matrix
#All Exogenous
A <- Selected_Data %>% select(nearc4,age,agesq,black,smsa76r,reg76r,smsa66r,reg662,reg663,reg664,reg665,reg666,reg667,reg668,reg669,Intercept) %>% as.matrix
#Included Endogenous
X <- Selected_Data %>% select(ed76,exp,expsq) %>%  as.matrix
#Included (all)
Z <- cbind(X,A_1)



#OLS
m1.ols <- K_class(0,A,Z,Y,n)

#IV
m1.tsls <- K_class(1,A,Z,Y,n)

#PULSE
m1.pulse <- PULSE(A,A_1,X,Y,p=0.05,N=10000,n)
#fuller4
dA <- ncol(A)
m1.fuller4 <-  K_class(FULLER_k(4,A,A_1,X,Y,n,dA),A,Z,Y,n)


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
m2.ols <- K_class(0,A,Z,Y,n)

#IV
m2.tsls <- K_class(1,A,Z,Y,n)

#PULSE
m2.pulse <- PULSE(A,A_1,X,Y,p=0.05,N=10000,n)
#fuller4
dA <- ncol(A)
m2.fuller4 <-  K_class(FULLER_k(4,A,A_1,X,Y,n,dA),A,Z,Y,n)

#LM and IVREG check
lm(lwage76~ed76+exp+expsq+black+smsa76r+reg76r+smsa66r+reg662+reg663+reg664+reg665+reg666+reg667+reg668+reg669+daded+momed+nodaded+nomomed +f1+f2+f3+f4+f5+f6+f7+f8+momdad14+sinmom14,data=Selected_Data)
ivfit<-ivreg(lwage76~ed76+exp+expsq+black+smsa76r+reg76r+smsa66r+reg662+reg663+reg664+reg665+reg666+reg667+reg668+reg669+daded+momed+nodaded+nomomed +f1+f2+f3+f4+f5+f6+f7+f8+momdad14+sinmom14|nearc4+age+agesq+black+smsa76r+reg76r+smsa66r+reg662+reg663+reg664+reg665+reg666+reg667+reg668+reg669+daded+momed+nodaded+nomomed +f1+f2+f3+f4+f5+f6+f7+f8+momdad14+sinmom14,data=Selected_Data,x=TRUE)
summary(ivfit)

anderson.rubin.ci(ivfit)



###############
# with kww test score and IQ as its instrument
###############
Selected_Data <- Proximity_Data %>%  filter(kww!=0,iq!=0)

n <- nrow(Selected_Data)
#Target
Y <- Selected_Data %>% select(lwage76)  %>% as.matrix
#Included Exogenous:
A_1 <- Selected_Data %>% select(black,smsa76r,reg76r,smsa66r,reg662,reg663,reg664,reg665,reg666,reg667,reg668,reg669,Intercept,daded,momed,nodaded,nomomed ,f1,f2,f3,f4,f5,f6,f7,f8,momdad14,sinmom14) %>%  as.matrix
#All Exogenous
A <- Selected_Data %>% select(nearc4,age,agesq,iq,black,smsa76r,reg76r,smsa66r,reg662,reg663,reg664,reg665,reg666,reg667,reg668,reg669,Intercept,daded,momed,nodaded,nomomed ,f1,f2,f3,f4,f5,f6,f7,f8,momdad14,sinmom14) %>% as.matrix
#Included Endogenous
X <- Selected_Data %>% select(ed76,exp,expsq,kww) %>%  as.matrix
#Included (all)
Z <- cbind(X,A_1)

#OLS
m3.ols <- K_class(0,A,Z,Y,n)

#IV
m3.tsls <- K_class(1,A,Z,Y,n)

#PULSE
m3.pulse <- PULSE(A,A_1,X,Y,p=0.05,N=10000,n)

#fuller4
dA <- ncol(A)
m3.fuller4 <-  K_class(FULLER_k(4,A,A_1,X,Y,n,dA),A,Z,Y,n)


#############
Table <- data.frame(
  OLS = c(m1.ols["ed76",],
          m2.ols["ed76",],
          m3.ols["ed76",]),
  TSLS = c(m1.tsls["ed76",],
           m2.tsls["ed76",],
           m3.tsls["ed76",]),
  FULLER4 = c(m1.fuller4["ed76",],
              m2.fuller4["ed76",],
              m3.fuller4["ed76",]),
  PULSE = c(m1.pulse["ed76","lwage76"],
            m2.pulse["ed76","lwage76"],
            m3.pulse["ed76","lwage76"]),
  message = c(m1.pulse["ed76","m"],
              m2.pulse["ed76","m"],
              m3.pulse["ed76","m"]),
  test = c(m1.pulse["ed76","t"],
           m2.pulse["ed76","t"],
           m3.pulse["ed76","t"]),
  threshold = c(m1.pulse["ed76","q"],
                m2.pulse["ed76","q"],
                m3.pulse["ed76","q"]))

Table

kbl(Table,format="latex",digits=4)

