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
library(kableExtra)

getwd()
source("../Estimators_Slow.R")

QOB_Data <- read.delim(file="Data/QOB",sep=" ",header=FALSE) %>% 
  as_tibble() %>% 
  rename(AGE=V1,
         AGEQ=V2,
         EDUC=V4,
         ENOCENT=V5,
         ESOCENT=V6,
         LWKLYWGE=V9,
         MARRIED=V10,
         MIDATL=V11,
         MT=V12,
         NEWENG=V13,
         CENSUS=V16,
         QOB=V18,
         RACE=V19,
         SMSA=V20,
         SOATL=V21,
         WNOCENT=V24,
         WSOCENT=V25,
         YOB=V27) %>% 
  select(-V8) %>% 
  mutate(COHORT = case_when( 20 <= YOB & YOB <=29 ~ "20.29",
                             30 <= YOB & YOB <=39 ~ "30.39",
                             40 <= YOB & YOB <=49 ~ "40.49",
                             30 <= YOB & YOB <=39 ~ "30.39",
                             TRUE ~ "NA"),
         AGEQ = ifelse(CENSUS==80,AGEQ-1900,AGEQ),
         AGEQSQ = AGEQ^2,
         YOBTEMP = YOB-floor(YOB/10)*10)  %>%  #CHECK IF FAILURE
  mutate(YR0 = ifelse(YOBTEMP==0,1,0),
         YR1 = ifelse(YOBTEMP==1,1,0),
         YR2 = ifelse(YOBTEMP==2,1,0),
         YR3 = ifelse(YOBTEMP==3,1,0),
         YR4 = ifelse(YOBTEMP==4,1,0),
         YR5 = ifelse(YOBTEMP==5,1,0),
         YR6 = ifelse(YOBTEMP==6,1,0),
         YR7 = ifelse(YOBTEMP==7,1,0),
         YR8 = ifelse(YOBTEMP==8,1,0),
         YR9 = ifelse(YOBTEMP==9,1,0),
         QTR1 = ifelse(QOB == 1, 1, 0),
         QTR2 = ifelse(QOB == 2, 1, 0),
         QTR3 = ifelse(QOB == 3, 1, 0),
         QTR4 = ifelse(QOB == 4, 1, 0),
         QTRxYR10= QTR1*YR0,
         QTRxYR11= QTR1*YR1,
         QTRxYR12= QTR1*YR2,
         QTRxYR13= QTR1*YR3,
         QTRxYR14= QTR1*YR4,
         QTRxYR15= QTR1*YR5,
         QTRxYR16= QTR1*YR6,
         QTRxYR17= QTR1*YR7,
         QTRxYR18= QTR1*YR8,
         QTRxYR19= QTR1*YR9,
         QTRxYR20= QTR2*YR0,
         QTRxYR21= QTR2*YR1,
         QTRxYR22= QTR2*YR2,
         QTRxYR23= QTR2*YR3,
         QTRxYR24= QTR2*YR4,
         QTRxYR25= QTR2*YR5,
         QTRxYR26= QTR2*YR6,
         QTRxYR27= QTR2*YR7,
         QTRxYR28= QTR2*YR8,
         QTRxYR29= QTR2*YR9,
         QTRxYR30= QTR3*YR0,
         QTRxYR31= QTR3*YR1,
         QTRxYR32= QTR3*YR2,
         QTRxYR33= QTR3*YR3,
         QTRxYR34= QTR3*YR4,
         QTRxYR35= QTR3*YR5,
         QTRxYR36= QTR3*YR6,
         QTRxYR37= QTR3*YR7,
         QTRxYR38= QTR3*YR8,
         QTRxYR39= QTR3*YR9) %>% 
  mutate(Intercept = 1)
#################################################################
#################################################################
#################################################################
######## TABLE IV: MEN BORN 1930-1939 - 1980 CENSUS DATA ######## 
#################################################################
#################################################################
#################################################################

Selected_Data <- QOB_Data  %>% filter(COHORT == "30.39") %>% 
  select(-c(V3,V7,V14,V15,V17,V22,V23)) 

############################################################################
#COLUMN 1 & 2 : LWKLYWGE ~ YR20-YR29 + EDUC (INSTRUMENTS = YR20-YR29 + QTR)#
############################################################################
n <- nrow(Selected_Data)
#Target
Y <- Selected_Data %>% select(LWKLYWGE) %>%  as.matrix
#Included Exogenous:
A_1 <- Selected_Data %>% select(starts_with("YR")) %>% as.matrix
#All Exogenous
A <- Selected_Data %>% select(starts_with("YR"),starts_with("QTRxYR")) %>% as.matrix
#Included Endogenous
X <- Selected_Data %>% select(EDUC) %>%  as.matrix
#Included (all)
Z <- cbind(X,A_1)

#OLS
m1.ols <- K_class(0,A,Z,Y,n)

#IV
m1.tsls <- K_class(1,A,Z,Y,n)

#PULSE
m1.pulse <- PULSE(A,A_1,X,Y,p=0.05,N=10000,n)

#LM and IVREG check
lm(LWKLYWGE~EDUC+YR0+YR1+YR2+YR3+YR4+YR5+YR6+YR7+YR8+YR9-1,data=Selected_Data)
ivfit<-ivreg(LWKLYWGE~EDUC+YR0+YR1+YR2+YR3+YR4+YR5+YR6+YR7+YR8+YR9-1|YR0+YR1+YR2+YR3+YR4+YR5+YR6+YR7+YR8+YR9+
QTRxYR10+QTRxYR11+QTRxYR12+QTRxYR13+QTRxYR14+QTRxYR15+QTRxYR16+QTRxYR17+QTRxYR18+QTRxYR19+QTRxYR20+QTRxYR21+QTRxYR22+QTRxYR23+QTRxYR24+QTRxYR25+QTRxYR26+QTRxYR27+QTRxYR28+QTRxYR29+QTRxYR30+QTRxYR31+QTRxYR32+QTRxYR33+QTRxYR34+QTRxYR35+QTRxYR36+QTRxYR37+QTRxYR38+QTRxYR39,data=Selected_Data,x=TRUE)
summary(ivfit)

anderson.rubin.ci(ivfit)


############################################################################
#COLUMN 3 & 4 : LWKLYWGE ~ YR20-YR29 + AGEQ + AGEQSQ + EDUC (INSTRUMENTS = YR20-YR29 + QTR + AGEQ + AGEQSQ)#
############################################################################
n <- nrow(Selected_Data)
#Target
Y <- Selected_Data %>% select(LWKLYWGE) %>%  as.matrix
#Included Exogenous:
A_1 <- Selected_Data %>% select(starts_with("YR"),AGEQ,AGEQSQ) %>% as.matrix
#All Exogenous
A <- Selected_Data %>% select(starts_with("YR"),AGEQ,AGEQSQ,starts_with("QTRxYR"),-QTRxYR38,-QTRxYR39) %>% as.matrix
#Included Endogenous
X <- Selected_Data %>% select(EDUC) %>%  as.matrix
#Included (all)
Z <- cbind(X,A_1)

#OLS
m2.ols <-K_class(0,A,Z,Y,n)
#IV
m2.tsls <- K_class(1,A,Z,Y,n)

#PULSE
m2.pulse <- PULSE(A,A_1,X,Y,p=0.05,N=1000,n)

#LM and IVREG check

############################################################################
#COLUMN 5 & 6 : LWKLYWGE ~ YR20-YR29 + RACE MARRIED SMSA NEWENG MIDATL ENOCENT WNOCENT SOATL ESOCENT WSOCENT MT + EDUC (INSTRUMENTS = YR20-YR29 + QTR + RACE MARRIED SMSA NEWENG MIDATL ENOCENT WNOCENT SOATL ESOCENT WSOCENT MT)#
############################################################################

n <- nrow(Selected_Data)
#Target
Y <- Selected_Data %>% select(LWKLYWGE) %>%  as.matrix
#Included Exogenous:
A_1 <- Selected_Data %>% select(starts_with("YR"),RACE,MARRIED,SMSA,NEWENG,MIDATL,ENOCENT,WNOCENT,SOATL,ESOCENT,WSOCENT,MT) %>% as.matrix
#All Exogenous
A <- Selected_Data %>% select(starts_with("YR"),RACE,MARRIED,SMSA,NEWENG,MIDATL,ENOCENT,WNOCENT,SOATL,ESOCENT,WSOCENT,MT,starts_with("QTRxYR")) %>% as.matrix
#Included Endogenous
X <- Selected_Data %>% select(EDUC) %>%  as.matrix
#Included (all)
Z <- cbind(X,A_1)

#OLS
m3.ols <- K_class(0,A,Z,Y,n)
#IV
m3.tsls <- K_class(1,A,Z,Y,n)

#PULSE
m3.pulse <- PULSE(A,A_1,X,Y,p=0.05,N=1000,n)

#LM and IVREG check
lm(LWKLYWGE~EDUC+RACE+MARRIED+SMSA+NEWENG+MIDATL+ENOCENT+WNOCENT+SOATL+ESOCENT+WSOCENT+MT+YR0+YR1+YR2+YR3+YR4+YR5+YR6+YR7+YR8+YR9-1,data=Selected_Data)
ivfit<-ivreg(LWKLYWGE~EDUC+RACE+MARRIED+SMSA+NEWENG+MIDATL+ENOCENT+WNOCENT+SOATL+ESOCENT+WSOCENT+MT+YR0+YR1+YR2+YR3+YR4+YR5+YR6+YR7+YR8+YR9-1|RACE+MARRIED+SMSA+NEWENG+MIDATL+ENOCENT+WNOCENT+SOATL+ESOCENT+WSOCENT+MT+YR0+YR1+YR2+YR3+YR4+YR5+YR6+YR7+YR8+YR9+QTRxYR10+QTRxYR11+QTRxYR12+QTRxYR13+QTRxYR14+QTRxYR15+QTRxYR16+QTRxYR17+QTRxYR18+QTRxYR19+QTRxYR20+QTRxYR21+QTRxYR22+QTRxYR23+QTRxYR24+QTRxYR25+QTRxYR26+QTRxYR27+QTRxYR28+QTRxYR29+QTRxYR30+QTRxYR31+QTRxYR32+QTRxYR33+QTRxYR34+QTRxYR35+QTRxYR36+QTRxYR37+QTRxYR38+QTRxYR39,data=Selected_Data,x=TRUE)
summary(ivfit)
anderson.rubin.ci(ivfit)

############################################################################
#COLUMN 7 & 8 : LWKLYWGE ~ YR20-YR29 + AGEQ + AGEQSQ RACE + MARRIED SMSA NEWENG MIDATL ENOCENT WNOCENT SOATL ESOCENT WSOCENT MT + EDUC (INSTRUMENTS = YR20-YR29 + QTR +  AGEQ + AGEQSQ RACE MARRIED SMSA NEWENG MIDATL ENOCENT WNOCENT SOATL ESOCENT WSOCENT MT)#
############################################################################

n <- nrow(Selected_Data)
#Target
Y <- Selected_Data %>% select(LWKLYWGE) %>%  as.matrix
#Included Exogenous:
A_1 <- Selected_Data %>% select(starts_with("YR"),AGEQ,AGEQSQ,RACE, MARRIED, SMSA, NEWENG, MIDATL, ENOCENT, WNOCENT, SOATL, ESOCENT, WSOCENT, MT) %>% as.matrix
#All Exogenous
A <- Selected_Data %>% select(starts_with("YR"),AGEQ,AGEQSQ,RACE, MARRIED, SMSA, NEWENG, MIDATL, ENOCENT, WNOCENT, SOATL, ESOCENT, WSOCENT, MT,starts_with("QTRxYR"),-QTRxYR38,-QTRxYR39) %>% as.matrix
#Included Endogenous
X <- Selected_Data %>% select(EDUC) %>%  as.matrix
#Included (all)
Z <- cbind(X,A_1)

#OLS
m4.ols <- K_class(0,A,Z,Y,n)
#IV
m4.tsls <- K_class(1,A,Z,Y,n)
#PULSE
m4.pulse <- PULSE(A,A_1,X,Y,p=0.05,N=10000,n)


lm(LWKLYWGE~EDUC+
RACE+MARRIED+SMSA+NEWENG+MIDATL+ENOCENT+WNOCENT+SOATL+ESOCENT+WSOCENT+MT+
YR0+YR1+YR2+YR3+YR4+YR5+YR6+YR7+YR8+YR9-1,data=Selected_Data)
ivfit<-ivreg(LWKLYWGE~EDUC+AGEQ+AGEQSQ+RACE+MARRIED+SMSA+NEWENG+MIDATL+ENOCENT+WNOCENT+SOATL+ESOCENT+WSOCENT+MT+
YR0+YR1+YR2+YR3+YR4+YR5+YR6+YR7+YR8+YR9-1|AGEQ+AGEQSQ+RACE+MARRIED+SMSA+NEWENG+MIDATL+ENOCENT+WNOCENT+SOATL+ESOCENT+WSOCENT+MT+YR0+YR1+YR2+YR3+YR4+YR5+YR6+YR7+YR8+YR9+QTRxYR10+QTRxYR11+QTRxYR12+QTRxYR13+QTRxYR14+QTRxYR15+QTRxYR16+QTRxYR17+QTRxYR18+QTRxYR19+QTRxYR20+QTRxYR21+QTRxYR22+QTRxYR23+QTRxYR24+QTRxYR25+QTRxYR26+QTRxYR27+QTRxYR28+QTRxYR29+QTRxYR30+QTRxYR31+QTRxYR32+QTRxYR33+QTRxYR34+QTRxYR35+QTRxYR36+QTRxYR37+QTRxYR38+QTRxYR39,data=Selected_Data,x=TRUE)
summary(ivfit)
anderson.rubin.ci(ivfit)




############################
###### Custom model ########
############################
Selected_Data_Centered <- Selected_Data# %>% mutate_at(c("LWKLYWGE","AGEQ","AGEQSQ","EDUC"),.funs=function(x){x-mean(x)})

n <- nrow(Selected_Data_Centered)
#Target
Y <- Selected_Data_Centered %>% select(LWKLYWGE) %>%  as.matrix
#Included Exogenous:
A_1 <- Selected_Data_Centered %>% select(AGEQ,AGEQSQ) %>% as.matrix
#All Exogenous
A <- Selected_Data_Centered %>% select(AGEQ,AGEQSQ,QTR1,QTR2,QTR3) %>% as.matrix
#Included Endogenous
X <- Selected_Data_Centered %>% select(EDUC) %>%  as.matrix
#Included (all)
Z <- cbind(X,A_1)

lm(X~A)

#OLS
m5.ols <- K_class(0,A,Z,Y,n)
#IV
m5.tsls <- K_class(1,A,Z,Y,n)
#PULSE
m5.pulse <- PULSE(A,A_1,X,Y,p=0.05,N=10000,n)

Table <- data.frame(
  OLS = c(m1.ols["EDUC",],
          m2.ols["EDUC",],
          m3.ols["EDUC",],
          m4.ols["EDUC",],
          m5.ols["EDUC",]),
  TSLS = c(m1.tsls["EDUC",],
           m2.tsls["EDUC",],
           m3.tsls["EDUC",],
           m4.tsls["EDUC",],
           m5.tsls["EDUC",]),
  PULSE = c(m1.pulse["EDUC","LWKLYWGE"],
            m2.pulse["EDUC","LWKLYWGE"],
            m3.pulse["EDUC","LWKLYWGE"],
            m4.pulse["EDUC","LWKLYWGE"],
            m5.pulse["EDUC","LWKLYWGE"]),
  message = c(m1.pulse["EDUC","m"],
              m2.pulse["EDUC","m"],
              m3.pulse["EDUC","m"],
              m4.pulse["EDUC","m"],
              m5.pulse["EDUC","m"]),
  test = c(m1.pulse["EDUC","t"],
           m2.pulse["EDUC","t"],
           m3.pulse["EDUC","t"],
           m4.pulse["EDUC","t"],
           m5.pulse["EDUC","t"]),
threshold = c(m1.pulse["EDUC","q"],
           m2.pulse["EDUC","q"],
           m3.pulse["EDUC","q"],
           m4.pulse["EDUC","q"],
           m5.pulse["EDUC","q"]))

kbl(Table,format="latex",digits=4)
