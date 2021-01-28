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


getwd()
source("../Estimators_Fast.R")

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
  mutate(YR20 = ifelse(YOBTEMP==0,1,0),
         YR21 = ifelse(YOBTEMP==1,1,0),
         YR22 = ifelse(YOBTEMP==2,1,0),
         YR23 = ifelse(YOBTEMP==3,1,0),
         YR24 = ifelse(YOBTEMP==4,1,0),
         YR25 = ifelse(YOBTEMP==5,1,0),
         YR26 = ifelse(YOBTEMP==6,1,0),
         YR27 = ifelse(YOBTEMP==7,1,0),
         YR28 = ifelse(YOBTEMP==8,1,0),
         YR29 = ifelse(YOBTEMP==9,1,0),
         QTR1 = ifelse(QOB == 1, 1, 0),
         QTR2 = ifelse(QOB == 2, 1, 0),
         QTR3 = ifelse(QOB == 3, 1, 0),
         QTR4 = ifelse(QOB == 4, 1, 0),
         QTR120= QTR1*YR20,
         QTR121= QTR1*YR21,
         QTR122= QTR1*YR22,
         QTR123= QTR1*YR23,
         QTR124= QTR1*YR24,
         QTR125= QTR1*YR25,
         QTR126= QTR1*YR26,
         QTR127= QTR1*YR27,
         QTR128= QTR1*YR28,
         QTR129= QTR1*YR29,
         QTR220= QTR2*YR20,
         QTR221= QTR2*YR21,
         QTR222= QTR2*YR22,
         QTR223= QTR2*YR23,
         QTR224= QTR2*YR24,
         QTR225= QTR2*YR25,
         QTR226= QTR2*YR26,
         QTR227= QTR2*YR27,
         QTR228= QTR2*YR28,
         QTR229= QTR2*YR29,
         QTR320= QTR3*YR20,
         QTR321= QTR3*YR21,
         QTR322= QTR3*YR22,
         QTR323= QTR3*YR23,
         QTR324= QTR3*YR24,
         QTR325= QTR3*YR25,
         QTR326= QTR3*YR26,
         QTR327= QTR3*YR27,
         QTR328= QTR3*YR28,
         QTR329= QTR3*YR29) 



Selected_Data <- QOB_Data  %>% filter(COHORT == "30.39") %>% 
  select(-c(V3,V7,V14,V15,V17,V22,V23))

Selected_Data_OLS <- Selected_Data %>% select(LWKLYWGE,EDUC, YR20, YR21, YR22, YR23, YR24, YR25, YR26, YR27, YR28, YR29) %>% pull

n <- nrow(Selected_Data)

Y <- Selected_Data %>% select(LWKLYWGE) %>%  as.matrix
A <- Selected_Data %>% select(YR20, YR21, YR22, YR23, YR24, YR25, YR26, YR27, YR28, YR29) %>% as.matrix
Z <- Selected_Data %>% select(EDUC,YR20, YR21, YR22, YR23, YR24, YR25, YR26, YR27, YR28, YR29) %>%  as.matrix

str(Y)
K_class(0,A,Z,Y,n)
K_class <- function(k,A,Z,Y,n){
  P_A <- A%*%solve(t(A)%*%A)%*%t(A)
  Wk <- t(Z) %*% ((1-k)*diag(n)+k* P_A )  %*%Z
  a <- solve(Wk) %*% t(Z) %*% ( (1-k)*diag(n) + k*P_A ) %*% Y 
  return(a)
}

Selected_DataSelected_Data

lm(LWKLYWGE ~ EDUC + YR20 + YR21+ YR22+ YR23+ YR24+ YR25+ YR26+ YR27+ YR28+ YR29,data=Selected_Data)
ivreg(LWKLYWGE ~ EDUC + YR20 + YR21+ YR22+ YR23+ YR24+ YR25+ YR26+ YR27+ YR28+ YR29 | YR20 + YR21+ YR22+ YR23+ YR24+ YR25+ YR26+ YR27+ YR28+ YR29,
      data= Selected_Data)
