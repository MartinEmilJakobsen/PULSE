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

COLONIAL_Data <- read.delim(file="Data/COLONIAL_T4.csv",sep=",",header=TRUE)  %>% mutate(Intercept = 1)
###################################
###################################
###################################
######## TABLE IV: COLONIAL #######
###################################
###################################
###################################



Selected_Data <- COLONIAL_Data  %>% filter(baseco == 1)

############################################################################
#COLUMN 1 : logpgp95 ~ avexpr (INSTRUMENTS = logem4)#
############################################################################
n <- nrow(Selected_Data)
#Target
Y <- Selected_Data %>% select(logpgp95) %>%  as.matrix
#Included Exogenous:
A_1 <- Selected_Data %>% select(Intercept) %>%  as.matrix
#All Exogenous
A <- Selected_Data %>% select(logem4,Intercept) %>% as.matrix
#Included Endogenous
X <- Selected_Data %>% select(avexpr) %>%  as.matrix
#Included (all)
Z <- cbind(X,A_1)


#OLS
K_class(0,A,Z,Y,n)

#IV
K_class(1,A,Z,Y,n)

#PULSE
PULSE(A,A_1,X,Y,p=0.05,N=10000,n)


K_class(FULLER_k(4,A,A_1,X,Y,n,2),A,Z,Y,n)


#LM and IVREG check
lm(logpgp95~avexpr,data=Selected_Data)
ivfit<-ivreg(logpgp95~avexpr|logem4,data=Selected_Data,x=TRUE)
summary(ivfit)

anderson.rubin.ci(ivfit)

############################################################################
#COLUMN 2 : logpgp95 ~ avexpr + lat_abst (INSTRUMENTS = logem4 + lat_abst)#
############################################################################
n <- nrow(Selected_Data)
#Target
Y <- Selected_Data %>% select(logpgp95) %>%  as.matrix
#Included Exogenous:
A_1 <- Selected_Data %>% select(lat_abst,Intercept) %>%  as.matrix
#All Exogenous
A <- Selected_Data %>% select(lat_abst,logem4,Intercept) %>% as.matrix
#Included Endogenous
X <- Selected_Data %>% select(avexpr) %>%  as.matrix
#Included (all)
Z <- cbind(X,A_1)


#OLS
K_class(0,A,Z,Y,n)

#IV
K_class(1,A,Z,Y,n)

#PULSE
PULSE(A,A_1,X,Y,p=0.05,N=10000,n)

#LM and IVREG check
lm(logpgp95~avexpr+lat_abst,data=Selected_Data)
ivfit<-ivreg(logpgp95~avexpr+lat_abst|logem4+lat_abst,data=Selected_Data,x=TRUE)
summary(ivfit)

anderson.rubin.ci(ivfit)


############################################################################
#COLUMN 3 : logpgp95 ~ avexpr (INSTRUMENTS = logem4) w/o Neo-Europes
############################################################################
Selected_Data_woNeoEuropes <- Selected_Data %>%  filter(rich4!= 1)
n <- nrow(Selected_Data_woNeoEuropes)
#Target
Y <- Selected_Data_woNeoEuropes %>% select(logpgp95) %>%  as.matrix
#Included Exogenous:
A_1 <- Selected_Data_woNeoEuropes %>% select(Intercept) %>%  as.matrix
#All Exogenous
A <- Selected_Data_woNeoEuropes %>% select(logem4,Intercept) %>% as.matrix
#Included Endogenous
X <- Selected_Data_woNeoEuropes %>% select(avexpr) %>%  as.matrix
#Included (all)
Z <- cbind(X,A_1)


#OLS
K_class(0,A,Z,Y,n)

#IV
K_class(1,A,Z,Y,n)

#PULSE
PULSE(A,A_1,X,Y,p=0.05,N=10000,n)

#LM and IVREG check
lm(logpgp95~avexpr,data=Selected_Data_woNeoEuropes)
ivfit<-ivreg(logpgp95~avexpr|logem4,data=Selected_Data_woNeoEuropes,x=TRUE)
summary(ivfit)

anderson.rubin.ci(ivfit)

############################################################################
#COLUMN 4 :  logpgp95 ~ avexpr + lat_abst (INSTRUMENTS = logem4 + lat_abst) w/o Neo-Europes
############################################################################
Selected_Data_woNeoEuropes <- Selected_Data %>%  filter(rich4!= 1)
n <- nrow(Selected_Data_woNeoEuropes)
#Target
Y <- Selected_Data_woNeoEuropes %>% select(logpgp95) %>%  as.matrix
#Included Exogenous:
A_1 <- Selected_Data_woNeoEuropes %>% select(lat_abst,Intercept) %>%  as.matrix
#All Exogenous
A <- Selected_Data_woNeoEuropes %>% select(lat_abst,logem4,Intercept) %>% as.matrix
#Included Endogenous
X <- Selected_Data_woNeoEuropes %>% select(avexpr) %>%  as.matrix
#Included (all)
Z <- cbind(X,A_1)

which(X<=16)
X[which(X<=16)]

#OLS
K_class(0,A,Z,Y,n)

#IV
K_class(1,A,Z,Y,n)

#PULSE
PULSE(A,A_1,X,Y,p=0.05,N=10000,n)

#LM and IVREG check
lm(logpgp95~avexpr+lat_abst,data=Selected_Data_woNeoEuropes)
ivfit<-ivreg(logpgp95~avexpr+lat_abst|logem4+lat_abst,data=Selected_Data_woNeoEuropes,x=TRUE)
summary(ivfit)

anderson.rubin.ci(ivfit)

############################################################################
#COLUMN 5 : logpgp95 ~ avexpr (INSTRUMENTS = logem4) w/o Africa
############################################################################
Selected_Data_woAfrica <- Selected_Data %>%  filter(africa!= 1)
n <- nrow(Selected_Data_woAfrica)
#Target
Y <- Selected_Data_woAfrica %>% select(logpgp95) %>%  as.matrix
#Included Exogenous:
A_1 <- Selected_Data_woAfrica %>% select(Intercept) %>%  as.matrix
#All Exogenous
A <- Selected_Data_woAfrica %>% select(logem4,Intercept) %>% as.matrix
#Included Endogenous
X <- Selected_Data_woAfrica %>% select(avexpr) %>%  as.matrix
#Included (all)
Z <- cbind(X,A_1)


#OLS
K_class(0,A,Z,Y,n)

#IV
K_class(1,A,Z,Y,n)

#PULSE
PULSE(A,A_1,X,Y,p=0.05,N=10000,n)

#LM and IVREG check
lm(logpgp95~avexpr,data=Selected_Data_woAfrica)
ivfit<-ivreg(logpgp95~avexpr|logem4,data=Selected_Data_woAfrica,x=TRUE)
summary(ivfit)

anderson.rubin.ci(ivfit)


############################################################################
#COLUMN 6 :  logpgp95 ~ avexpr + lat_abst (INSTRUMENTS = logem4 + lat_abst) w/o Africa
############################################################################
Selected_Data_woAfrica <- Selected_Data %>%  filter(africa!= 1)
n <- nrow(Selected_Data_woAfrica)
#Target
Y <- Selected_Data_woAfrica %>% select(logpgp95) %>%  as.matrix
#Included Exogenous:
A_1 <- Selected_Data_woAfrica %>% select(lat_abst,Intercept) %>%  as.matrix
#All Exogenous
A <- Selected_Data_woAfrica %>% select(lat_abst,logem4,Intercept) %>% as.matrix
#Included Endogenous
X <- Selected_Data_woAfrica %>% select(avexpr) %>%  as.matrix
#Included (all)
Z <- cbind(X,A_1)


#OLS
K_class(0,A,Z,Y,n)

#IV
K_class(1,A,Z,Y,n)

#PULSE
PULSE(A,A_1,X,Y,p=0.05,N=10000,n)

#LM and IVREG check
lm(logpgp95~avexpr+lat_abst,data=Selected_Data_woAfrica)
ivfit<-ivreg(logpgp95~avexpr+lat_abst|logem4+lat_abst,data=Selected_Data_woAfrica,x=TRUE)
summary(ivfit)

anderson.rubin.ci(ivfit)


############################################################################
#COLUMN 7 :  logpgp95 ~ avexpr + africa + asia + other_cont (INSTRUMENTS = africa + asia + other_cont+ logem4 ) w/continentIndicator
############################################################################
Selected_Data_wContinent <- Selected_Data %>%  mutate(other_cont = ifelse(shortnam %in% c("AUS","MLT","NZL"),1,0))


n <- nrow(Selected_Data_wContinent)
#Target
Y <- Selected_Data_wContinent %>% select(logpgp95) %>%  as.matrix
#Included Exogenous:
A_1 <- Selected_Data_wContinent %>% select(Intercept,africa,asia,other_cont) %>%  as.matrix
#All Exogenous
A <- Selected_Data_wContinent %>% select(Intercept,africa,asia,other_cont,logem4) %>% as.matrix
#Included Endogenous
X <- Selected_Data_wContinent %>% select(avexpr) %>%  as.matrix
#Included (all)
Z <- cbind(X,A_1)


#OLS
K_class(0,A,Z,Y,n)

#IV
K_class(1,A,Z,Y,n)

#PULSE
PULSE(A,A_1,X,Y,p=0.05,N=10000,n)

#LM and IVREG check
lm(logpgp95~avexpr + africa + asia + other_cont,data=Selected_Data_wContinent)
ivfit<-ivreg(logpgp95~avexpr + africa + asia + other_cont|logem4+africa + asia + other_cont,data=Selected_Data_wContinent,x=TRUE)
summary(ivfit)

anderson.rubin.ci(ivfit)

############################################################################
#COLUMN 8 :  logpgp95 ~ avexpr + lat_abst + africa + asia + other_cont (INSTRUMENTS = lat_abst + africa + asia + other_cont ) w/continentIndicator
############################################################################
Selected_Data_wContinent <- Selected_Data %>%  mutate(other_cont = ifelse(shortnam %in% c("AUS","MLT","NZL"),1,0))


n <- nrow(Selected_Data_wContinent)
#Target
Y <- Selected_Data_wContinent %>% select(logpgp95) %>%  as.matrix
#Included Exogenous:
A_1 <- Selected_Data_wContinent %>% select(Intercept,africa,asia,other_cont,lat_abst) %>%  as.matrix
#All Exogenous
A <- Selected_Data_wContinent %>% select(Intercept,africa,asia,other_cont,logem4,lat_abst) %>% as.matrix
#Included Endogenous
X <- Selected_Data_wContinent %>% select(avexpr) %>%  as.matrix
#Included (all)
Z <- cbind(X,A_1)


#OLS
K_class(0,A,Z,Y,n)

#IV
K_class(1,A,Z,Y,n)

#PULSE
PULSE(A,A_1,X,Y,p=0.05,N=10000,n)

#LM and IVREG check
lm(logpgp95~avexpr + africa + asia + other_cont+lat_abst,data=Selected_Data_wContinent)
ivfit<-ivreg(logpgp95~avexpr + africa + asia + other_cont+lat_abst| africa + asia + other_cont+logem4+lat_abst,data=Selected_Data_wContinent,x=TRUE)
summary(ivfit)

anderson.rubin.ci(ivfit)

############################################################################
#COLUMN 9 :  loghjypl ~ avexpr (INSTRUMENTS = logem4 )
############################################################################

#removing NA log output per worker samples
Selected_Data_NAremoved <- Selected_Data %>%  filter(!is.na(loghjypl))

n <- nrow(Selected_Data_NAremoved)
#Target
Y <- Selected_Data_NAremoved %>% select(loghjypl) %>%  as.matrix
#Included Exogenous:
A_1 <- Selected_Data_NAremoved %>% select(Intercept) %>% as.matrix
#All Exogenous
A <- Selected_Data_NAremoved %>% select(Intercept,logem4) %>% as.matrix
#Included Endogenous
X <- Selected_Data_NAremoved %>% select(avexpr) %>%  as.matrix
#Included (all)
Z <- cbind(X,A_1)


#OLS
K_class(0,A,Z,Y,n)

#IV
K_class(1,A,Z,Y,n)

#PULSE
PULSE(A,A_1,X,Y,p=0.05,N=10000,n)

#LM and IVREG check
lm(loghjypl~avexpr,data=Selected_Data_NAremoved)
ivfit<-ivreg(loghjypl~avexpr | logem4,data=Selected_Data_NAremoved,x=TRUE)
summary(ivfit)

anderson.rubin.ci(ivfit)
