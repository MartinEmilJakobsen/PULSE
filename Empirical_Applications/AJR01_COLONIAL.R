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

`%notin%` <- Negate(`%in%`)

source("../Estimators_Slow.R")

COLONIAL_Data <- read.delim(file="Data/AJR01_COLONIAL.csv",sep=",",header=TRUE)  %>% mutate(Intercept = 1)
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
m1.ols <- K_class(0,A,Z,Y,n)

#IV
m1.tsls <- K_class(1,A,Z,Y,n)

#PULSE
m1.pulse <- PULSE(A,A_1,X,Y,p=0.05,N=10000,n)

#fuller4
dA <- ncol(A)
m1.fuller4 <-  K_class(FULLER_k(4,A,A_1,X,Y,n,dA),A,Z,Y,n)

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
m2.ols <- K_class(0,A,Z,Y,n)

#IV
m2.tsls <- K_class(1,A,Z,Y,n)

#PULSE
m2.pulse <- PULSE(A,A_1,X,Y,p=0.05,N=10000,n)

#fuller4
dA <- ncol(A)
m2.fuller4 <-  K_class(FULLER_k(4,A,A_1,X,Y,n,dA),A,Z,Y,n)


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
m3.ols <- K_class(0,A,Z,Y,n)

#IV
m3.tsls <- K_class(1,A,Z,Y,n)

#PULSE
m3.pulse <- PULSE(A,A_1,X,Y,p=0.05,N=10000,n)

#fuller4
dA <- ncol(A)
m3.fuller4 <-  K_class(FULLER_k(4,A,A_1,X,Y,n,dA),A,Z,Y,n)


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

#OLS
m4.ols <- K_class(0,A,Z,Y,n)

#IV
m4.tsls <- K_class(1,A,Z,Y,n)

#PULSE
m4.pulse <- PULSE(A,A_1,X,Y,p=0.05,N=10000,n)

#fuller4
dA <- ncol(A)
m4.fuller4 <-  K_class(FULLER_k(4,A,A_1,X,Y,n,dA),A,Z,Y,n)


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
m5.ols <- K_class(0,A,Z,Y,n)

#IV
m5.tsls <- K_class(1,A,Z,Y,n)

#PULSE
m5.pulse <- PULSE(A,A_1,X,Y,p=0.05,N=10000,n)

#fuller4
dA <- ncol(A)
m5.fuller4 <-  K_class(FULLER_k(4,A,A_1,X,Y,n,dA),A,Z,Y,n)


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
m6.ols <- K_class(0,A,Z,Y,n)

#IV
m6.tsls <- K_class(1,A,Z,Y,n)

#PULSE
m6.pulse <- PULSE(A,A_1,X,Y,p=0.05,N=10000,n)

#fuller4
dA <- ncol(A)
m6.fuller4 <-  K_class(FULLER_k(4,A,A_1,X,Y,n,dA),A,Z,Y,n)


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
m7.ols <-K_class(0,A,Z,Y,n)

#IV
m7.tsls <- K_class(1,A,Z,Y,n)

#PULSE
m7.pulse <- PULSE(A,A_1,X,Y,p=0.05,N=10000,n)

#fuller4
dA <- ncol(A)
m7.fuller4 <-  K_class(FULLER_k(4,A,A_1,X,Y,n,dA),A,Z,Y,n)


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
m8.ols <- K_class(0,A,Z,Y,n)

#IV
m8.tsls <- K_class(1,A,Z,Y,n)

#PULSE
m8.pulse <- PULSE(A,A_1,X,Y,p=0.05,N=10000,n)

#fuller4
dA <- ncol(A)
m8.fuller4 <-  K_class(FULLER_k(4,A,A_1,X,Y,n,dA),A,Z,Y,n)


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


#############
Table <- data.frame(
  OLS = c(m1.ols["avexpr",],
          m2.ols["avexpr",],
          m3.ols["avexpr",],
          m4.ols["avexpr",],
          m5.ols["avexpr",],
          m6.ols["avexpr",],
          m7.ols["avexpr",],
          m8.ols["avexpr",]),
  TSLS = c(m1.tsls["avexpr",],
           m2.tsls["avexpr",],
           m3.tsls["avexpr",],
           m4.tsls["avexpr",],
           m5.tsls["avexpr",],
           m6.tsls["avexpr",],
           m7.tsls["avexpr",],
           m8.tsls["avexpr",]),
  FULLER4 = c(m1.fuller4["avexpr",],
           m2.fuller4["avexpr",],
           m3.fuller4["avexpr",],
           m4.fuller4["avexpr",],
           m5.fuller4["avexpr",],
           m6.fuller4["avexpr",],
           m7.fuller4["avexpr",],
           m8.fuller4["avexpr",]),
  PULSE = c(m1.pulse["avexpr","logpgp95"],
            m2.pulse["avexpr","logpgp95"],
            m3.pulse["avexpr","logpgp95"],
            m4.pulse["avexpr","logpgp95"],
            m5.pulse["avexpr","logpgp95"],
            m6.pulse["avexpr","logpgp95"],
            m7.pulse["avexpr","logpgp95"],
            m8.pulse["avexpr","logpgp95"]),
  message = c(m1.pulse["avexpr","m"],
              m2.pulse["avexpr","m"],
              m3.pulse["avexpr","m"],
              m4.pulse["avexpr","m"],
              m5.pulse["avexpr","m"],
              m6.pulse["avexpr","m"],
              m7.pulse["avexpr","m"],
              m8.pulse["avexpr","m"]),
  test = c(m1.pulse["avexpr","t"],
           m2.pulse["avexpr","t"],
           m3.pulse["avexpr","t"],
           m4.pulse["avexpr","t"],
           m5.pulse["avexpr","t"],
           m6.pulse["avexpr","t"],
           m7.pulse["avexpr","t"],
           m8.pulse["avexpr","t"]),
  threshold = c(m1.pulse["avexpr","q"],
                m2.pulse["avexpr","q"],
                m3.pulse["avexpr","q"],
                m4.pulse["avexpr","q"],
                m5.pulse["avexpr","q"],
                m6.pulse["avexpr","q"],
                m7.pulse["avexpr","q"],
                m8.pulse["avexpr","q"]))

kbl(Table,format="latex",digits=4)

