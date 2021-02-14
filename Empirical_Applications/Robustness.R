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


`%notin%` <- Negate(`%in%`)


getwd()
source("../Estimators_Slow.R")


################################################
################################################
################################################
######## COLONIAL ROBUSTNESS EXPERIMENTS #######
################################################
################################################
################################################

COLONIAL_Data <- read.delim(file="Data/COLONIAL_T4.csv",sep=",",header=TRUE)



Selected_Data <- COLONIAL_Data  %>% filter(baseco == 1) %>% mutate_at(c("lat_abst","avexpr","logpgp95","logem4"),.funs=function(x){x-mean(x)}) %>% mutate(Intercept=1)

 #############################################################################################
# MODEL M1 WITHOUT INTERCEPT: logpgp95 ~ avexpr  (INSTRUMENTS = logem4)#
#############################################################################################
holdoutextremes <- function(Data,noExtremes){

    HoldOutDataset <- Data %>% 
    as_tibble() %>%  
    arrange(logem4) %>% 
    head(-noExtremes/2) %>% 
    tail(-noExtremes/2)
  HeldOutSamples <- Data %>% 
    filter(shortnam %notin% HoldOutDataset$shortnam)
  
  n <- nrow(HoldOutDataset)
  #Target
  Y <- HoldOutDataset %>% select(logpgp95) %>%  as.matrix
  #Included Exogenous:
  A_1 <- HoldOutDataset %>% select() %>%  as.matrix
  #All Exogenous
  A <- HoldOutDataset %>% select(logem4) %>% as.matrix
  #Included Endogenous
  X <- HoldOutDataset %>% select(avexpr) %>%  as.matrix
  #Included (all)
  Z <- cbind(X,A_1)
  
  #OLS
  ols <- K_class(0,A,Z,Y,n)
  #IV
  iv <- K_class(1,A,Z,Y,n)
  #PULSE
  pulse <- PULSE(A,A_1,X,Y,p=0.05,N=10000,n)
  
  
  MSPE <- HeldOutSamples %>%  
    select(logpgp95,avexpr,Intercept) %>% 
    mutate(resOLS = logpgp95- (avexpr*ols["avexpr",]),
           resIV = logpgp95- (avexpr*iv["avexpr",]),
           resPULSE = logpgp95- (avexpr*pulse["avexpr","logpgp95"])) %>% 
    summarize(OLS = sum(resOLS^2)/n(),
              IV = sum(resIV^2)/n(),
              PULSE = sum(resPULSE^2)/n(),
              ORDER = case_when(OLS <= IV & IV <= PULSE ~ "OLS<IV<PULSE",
                                OLS < PULSE & PULSE <= IV ~ "OLS<PULSE<IV",
                                OLS == PULSE & PULSE <= IV ~ "PULSE=OLS<IV",
                                IV <= PULSE & PULSE < OLS ~ "IV<PULSE<OLS",
                                IV <= PULSE & PULSE == OLS ~ "IV<PULSE=OLS",
                                IV <= OLS & OLS <= PULSE ~ "IV<OLS<PULSE",
                                IV <= OLS & OLS == PULSE ~ "IV<OLS=PULSE",
                                PULSE < OLS & OLS <= IV ~ "PULSE<OLS<IV",
                                PULSE == OLS & OLS <= IV ~ "PULSE=OLS<IV",
                                PULSE <= IV & IV <= OLS ~ "PULSE<IV<OLS"
              ),
              OLSmIV = OLS-IV,
              OLSmPULSE = OLS-PULSE,
              PULSEmIV = PULSE -IV,
              coef.OLS = ols["avexpr",],
              coef.IV = iv["avexpr",],
              coef.PULSE = pulse["avexpr","logpgp95"]) %>% 
    mutate(t=pulse["avexpr","t"],q=pulse["avexpr","q"],l=pulse["avexpr","l"],k=pulse["avexpr","k"])
  return(MSPE)
}

summarydat <-data.frame(n=seq(2,40,2)) %>% 
  rowwise() %>% 
  mutate(order=list(holdoutextremes(Selected_Data,n)))%>%
  unnest(cols=c(order))

summarydat %>%  select(n,OLS,IV,PULSE,coef.OLS,coef.IV,coef.PULSE,ORDER) %>% 
  kbl(.,format="latex",digits=4)


#################
#Hold out random#
#################

set.seed(1)

holdoutrandom <- function(Data,percentage){
  HoldOutDataset <- Data %>% 
    as_tibble() %>%  
    sample_frac(percentage)
  HeldOutSamples <- Data %>% 
    filter(shortnam %notin% HoldOutDataset$shortnam)
  
  n <- nrow(HoldOutDataset)
  #Target
  Y <- HoldOutDataset %>% select(logpgp95) %>%  as.matrix
  #Included Exogenous:
  A_1 <- HoldOutDataset %>% select() %>%  as.matrix
  #All Exogenous
  A <- HoldOutDataset %>% select(logem4) %>% as.matrix
  #Included Endogenous
  X <- HoldOutDataset %>% select(avexpr) %>%  as.matrix
  #Included (all)
  Z <- cbind(X,A_1)
  
  #OLS
  ols <- K_class(0,A,Z,Y,n)
  #IV
  iv <- K_class(1,A,Z,Y,n)
  #PULSE
  pulse <- PULSE(A,A_1,X,Y,p=0.05,N=10000,n)
  
  
  MSPE <- HeldOutSamples %>%  
    select(logpgp95,avexpr,Intercept) %>% 
    mutate(resOLS = logpgp95- (avexpr*ols["avexpr",]),
           resIV = logpgp95- (avexpr*iv["avexpr",]),
           resPULSE = logpgp95- (avexpr*pulse["avexpr","logpgp95"])) %>% 
    summarize(OLS = sum(resOLS^2)/n(),
              IV = sum(resIV^2)/n(),
              PULSE = sum(resPULSE^2)/n(),
              ORDER = case_when(pulse["avexpr","l"]== 0 & PULSE <= IV ~ "PULSE=OLS<IV",
                                pulse["avexpr","l"]== 0 & OLS <= IV ~ "PULSE=OLS<IV",
                                IV <= OLS & pulse["avexpr","l"] == 0 ~ "IV<OLS=PULSE",
                                IV <= PULSE & pulse["avexpr","l"] == 0 ~ "IV<PULSE=OLS",
                                OLS <= IV & IV <= PULSE ~ "OLS<IV<PULSE",
                                OLS < PULSE & PULSE <= IV ~ "OLS<PULSE<IV",
                                IV <= PULSE & PULSE < OLS ~ "IV<PULSE<OLS",
                                IV <= OLS & OLS <= PULSE ~ "IV<OLS<PULSE",
                                PULSE < OLS & OLS <= IV ~ "PULSE<OLS<IV",
                                PULSE <= IV & IV <= OLS ~ "PULSE<IV<OLS"
              ),
              OLSmIV = OLS-IV,
              OLSmPULSE = OLS-PULSE,
              PULSEmIV = PULSE -IV) %>% 
    mutate(t=pulse["avexpr","t"],q=pulse["avexpr","q"],l=pulse["avexpr","l"],k=pulse["avexpr","k"])
  return(MSPE)
}

norep <- 1000

summarydat <-data.frame(n=seq(1,norep,1)) %>% 
  rowwise() %>% 
  mutate(order=list(holdoutrandom(Selected_Data,0.9)))


summarydat %>%
  unnest(cols=c(order)) %>% 
  mutate(OLSeqPULSE = I(t!=q)) %>% 
  group_by(ORDER,OLSeqPULSE) %>% 
  summarise(count= n(),
            WMSPE.OLS = max(OLS),
            WMSPE.iv = max(IV),
            WMSPE.pulse = max(PULSE),
            OLS = mean(OLS),
            PULSE = mean(PULSE),
            IV = mean(IV)
            ) %>% 
  mutate(FirstDiff = case_when(ORDER == "OLS<IV<PULSE" ~ IV-OLS,
                               ORDER == "OLS<PULSE<IV" ~  PULSE- OLS,
                               ORDER == "OLS=PULSE<IV" ~ PULSE - OLS,
                               ORDER == "IV<PULSE<OLS" ~ PULSE - IV,
                               ORDER == "IV<PULSE=OLS" ~ PULSE -IV,
                               ORDER == "IV<OLS<PULSE" ~ OLS - IV,
                               ORDER == "IV<OLS=PULSE" ~ OLS - IV,
                               ORDER == "PULSE<OLS<IV" ~ OLS - PULSE,
                               ORDER == "PULSE=OLS<IV" ~ OLS -PULSE,
                               ORDER == "PULSE<IV<OLS" ~ IV-PULSE),
         SecondDiff = case_when(ORDER == "OLS<IV<PULSE" ~ PULSE -IV,
                                ORDER == "OLS<PULSE<IV" ~  IV -PULSE,
                                ORDER == "OLS=PULSE<IV" ~ IV-PULSE,
                                ORDER == "IV<PULSE<OLS" ~ OLS-PULSE,
                                ORDER == "IV<PULSE=OLS" ~ OLS-PULSE,
                                ORDER == "IV<OLS<PULSE" ~ PULSE-OLS,
                                ORDER == "IV<OLS=PULSE" ~ PULSE-OLS,
                                ORDER == "PULSE<OLS<IV" ~ IV-OLS,
                                ORDER == "PULSE=OLS<IV" ~ IV-OLS,
                                ORDER == "PULSE<IV<OLS" ~ OLS-IV)) %>% 
  
  mutate(Percentage = 100*count/norep) %>% 
  select(ORDER,Percentage,OLS,IV,PULSE,WMSPE.OLS,WMSPE.iv,WMSPE.pulse) %>% 
  kbl(.,format="latex",digits=3)

################################################
################################################
################################################
########   ANGRIST AND KRUGER 1991 - QOB #######
################################################
################################################
################################################


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

Selected_Data <- QOB_Data  %>% filter(COHORT == "30.39") %>% 
  select(-c(V3,V7,V14,V15,V17,V22,V23))  %>% mutate(ID=1:n())

Selected_Data_Centered <- Selected_Data %>% mutate_at(c("LWKLYWGE","AGEQ","AGEQSQ","EDUC"),.funs=function(x){x-mean(x)}) %>%  mutate(ID=1:n())


##### CUSTOM MODEL M4 #### 
# MAYBE DELETE THIS.

holdoutextremes <- function(Data,noExtremes){
  HoldOutDataset <- Data %>% 
    as_tibble() %>%  
    arrange(AGEQ) %>% 
    head(-noExtremes/2) %>% 
    tail(-noExtremes/2)
  HeldOutSamples <- Data %>% 
    filter(ID %notin% HoldOutDataset$ID)

  
    n <- nrow(HoldOutDataset)
  #Target
  Y <- HoldOutDataset %>% select(LWKLYWGE) %>%  as.matrix
  #Included Exogenous:
  A_1 <- HoldOutDataset %>% select(AGEQ,AGEQSQ) %>% as.matrix
  #All Exogenous
  A <- HoldOutDataset %>% select(AGEQ,AGEQSQ,QTR1,QTR2,QTR3) %>% as.matrix
  #Included Endogenous
  X <- HoldOutDataset %>% select(EDUC) %>%  as.matrix
  #Included (all)
  Z <- cbind(X,A_1)
  
  
  #OLS
  ols <- K_class(0,A,Z,Y,n)
  #IV
  iv <- K_class(1,A,Z,Y,n)
  #PULSE
  pulse <- PULSE(A,A_1,X,Y,p=0.05,N=10000,n)
  
  
  MSPE <- HeldOutSamples %>%  
    select(LWKLYWGE,AGEQ,AGEQSQ,QTR1,QTR2,QTR3,EDUC,Intercept) %>% 
    mutate(resOLS = LWKLYWGE - (EDUC*ols["EDUC",]+AGEQ*ols["AGEQ",]+AGEQSQ*ols["AGEQSQ",]),
           resIV = LWKLYWGE - (EDUC*iv["EDUC",]+AGEQ*iv["AGEQ",]+AGEQSQ*iv["AGEQSQ",]),
           resPULSE = LWKLYWGE - (EDUC*pulse["EDUC","LWKLYWGE"]+AGEQ*pulse["AGEQ","LWKLYWGE"]+AGEQSQ*pulse["AGEQSQ","LWKLYWGE"])) %>% 
    summarize(OLS = sum(resOLS^2)/n(),
              IV = sum(resIV^2)/n(),
              PULSE = sum(resPULSE^2)/n(),
              ORDER = case_when(OLS <= IV & IV <= PULSE ~ "OLS<IV<PULSE",
                                OLS < PULSE & PULSE <= IV ~ "OLS<PULSE<IV",
                                OLS == PULSE & PULSE <= IV ~ "PULSE=OLS<IV",
                                IV <= PULSE & PULSE < OLS ~ "IV<PULSE<OLS",
                                IV <= PULSE & PULSE == OLS ~ "IV<PULSE=OLS",
                                IV <= OLS & OLS <= PULSE ~ "IV<OLS<PULSE",
                                IV <= OLS & OLS == PULSE ~ "IV<OLS=PULSE",
                                PULSE < OLS & OLS <= IV ~ "PULSE<OLS<IV",
                                PULSE == OLS & OLS <= IV ~ "PULSE=OLS<IV",
                                PULSE <= IV & IV <= OLS ~ "PULSE<IV<OLS"
              ),
              OLSmIV = OLS-IV,
              OLSmPULSE = OLS-PULSE,
              PULSEmIV = PULSE -IV,
              coef.OLS = ols["EDUC",],
              coef.IV = iv["EDUC",],
              coef.PULSE = pulse["EDUC","LWKLYWGE"])%>% 
    mutate(t=pulse["EDUC","t"],q=pulse["EDUC","q"])
  return(MSPE)
}

summarydat <-data.frame(n=c(50000,100000)) %>% 
  rowwise() %>% 
  mutate(order=list(holdoutextremes(Selected_Data_Centered,n)))%>%
  unnest(cols=c(order))

summarydat %>%  select(n,OLS,IV,PULSE,ORDER,coef.OLS,coef.IV,coef.PULSE) %>% 
  kbl(.,format="latex",digits=4)


######## Hold out random #########



set.seed(1)

holdoutrandom <- function(Data,percentage){
  HoldOutDataset <- Data %>% 
    as_tibble() %>%  
    sample_frac(percentage)
  HeldOutSamples <- Data %>% 
    filter(ID %notin% HoldOutDataset$ID)
  
 
  n <- nrow(HoldOutDataset)
  #Target
  Y <- HoldOutDataset %>% select(LWKLYWGE) %>%  as.matrix
  #Included Exogenous:
  A_1 <- HoldOutDataset %>% select(starts_with("YR")) %>% as.matrix
  #All Exogenous
  A <- HoldOutDataset %>% select(starts_with("YR"),starts_with("QTRxYR")) %>% as.matrix
  #Included Endogenous
  X <- HoldOutDataset %>% select(EDUC) %>%  as.matrix
  #Included (all)
  Z <- cbind(X,A_1)
  
  #OLS
  ols <- K_class(0,A,Z,Y,n)
  #IV
  iv <- K_class(1,A,Z,Y,n)
  #PULSE
  pulse <- PULSE(A,A_1,X,Y,p=0.05,N=10000,n)
  
  
  MSPE <- HeldOutSamples %>%  
    select(LWKLYWGE,EDUC,YR0 ,  YR1 ,  YR2 ,  YR3 ,  YR4 ,  YR5  , YR6  , YR7 ,  YR8,   YR9) %>% 
    mutate(resOLS = LWKLYWGE- (EDUC*ols["EDUC",]+ YR0*ols["YR0",]+YR1*ols["YR1",]+YR2*ols["YR2",]+YR3*ols["YR3",]+YR4*ols["YR4",]+YR5*ols["YR5",]+YR6*ols["YR6",]+YR7*ols["YR7",]+YR8*ols["YR8",]+YR9*ols["YR9",]),
           resIV = LWKLYWGE- (EDUC*iv["EDUC",]+ YR0*iv["YR0",]+YR1*iv["YR1",]+YR2*iv["YR2",]+YR3*iv["YR3",]+YR4*iv["YR4",]+YR5*iv["YR5",]+YR6*iv["YR6",]+YR7*iv["YR7",]+YR8*iv["YR8",]+YR9*iv["YR9",]),
           resPULSE = LWKLYWGE- (EDUC*pulse["EDUC","LWKLYWGE"]+ YR0*pulse["YR0","LWKLYWGE"]+YR1*pulse["YR1","LWKLYWGE"]+YR2*pulse["YR2","LWKLYWGE"]+YR3*pulse["YR3","LWKLYWGE"]+YR4*pulse["YR4","LWKLYWGE"]+YR5*pulse["YR5","LWKLYWGE"]+YR6*pulse["YR6","LWKLYWGE"]+YR7*pulse["YR7","LWKLYWGE"]+YR8*pulse["YR8","LWKLYWGE"]+YR9*pulse["YR9","LWKLYWGE"])) %>% 
    summarize(OLS = sum(resOLS^2)/n(),
              IV = sum(resIV^2)/n(),
              PULSE = sum(resPULSE^2)/n(),
              ORDER = case_when(pulse["EDUC","l"]== 0 & PULSE <= IV ~ "PULSE=OLS<IV",
                                pulse["EDUC","l"]== 0 & OLS <= IV ~ "PULSE=OLS<IV",
                                IV <= OLS & pulse["EDUC","l"] == 0 ~ "IV<OLS=PULSE",
                                IV <= PULSE & pulse["EDUC","l"] == 0 ~ "IV<PULSE=OLS",
                                OLS <= IV & IV <= PULSE ~ "OLS<IV<PULSE",
                                OLS < PULSE & PULSE <= IV ~ "OLS<PULSE<IV",
                                IV <= PULSE & PULSE < OLS ~ "IV<PULSE<OLS",
                                IV <= OLS & OLS <= PULSE ~ "IV<OLS<PULSE",
                                PULSE < OLS & OLS <= IV ~ "PULSE<OLS<IV",
                                PULSE <= IV & IV <= OLS ~ "PULSE<IV<OLS"
              ),
              OLSmIV = OLS-IV,
              OLSmPULSE = OLS-PULSE,
              PULSEmIV = PULSE -IV) %>% 
    mutate(t=pulse["EDUC","t"],q=pulse["EDUC","q"],l=pulse["EDUC","l"],k=pulse["EDUC","k"])
  return(MSPE)
}

norep <- 1000

Selected_Data_subset <- Selected_Data %>% sample_frac(0.01)

summarydat <-data.frame(n=seq(1,norep,1)) %>% 
  rowwise() %>% 
  mutate(order=list(holdoutrandom(Selected_Data_subset  ,0.9)))


AverageStats <- summarydat %>%
  unnest(cols=c(order)) %>% 
  group_by(ORDER) %>% 
  summarise(count= n(),
            OLS = mean(OLS),
            PULSE = mean(PULSE),
            IV = mean(IV)  ) %>% 
  mutate(Percentage = 100*count/norep)
 

WorstCaseAverageStats <- bind_rows(
summarydat %>%
  unnest(cols=c(order)) %>%
  filter(IV <= PULSE & IV <= OLS) %>% 
  arrange(desc(IV)) %>% 
  filter(row_number() <= n()*0.1) %>% 
  group_by(ORDER) %>% 
  summarise(OLS = mean(OLS),
            PULSE = mean(PULSE),
            IV = mean(IV))
,
summarydat %>%
  unnest(cols=c(order)) %>%
  filter(OLS <= IV ) %>% 
  arrange(desc(IV)) %>% 
  filter(row_number() <= n()*0.1) %>% 
  group_by(ORDER) %>% 
  summarise(OLS = mean(OLS),
            PULSE = mean(PULSE),
            IV = mean(IV)))

left_join(AverageStats,WorstCaseAverageStats,by="ORDER") %>% 
select(ORDER,Percentage,OLS.x,IV.x,PULSE.x,OLS.y,IV.y,PULSE.y) %>% 
  kbl(.,format="latex",digits=3)

################################################
################################################
################################################
######   CARD 1995 - Proximity to college ######
################################################
################################################
################################################


#### HOLD OUT EXTREME AGES ####

getwd()
source("../Estimators_Slow.R")

Proximity_Data <- read.delim(file="Data/nls.dat",sep="",header=FALSE)  %>% rename(
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
         f8 = ifelse(famed==8,1,0)) %>% 
  mutate(ID = 1:n())




holdoutextremes <- function(Data,noExtremes){
  
  HoldOutDataset <- Data %>% 
    as_tibble() %>%  
    arrange(age) %>% 
    head(-noExtremes/2) %>% 
    tail(-noExtremes/2)
  HeldOutSamples <- Data %>% 
    filter(ID %notin% HoldOutDataset$ID)
  
  
  
  n <- nrow(HoldOutDataset)
  #Target
  Y <- HoldOutDataset %>% select(lwage76)  %>% as.matrix
  #Included Exogenous:
  A_1 <- HoldOutDataset %>% select(Intercept) %>%  as.matrix
  #All Exogenous
  A <- HoldOutDataset %>% select(nearc4,age,agesq,Intercept) %>% as.matrix
  #Included Endogenous
  X <- HoldOutDataset %>% select(ed76,exp,expsq) %>%  as.matrix
  #Included (all)
  Z <- cbind(X,A_1)
  
  
  #OLS
  ols <- K_class(0,A,Z,Y,n)
  #IV
  iv <- K_class(1,A,Z,Y,n)
  #PULSE
  pulse <- PULSE(A,A_1,X,Y,p=0.05,N=10000,n)
  
  
  MSPE <- HeldOutSamples %>%  
    select(lwage76,Intercept,ed76,exp,expsq) %>% 
    mutate(resOLS = lwage76 - (ed76*ols["ed76",]+exp*ols["exp",]+expsq*ols["expsq",]+Intercept*ols["Intercept",]),
           resIV = lwage76 - (ed76*iv["ed76",]+exp*iv["exp",]+expsq*iv["expsq",]+Intercept*iv["Intercept",]),
           resPULSE = lwage76 - (ed76*pulse["ed76","lwage76"]+exp*pulse["exp","lwage76"]+expsq*pulse["expsq","lwage76"]+Intercept*pulse["Intercept","lwage76"])) %>% 
    summarize(OLS = sum(resOLS^2)/n(),
              IV = sum(resIV^2)/n(),
              PULSE = sum(resPULSE^2)/n(),
              ORDER = case_when(OLS <= IV & IV <= PULSE ~ "OLS<IV<PULSE",
                                OLS < PULSE & PULSE <= IV ~ "OLS<PULSE<IV",
                                OLS == PULSE & PULSE <= IV ~ "PULSE=OLS<IV",
                                IV <= PULSE & PULSE < OLS ~ "IV<PULSE<OLS",
                                IV <= PULSE & PULSE == OLS ~ "IV<PULSE=OLS",
                                IV <= OLS & OLS <= PULSE ~ "IV<OLS<PULSE",
                                IV <= OLS & OLS == PULSE ~ "IV<OLS=PULSE",
                                PULSE < OLS & OLS <= IV ~ "PULSE<OLS<IV",
                                PULSE == OLS & OLS <= IV ~ "PULSE=OLS<IV",
                                PULSE <= IV & IV <= OLS ~ "PULSE<IV<OLS"
              ),
              OLSmIV = OLS-IV,
              OLSmPULSE = OLS-PULSE,
              PULSEmIV = PULSE -IV,
              coef.OLS = ols["ed76",],
              coef.IV = iv["ed76",],
              coef.PULSE = pulse["ed76","lwage76"]) %>% 
    mutate(t=pulse["Intercept","t"],q=pulse["Intercept","q"])
  return(MSPE)
}

summarydat <-data.frame(n= c(250,500,750,1000,1250,1500,1750,2000,2250,2300)) %>%  #seq(500,2500 ,500)) %>% 
  rowwise() %>% 
  mutate(order=list(holdoutextremes(Proximity_Data,n)))%>%
  unnest(cols=c(order))

summarydat %>%  select(n,OLS,IV,PULSE,ORDER,coef.OLS,coef.IV,coef.PULSE) %>% 
  kbl(.,format="latex",digits=4)

#### HOLD OUT EXTREME IQs ####


##############################################################
# MODEL M3 # with kww test score and IQ as its instrument
##############################################################
holdoutextremes <- function(Data,noExtremes){
  
  HoldOutDataset <- Data %>% 
    as_tibble() %>%  
    arrange(iq) %>% 
    head(-noExtremes/2) %>% 
    tail(-noExtremes/2)
  HeldOutSamples <- Data %>% 
    filter(ID %notin% HoldOutDataset$ID)
  

  
  n <- nrow(HoldOutDataset)
  #Target
  Y <- HoldOutDataset %>% select(lwage76)  %>% as.matrix
  #Included Exogenous:
  A_1 <- HoldOutDataset %>% select(black,smsa76r,reg76r,smsa66r,reg662,reg663,reg664,reg665,reg666,reg667,reg668,reg669,Intercept,daded,momed,nodaded,nomomed ,f1,f2,f3,f4,f5,f6,f7,f8,momdad14,sinmom14) %>%  as.matrix
  #All Exogenous
  A <- HoldOutDataset %>% select(nearc4,age,agesq,iq,black,smsa76r,reg76r,smsa66r,reg662,reg663,reg664,reg665,reg666,reg667,reg668,reg669,Intercept,daded,momed,nodaded,nomomed ,f1,f2,f3,f4,f5,f6,f7,f8,momdad14,sinmom14) %>% as.matrix
  #Included Endogenous
  X <- HoldOutDataset %>% select(ed76,exp,expsq,kww) %>%  as.matrix
  #Included (all)
  Z <- cbind(X,A_1)
  
  #OLS
  ols <- K_class(0,A,Z,Y,n)
  #IV
  iv <- K_class(1,A,Z,Y,n)
  #PULSE
  pulse <- PULSE(A,A_1,X,Y,p=0.05,N=10000,n)
  
  
  MSPE <- HeldOutSamples %>%  
    select(lwage76,black,smsa76r,reg76r,smsa66r,reg662,reg663,reg664,reg665,reg666,reg667,reg668,reg669,Intercept,ed76,exp,expsq,kww) %>% 
    mutate(resOLS = lwage76 - (ed76*ols["ed76",]+exp*ols["exp",]+expsq*ols["expsq",]+black*ols["black",]+smsa76r*ols["smsa76r",]+reg76r*ols["reg76r",]+smsa66r*ols["smsa66r",]+reg662*ols["reg662",]+reg663*ols["reg663",]+reg664*ols["reg664",]+reg665*ols["reg665",]+reg666*ols["reg666",]+reg667*ols["reg667",]+reg668*ols["reg668",]+reg669*ols["reg669",]+Intercept*ols["Intercept",]+kww*ols["kww",]),
           resIV = lwage76 - (ed76*iv["ed76",]+exp*iv["exp",]+expsq*iv["expsq",]+black*iv["black",]+smsa76r*iv["smsa76r",]+reg76r*iv["reg76r",]+smsa66r*iv["smsa66r",]+reg662*iv["reg662",]+reg663*iv["reg663",]+reg664*iv["reg664",]+reg665*iv["reg665",]+reg666*iv["reg666",]+reg667*iv["reg667",]+reg668*iv["reg668",]+reg669*iv["reg669",]+Intercept*iv["Intercept",]+kww*iv["kww",]),
           resPULSE = lwage76 - (ed76*pulse["ed76","lwage76"]+exp*pulse["exp","lwage76"]+expsq*pulse["expsq","lwage76"]+black*pulse["black","lwage76"]+smsa76r*pulse["smsa76r","lwage76"]+reg76r*pulse["reg76r","lwage76"]+smsa66r*pulse["smsa66r","lwage76"]+reg662*pulse["reg662","lwage76"]+reg663*pulse["reg663","lwage76"]+reg664*pulse["reg664","lwage76"]+reg665*pulse["reg665","lwage76"]+reg666*pulse["reg666","lwage76"]+reg667*pulse["reg667","lwage76"]+reg668*pulse["reg668","lwage76"]+reg669*pulse["reg669","lwage76"]+Intercept*pulse["Intercept","lwage76"]++kww*pulse["kww","lwage76"])) %>% 
    summarize(OLS = sum(resOLS^2)/n(),
              IV = sum(resIV^2)/n(),
              PULSE = sum(resPULSE^2)/n(),
              ORDER = case_when(OLS <= IV & IV <= PULSE ~ "OLS<IV<PULSE",
                                OLS < PULSE & PULSE <= IV ~ "OLS<PULSE<IV",
                                OLS == PULSE & PULSE <= IV ~ "PULSE=OLS<IV",
                                IV <= PULSE & PULSE < OLS ~ "IV<PULSE<OLS",
                                IV <= PULSE & PULSE == OLS ~ "IV<PULSE=OLS",
                                IV <= OLS & OLS <= PULSE ~ "IV<OLS<PULSE",
                                IV <= OLS & OLS == PULSE ~ "IV<OLS=PULSE",
                                PULSE < OLS & OLS <= IV ~ "PULSE<OLS<IV",
                                PULSE == OLS & OLS <= IV ~ "PULSE=OLS<IV",
                                PULSE <= IV & IV <= OLS ~ "PULSE<IV<OLS"
              ),
              OLSmIV = OLS-IV,
              OLSmPULSE = OLS-PULSE,
              PULSEmIV = PULSE -IV,
              coef.OLS = ols["ed76",],
              coef.IV = iv["ed76",],
              coef.PULSE = pulse["ed76","lwage76"]) %>% 
    mutate(t=pulse["ed76","t"],q=pulse["ed76","q"],l=pulse["ed76","l"],k=pulse["ed76","k"])
  return(MSPE)
}

summarydat <-data.frame(n=c(100,200,300,400,500,600,700,800,900,1000,1100,1200,1300,1400,1500,1600,1700,1800)) %>% 
  rowwise() %>% 
  mutate(order=list(holdoutextremes(Proximity_Data %>%  filter(kww!=0,iq!=0),n)))%>%
  unnest(cols=c(order))

summarydat %>%  select(n,OLS,IV,PULSE,coef.OLS,coef.IV,coef.PULSE,ORDER) %>% 
  kbl(.,format="latex",digits=4)



########################################
### FINDING KAPPAS FOR HELDOUT DATA  ###
########################################

fff <-function(Data,no){
  HoldOutDataset <-  Data %>% 
    as_tibble() %>%  
    arrange(iq) %>% 
    head(-no/2) %>% 
    tail(-no/2)
  HeldOutSamples <- Data %>% 
    filter(ID %notin% HoldOutDataset$ID)
  
  d <- data.frame(HoldOutVar = HoldOutDataset %>% summarize(var(iq) ) %>% pull,HeldOutVar = HeldOutSamples %>%  summarize(var(iq)) %>% pull)
  return(d)
}
data.frame(no=seq(100,1800,100)) %>% rowwise() %>%  mutate( g= list(fff(Proximity_Data %>%  filter(kww!=0,iq!=0),no))) %>% unnest(cols=c(g)) %>% mutate(KappaForEquality = -(HoldOutVar-HeldOutVar)/HeldOutVar) %>% mutate(HeldOutOverHoldOut= HeldOutVar/HoldOutVar)





#### HOLD OUT RANDOM ####

set.seed(1)

holdoutrandom <- function(Data,percentage){
  HoldOutDataset <- Data %>% 
    as_tibble() %>%  
    sample_frac(percentage)
  HeldOutSamples <- Data %>% 
    filter(ID %notin% HoldOutDataset$ID)
  
  
  n <- nrow(HoldOutDataset)
  #Target
  Y <- HoldOutDataset %>% select(lwage76)  %>% as.matrix
  #Included Exogenous:
  A_1 <- HoldOutDataset %>% select(black,smsa76r,reg76r,smsa66r,reg662,reg663,reg664,reg665,reg666,reg667,reg668,reg669,Intercept) %>%  as.matrix
  #All Exogenous
  A <- HoldOutDataset %>% select(nearc4,age,agesq,black,smsa76r,reg76r,smsa66r,reg662,reg663,reg664,reg665,reg666,reg667,reg668,reg669,Intercept) %>% as.matrix
  #Included Endogenous
  X <- HoldOutDataset %>% select(ed76,exp,expsq) %>%  as.matrix
  #Included (all)
  Z <- cbind(X,A_1)
  
  
  #OLS
  ols <- K_class(0,A,Z,Y,n)
  #IV
  iv <- K_class(1,A,Z,Y,n)
  #PULSE
  pulse <- PULSE(A,A_1,X,Y,p=0.05,N=10000,n)
  
  #OLS
  ols <- K_class(0,A,Z,Y,n)
  #IV
  iv <- K_class(1,A,Z,Y,n)
  #PULSE
  pulse <- PULSE(A,A_1,X,Y,p=0.05,N=10000,n)
  
  
MSPE <- HeldOutSamples %>%  
    select(lwage76,black,smsa76r,reg76r,smsa66r,reg662,reg663,reg664,reg665,reg666,reg667,reg668,reg669,Intercept,ed76,exp,expsq) %>% 
    mutate(resOLS = lwage76 - (ed76*ols["ed76",]+exp*ols["exp",]+expsq*ols["expsq",]+black*ols["black",]+smsa76r*ols["smsa76r",]+reg76r*ols["reg76r",]+smsa66r*ols["smsa66r",]+reg662*ols["reg662",]+reg663*ols["reg663",]+reg664*ols["reg664",]+reg665*ols["reg665",]+reg666*ols["reg666",]+reg667*ols["reg667",]+reg668*ols["reg668",]+reg669*ols["reg669",]+Intercept*ols["Intercept",]),
           resIV = lwage76 - (ed76*iv["ed76",]+exp*iv["exp",]+expsq*iv["expsq",]+black*iv["black",]+smsa76r*iv["smsa76r",]+reg76r*iv["reg76r",]+smsa66r*iv["smsa66r",]+reg662*iv["reg662",]+reg663*iv["reg663",]+reg664*iv["reg664",]+reg665*iv["reg665",]+reg666*iv["reg666",]+reg667*iv["reg667",]+reg668*iv["reg668",]+reg669*iv["reg669",]+Intercept*iv["Intercept",]),
           resPULSE = lwage76 - (ed76*pulse["ed76","lwage76"]+exp*pulse["exp","lwage76"]+expsq*pulse["expsq","lwage76"]+black*pulse["black","lwage76"]+smsa76r*pulse["smsa76r","lwage76"]+reg76r*pulse["reg76r","lwage76"]+smsa66r*pulse["smsa66r","lwage76"]+reg662*pulse["reg662","lwage76"]+reg663*pulse["reg663","lwage76"]+reg664*pulse["reg664","lwage76"]+reg665*pulse["reg665","lwage76"]+reg666*pulse["reg666","lwage76"]+reg667*pulse["reg667","lwage76"]+reg668*pulse["reg668","lwage76"]+reg669*pulse["reg669","lwage76"]+Intercept*pulse["Intercept","lwage76"])) %>% 
    summarize(OLS = sum(resOLS^2)/n(),
              IV = sum(resIV^2)/n(),
              PULSE = sum(resPULSE^2)/n(),
              ORDER = case_when(pulse["ed76","l"]== 0 & PULSE <= IV ~ "PULSE=OLS<IV",
                                pulse["ed76","l"]== 0 & OLS <= IV ~ "PULSE=OLS<IV",
                                IV <= OLS & pulse["ed76","l"] == 0 ~ "IV<OLS=PULSE",
                                IV <= PULSE & pulse["ed76","l"] == 0 ~ "IV<PULSE=OLS",
                                OLS <= IV & IV <= PULSE ~ "OLS<IV<PULSE",
                                OLS < PULSE & PULSE <= IV ~ "OLS<PULSE<IV",
                                IV <= PULSE & PULSE < OLS ~ "IV<PULSE<OLS",
                                IV <= OLS & OLS <= PULSE ~ "IV<OLS<PULSE",
                                PULSE < OLS & OLS <= IV ~ "PULSE<OLS<IV",
                                PULSE <= IV & IV <= OLS ~ "PULSE<IV<OLS"
              ),
              OLSmIV = OLS-IV,
              OLSmPULSE = OLS-PULSE,
              PULSEmIV = PULSE -IV) %>% 
    mutate(t=pulse["ed76","t"],q=pulse["ed76","q"],l=pulse["ed76","l"],k=pulse["ed76","k"])
  return(MSPE)
}

norep <- 1000


summarydat <-data.frame(n=seq(1,norep,1)) %>% 
  rowwise() %>% 
  mutate(order=list(holdoutrandom(Proximity_Data  ,0.9)))


AverageStats <- summarydat %>%
  unnest(cols=c(order)) %>% 
  group_by(ORDER) %>% 
  summarise(count= n(),
            OLS = mean(OLS),
            PULSE = mean(PULSE),
            IV = mean(IV)  ) %>% 
  mutate(Percentage = 100*count/norep)


WorstCaseAverageStats <- bind_rows(
  summarydat %>%
    unnest(cols=c(order)) %>%
    filter(IV <= PULSE & IV <= OLS) %>% 
    arrange(desc(IV)) %>% 
    filter(row_number() <= n()*0.1) %>% 
    group_by(ORDER) %>% 
    summarise(OLS = mean(OLS),
              PULSE = mean(PULSE),
              IV = mean(IV))
  ,
  summarydat %>%
    unnest(cols=c(order)) %>%
    filter(OLS <= IV ) %>% 
    arrange(desc(IV)) %>% 
    filter(row_number() <= n()*0.1) %>% 
    group_by(ORDER) %>% 
    summarise(OLS = mean(OLS),
              PULSE = mean(PULSE),
              IV = mean(IV)))

left_join(AverageStats,WorstCaseAverageStats,by="ORDER") %>% 
  select(ORDER,Percentage,OLS.x,IV.x,PULSE.x,OLS.y,IV.y,PULSE.y) %>% 
  kbl(.,format="latex",digits=3)
