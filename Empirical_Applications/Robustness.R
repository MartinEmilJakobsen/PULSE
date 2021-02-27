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
library(kableExtra)

`%notin%` <- Negate(`%in%`)


source("../Estimators_Slow.R")


################################################
################################################
################################################
######## COLONIAL ROBUSTNESS EXPERIMENTS #######
################################################
################################################
################################################

COLONIAL_Data <- read.delim(file="Data/AJR01_COLONIAL.csv",sep=",",header=TRUE)
Selected_Data <- COLONIAL_Data  %>% filter(baseco == 1) %>% mutate_at(c("lat_abst","avexpr","logpgp95","logem4"),.funs=function(x){x-mean(x)}) %>% mutate(Intercept=1)

#############################################################################################
# Hold out extremes: MODEL M1 WITHOUT INTERCEPT: logpgp95 ~ avexpr  (INSTRUMENTS = logem4)#
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
  #Fuller4
  fuller4 <-  K_class(FULLER_k(4,A,A_1,X,Y,n,1),A,Z,Y,n)

  
  MSPE <- HeldOutSamples %>%  
    select(logpgp95,avexpr,Intercept) %>% 
    mutate(resOLS = logpgp95- (avexpr*ols["avexpr",]),
           resIV = logpgp95- (avexpr*iv["avexpr",]),
           resPULSE = logpgp95- (avexpr*pulse["avexpr","logpgp95"]),
           resFULLER4 = logpgp95- (avexpr*fuller4["avexpr","logpgp95"])) %>% 
    summarize(OLS = sum(resOLS^2)/n(),
              IV = sum(resIV^2)/n(),
              PULSE = sum(resPULSE^2)/n(),
              FULLER4 = sum(resFULLER4^2)/n(),
              WOLS = max(resOLS^2),
              WIV = max(resIV^2),
              WPULSE = max(resPULSE^2),
              WFULLER4 = max(resFULLER4^2),
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
              coef.PULSE = pulse["avexpr","logpgp95"],
              coef.FULLER4 = fuller4["avexpr","logpgp95"],
              kappa.FULLER4 = FULLER_k(4,A,A_1,X,Y,n,1)) %>% 
    mutate(t=pulse["avexpr","t"],q=pulse["avexpr","q"],l=pulse["avexpr","l"],k=pulse["avexpr","k"])
  return(MSPE)
}

summarydat <-data.frame(n=seq(4,32,2)) %>% 
  rowwise() %>% 
  mutate(order=list(holdoutextremes(Selected_Data,n)))%>%
  unnest(cols=c(order)) %>% 
  mutate(kappa.PULSE = l/(1+l))





summarydat %>%  select(n,coef.OLS,coef.IV,coef.PULSE,coef.FULLER4,kappa.PULSE,kappa.FULLER4,OLS,IV,PULSE,FULLER4) %>% 
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
  #fuller4
  fuller4 <-  K_class(FULLER_k(4,A,A_1,X,Y,n,1),A,Z,Y,n)
  
  
  MSPE <- HeldOutSamples %>%  
    select(logpgp95,avexpr,Intercept) %>% 
    mutate(resOLS = logpgp95- (avexpr*ols["avexpr",]),
           resIV = logpgp95- (avexpr*iv["avexpr",]),
           resPULSE = logpgp95- (avexpr*pulse["avexpr","logpgp95"]),
           resFULLER4 = logpgp95- (avexpr*fuller4["avexpr","logpgp95"])) %>% 
    summarize(OLS = sum(resOLS^2)/n(),
              IV = sum(resIV^2)/n(),
              PULSE = sum(resPULSE^2)/n(),
              FULLER4 = sum(resFULLER4^2)/n()) %>% 
    mutate(t=pulse["avexpr","t"],q=pulse["avexpr","q"],l=pulse["avexpr","l"],k=pulse["avexpr","k"])
  return(MSPE)
}

norep <- 1000

summarydat <-data.frame(n=seq(1,norep,1)) %>% 
  rowwise() %>% 
  mutate(dat=list(holdoutrandom(Selected_Data,0.9)))

#OLD table

MSPEORDERS <- summarydat %>%
  unnest(cols=c(dat)) %>%  
  gather(key="Type",value="MSPE",c(OLS,IV,PULSE,FULLER4)) %>% 
  arrange(MSPE) %>% 
  group_by(n) %>% 
  summarise(n=max(n),
            t=max(t),
            q=max(q),
            l=max(l),
            k=max(k),
            ORDER = paste0(Type,collapse="<")) 
  
  summarydat %>%
  unnest(cols=c(dat)) %>% 
  left_join(.,MSPEORDERS %>% select(n,ORDER),by="n") %>% 
  group_by(ORDER) %>% 
  summarise(count= n(),
            WMSPE.OLS = max(OLS),
            WMSPE.iv = max(IV),
            WMSPE.pulse = max(PULSE),
            WMSPE.fuller4 = max(FULLER4),
            OLS = mean(OLS),
            PULSE = mean(PULSE),
            IV = mean(IV),
            FULLER4 = mean(FULLER4)
            ) %>% 
  mutate(Percentage = 100*count/norep) %>% 
  select(ORDER,Percentage,OLS,IV,PULSE,FULLER4,WMSPE.OLS,WMSPE.iv,WMSPE.pulse,WMSPE.fuller4) %>% 
  kbl(.,format="latex",digits=3)
  
  
  
#Table in paper: 
summarydat %>%
    unnest(cols=c(dat))  %>% 
    mutate(OT = I(OLS <= IV),
           OP = I(OLS <=PULSE),
           OF = I(OLS<=FULLER4),
           TO = I(IV<=OLS),
           TP = I(IV<=PULSE),
           TF = I(IV<=FULLER4),
           PO = I(PULSE <= OLS),
           PT = I(PULSE<= IV),
           PF = I(PULSE<=FULLER4),
           FO = I(FULLER4<=OLS),
           FT = I(FULLER4<=IV),
           FP = I(FULLER4<=PULSE)) %>% 
    summarise(OP= sum(OP)/n(),
              OF= sum(OF)/n(),
              OT= sum(OT)/n(),
              PO= sum(PO)/n(),
              PF= sum(PF)/n(),
              PT= sum(PT)/n(),
              FO= sum(FO)/n(),
              FP= sum(FP)/n(),
              FT= sum(FT)/n(),
              TO= sum(TO)/n(),
              TP= sum(TP)/n(),
              TF= sum(TF)/n()
              ) %>%  
    unlist(., use.names=FALSE) %>% 
    sapply(.,function(x){x*100}) %>% 
    matrix(.,ncol=3,byrow=TRUE) %>% 
    as.data.frame(.,row.names = c("OLS","TSLS","PULSE","FUL")) %>%
    mutate_at(c("V1","V2","V3"),.funs=function(x){round(x,digits=1) %>% paste0(.,'%')})  %>% 
    kbl(,format="latex",digits=3,row.names = TRUE)
  

################################################
################################################
################################################
########   ANGRIST AND KRUGER 1991 - QOB #######
################################################
################################################
################################################


QOB_Data <- read.delim(file="Data/AK91_QOB.txt",sep=" ",header=FALSE) %>% 
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

#################
#Hold out random#
#################


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
  #fuller4
  fuller4 <-  K_class(FULLER_k(4,A,A_1,X,Y,n,1),A,Z,Y,n)
  
  MSPE <- HeldOutSamples %>%  
    select(LWKLYWGE,EDUC,YR0 ,  YR1 ,  YR2 ,  YR3 ,  YR4 ,  YR5  , YR6  , YR7 ,  YR8,   YR9) %>% 
    mutate(
      resOLS = LWKLYWGE- (EDUC*ols["EDUC",]+ YR0*ols["YR0",]+YR1*ols["YR1",]+YR2*ols["YR2",]+YR3*ols["YR3",]+YR4*ols["YR4",]+YR5*ols["YR5",]+YR6*ols["YR6",]+YR7*ols["YR7",]+YR8*ols["YR8",]+YR9*ols["YR9",]),
      resIV = LWKLYWGE- (EDUC*iv["EDUC",]+ YR0*iv["YR0",]+YR1*iv["YR1",]+YR2*iv["YR2",]+YR3*iv["YR3",]+YR4*iv["YR4",]+YR5*iv["YR5",]+YR6*iv["YR6",]+YR7*iv["YR7",]+YR8*iv["YR8",]+YR9*iv["YR9",]),
      resPULSE = LWKLYWGE- (EDUC*pulse["EDUC","LWKLYWGE"]+ YR0*pulse["YR0","LWKLYWGE"]+YR1*pulse["YR1","LWKLYWGE"]+YR2*pulse["YR2","LWKLYWGE"]+YR3*pulse["YR3","LWKLYWGE"]+YR4*pulse["YR4","LWKLYWGE"]+YR5*pulse["YR5","LWKLYWGE"]+YR6*pulse["YR6","LWKLYWGE"]+YR7*pulse["YR7","LWKLYWGE"]+YR8*pulse["YR8","LWKLYWGE"]+YR9*pulse["YR9","LWKLYWGE"]),
      resFULLER4 = LWKLYWGE- (EDUC*fuller4["EDUC","LWKLYWGE"]+ YR0*fuller4["YR0","LWKLYWGE"]+YR1*fuller4["YR1","LWKLYWGE"]+YR2*fuller4["YR2","LWKLYWGE"]+YR3*fuller4["YR3","LWKLYWGE"]+YR4*fuller4["YR4","LWKLYWGE"]+YR5*fuller4["YR5","LWKLYWGE"]+YR6*fuller4["YR6","LWKLYWGE"]+YR7*fuller4["YR7","LWKLYWGE"]+YR8*fuller4["YR8","LWKLYWGE"]+YR9*fuller4["YR9","LWKLYWGE"])) %>% 
    summarize(OLS = sum(resOLS^2)/n(),
              IV = sum(resIV^2)/n(),
              PULSE = sum(resPULSE^2)/n(),
              FULLER4 = sum(resFULLER4^2)/n()) %>% 
    mutate(t=pulse["EDUC","t"],q=pulse["EDUC","q"],l=pulse["EDUC","l"],k=pulse["EDUC","k"])
  return(MSPE)
}

norep <- 1000

Selected_Data_subset <- Selected_Data_Centered %>% sample_frac(0.01)

summarydat <-data.frame(n=seq(1,norep,1)) %>% 
  rowwise() %>% 
  mutate(order=list(holdoutrandom(Selected_Data_subset  ,0.9)))


MSPEORDERS <- summarydat %>%
  unnest(cols=c(order)) %>%  
  gather(key="Type",value="MSPE",c(OLS,IV,PULSE,FULLER4)) %>% 
  arrange(MSPE) %>% 
  group_by(n) %>% 
  summarise(n=max(n),
            t=max(t),
            q=max(q),
            l=max(l),
            k=max(k),
            ORDER = paste0(Type,collapse="<")) %>% 
  mutate(ORDER=case_when(ORDER == "PULSE<OLS<IV<FULLER4" ~ "PULSE=OLS<IV<FULLER4",
                         ORDER == "OLS<PULSE<IV<FULLER4" ~ "PULSE=OLS<IV<FULLER4",
                         ORDER == "IV<FULLER4<PULSE<OLS" ~ "IV<FULLER4<PULSE=OLS",
                         ORDER == "IV<OLS<PULSE<FULLER4" ~ "IV<PULSE=OLS<FULLER4",
                         ORDER == "IV<PULSE<OLS<FULLER4" ~ "IV<PULSE=OLS<FULLER4",
                         ORDER == "FULLER4<IV<OLS<PULSE" ~ "FULLER4<IV<PULSE=OLS",
                         ORDER == "IV<FULLER4<OLS<PULSE" ~ "IV<FULLER4<PULSE=OLS",
                         
                         ))


summarydat %>%
  unnest(cols=c(order)) %>% 
  left_join(.,MSPEORDERS %>% select(n,ORDER),by="n") %>% 
  ungroup %>% 
  group_by(ORDER) %>% 
  summarise(count= n(),
            WMSPE.OLS = quantile(OLS,probs=0.9),
            WMSPE.iv = quantile(IV,probs=0.9),
            WMSPE.pulse = quantile(PULSE,probs=0.9),
            WMSPE.fuller4 = quantile(FULLER4,probs=0.9),
            OLS = mean(OLS),
            PULSE = mean(PULSE),
            IV = mean(IV),
            FULLER4 = mean(FULLER4)
  ) %>% 
  mutate(Percentage = 100*count/norep) %>% 
  select(ORDER,Percentage,OLS,IV,PULSE,FULLER4,WMSPE.OLS,WMSPE.iv,WMSPE.pulse,WMSPE.fuller4) %>% 
  kbl(.,format="latex",digits=3)

################################################
################################################
################################################
######   CARD 1995 - Proximity to college ######
################################################
################################################
################################################


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
         f8 = ifelse(famed==8,1,0)) %>% 
  mutate(ID = 1:n())

#################
#Hold out random#
#################
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
  #fuller4
  fuller4 <-  K_class(FULLER_k(4,A,A_1,X,Y,n,1),A,Z,Y,n)
 
  
  
MSPE <- HeldOutSamples %>%  
    select(lwage76,black,smsa76r,reg76r,smsa66r,reg662,reg663,reg664,reg665,reg666,reg667,reg668,reg669,Intercept,ed76,exp,expsq) %>% 
    mutate(resOLS = lwage76 - (ed76*ols["ed76",]+exp*ols["exp",]+expsq*ols["expsq",]+black*ols["black",]+smsa76r*ols["smsa76r",]+reg76r*ols["reg76r",]+smsa66r*ols["smsa66r",]+reg662*ols["reg662",]+reg663*ols["reg663",]+reg664*ols["reg664",]+reg665*ols["reg665",]+reg666*ols["reg666",]+reg667*ols["reg667",]+reg668*ols["reg668",]+reg669*ols["reg669",]+Intercept*ols["Intercept",]),
           resIV = lwage76 - (ed76*iv["ed76",]+exp*iv["exp",]+expsq*iv["expsq",]+black*iv["black",]+smsa76r*iv["smsa76r",]+reg76r*iv["reg76r",]+smsa66r*iv["smsa66r",]+reg662*iv["reg662",]+reg663*iv["reg663",]+reg664*iv["reg664",]+reg665*iv["reg665",]+reg666*iv["reg666",]+reg667*iv["reg667",]+reg668*iv["reg668",]+reg669*iv["reg669",]+Intercept*iv["Intercept",]),
           resPULSE = lwage76 - (ed76*pulse["ed76","lwage76"]+exp*pulse["exp","lwage76"]+expsq*pulse["expsq","lwage76"]+black*pulse["black","lwage76"]+smsa76r*pulse["smsa76r","lwage76"]+reg76r*pulse["reg76r","lwage76"]+smsa66r*pulse["smsa66r","lwage76"]+reg662*pulse["reg662","lwage76"]+reg663*pulse["reg663","lwage76"]+reg664*pulse["reg664","lwage76"]+reg665*pulse["reg665","lwage76"]+reg666*pulse["reg666","lwage76"]+reg667*pulse["reg667","lwage76"]+reg668*pulse["reg668","lwage76"]+reg669*pulse["reg669","lwage76"]+Intercept*pulse["Intercept","lwage76"]),
           resFULLER4 = lwage76 - (ed76*fuller4["ed76","lwage76"]+exp*fuller4["exp","lwage76"]+expsq*fuller4["expsq","lwage76"]+black*fuller4["black","lwage76"]+smsa76r*fuller4["smsa76r","lwage76"]+reg76r*fuller4["reg76r","lwage76"]+smsa66r*fuller4["smsa66r","lwage76"]+reg662*fuller4["reg662","lwage76"]+reg663*fuller4["reg663","lwage76"]+reg664*fuller4["reg664","lwage76"]+reg665*fuller4["reg665","lwage76"]+reg666*fuller4["reg666","lwage76"]+reg667*fuller4["reg667","lwage76"]+reg668*fuller4["reg668","lwage76"]+reg669*fuller4["reg669","lwage76"]+Intercept*fuller4["Intercept","lwage76"])) %>% 
    summarize(OLS = sum(resOLS^2)/n(),
              IV = sum(resIV^2)/n(),
              PULSE = sum(resPULSE^2)/n(),
              FULLER4 = sum(resFULLER4^2)/n()) %>% 
    mutate(t=pulse["ed76","t"],q=pulse["ed76","q"],l=pulse["ed76","l"],k=pulse["ed76","k"])
  return(MSPE)
}

norep <- 1000


summarydat <-data.frame(n=seq(1,norep,1)) %>% 
  rowwise() %>% 
  mutate(order=list(holdoutrandom(Proximity_Data  ,0.9)))


MSPEORDERS <- summarydat %>%
  unnest(cols=c(order)) %>%  
  gather(key="Type",value="MSPE",c(OLS,IV,PULSE,FULLER4)) %>% 
  arrange(MSPE) %>% 
  group_by(n) %>% 
  summarise(n=max(n),
            t=max(t),
            q=max(q),
            l=max(l),
            k=max(k),
            ORDER = paste0(Type,collapse="<")) %>% 
  mutate(ORDER=case_when(ORDER == "PULSE<OLS<FULLER4<IV" ~ "PULSE=OLS<FULLER4<IV",
                         ORDER == "PULSE<OLS<IV<FULLER4" ~ "PULSE=OLS<IV<FULLER4",
                         ORDER == "OLS<PULSE<FULLER4<IV" ~ "PULSE=OLS<FULLER4<IV",
                         ORDER == "IV<FULLER4<OLS<PULSE" ~ "IV<FULLER4<PULSE=OLS",
                         ORDER == "OLS<PULSE<IV<FULLER4" ~ "PULSE=OLS<IV<FULLER4",
                         ORDER == "FULLER4<OLS<PULSE<IV" ~ "FULLER4<PULSE=OLS<IV",
                         ORDER == "FULLER4<PULSE<OLS<IV" ~ "FULLER4<PULSE=OLS<IV",
                         ORDER == "FULLER4<IV<PULSE<OLS" ~ "FULLER4<IV<PULSE=OLS",
                         ORDER == "IV<FULLER4<PULSE<OLS" ~ "IV<FULLER4<PULSE=OLS",
                         ORDER == "FULLER4<IV<OLS<PULSE" ~ "FULLER4<IV<PULSE=OLS",
                         ORDER == "IV<OLS<PULSE<FULLER4" ~ "IV<PULSE=OLS<FULLER4",
                         
  ))


summarydat %>%
  unnest(cols=c(order)) %>% 
  left_join(.,MSPEORDERS %>% select(n,ORDER),by="n") %>% 
  ungroup %>% 
  group_by(ORDER) %>% 
  summarise(count= n(),
            WMSPE.OLS = quantile(OLS,probs=0.9),
            WMSPE.iv = quantile(IV,probs=0.9),
            WMSPE.pulse = quantile(PULSE,probs=0.9),
            WMSPE.fuller4 = quantile(FULLER4,probs=0.9),
            OLS = mean(OLS),
            PULSE = mean(PULSE),
            IV = mean(IV),
            FULLER4 = mean(FULLER4)
  ) %>% 
  mutate(Percentage = 100*count/norep) %>% 
  select(ORDER,Percentage,OLS,IV,PULSE,FULLER4,WMSPE.OLS,WMSPE.iv,WMSPE.pulse,WMSPE.fuller4) %>% 
  kbl(.,format="latex",digits=3)
