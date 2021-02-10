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


Selected_Data %>% arrange(logem4) 
#############################################################################################
# MODEL M1 WITH INTERCEPT: logpgp95 ~ avexpr  (INSTRUMENTS = logem4 )#
#############################################################################################
holdoutextremes <- function(Data,noExtremes){
  
  
  # Data <- Selected_Data
  # noExtremes <- 2
  # instrumentordering <- "logem4"
  
  
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
  A_1 <- HoldOutDataset %>% select(Intercept) %>%  as.matrix
  #All Exogenous
  A <- HoldOutDataset %>% select(logem4,Intercept) %>% as.matrix
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
    mutate(resOLS = logpgp95- (avexpr*ols["avexpr",]+ols["Intercept",]),
           resIV = logpgp95- (avexpr*iv["avexpr",]+iv["Intercept",]),
           resPULSE = logpgp95- (avexpr*pulse["avexpr","logpgp95"]+pulse["Intercept","logpgp95"])) %>% 
    summarize(OLS = round(sum(resOLS^2)/n(),3),
              IV = round(sum(resIV^2)/n(),3),
              PULSE = round(sum(resPULSE^2)/n(),3),
              ORDER = case_when(OLS <= IV & IV <= PULSE ~ "OLS<IV<PULSE",
                                OLS < PULSE & PULSE <= IV ~ "OLS<PULSE<IV",
                                OLS == PULSE & PULSE <= IV ~ "OLS=PULSE<IV",
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
              PULSEmIV = PULSE -IV)
  return(MSPE)
}

summarydat <-data.frame(n=seq(1,40,2)) %>% 
  rowwise() %>% 
  mutate(order=list(holdoutextremes(Selected_Data,n)))%>%
  unnest(cols=c(order))

summarydat %>%  print(n=60)


#############################################################################################
# MODEL M1 WITHOUT INTERCEPT: logpgp95 ~ avexpr  (INSTRUMENTS = logem4)#
#############################################################################################
holdoutextremes <- function(Data,noExtremes){
  
  
  # Data <- Selected_Data
  # noExtremes <- 2
  # instrumentordering <- "logem4"
  
  
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
    summarize(OLS = round(sum(resOLS^2)/n(),3),
              IV = round(sum(resIV^2)/n(),3),
              PULSE = round(sum(resPULSE^2)/n(),3),
              ORDER = case_when(OLS <= IV & IV <= PULSE ~ "OLS<IV<PULSE",
                                OLS < PULSE & PULSE <= IV ~ "OLS<PULSE<IV",
                                OLS == PULSE & PULSE <= IV ~ "OLS=PULSE<IV",
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
              PULSEmIV = PULSE -IV)
  return(MSPE)
}

summarydat <-data.frame(n=seq(2,40,2)) %>% 
  rowwise() %>% 
  mutate(order=list(holdoutextremes(Selected_Data,n)))%>%
  unnest(cols=c(order))

summarydat %>%  print(n=60)


#############################################################################################
# MODEL M1 WITHOUT INTERCEPT: logpgp95 ~ avexpr  (INSTRUMENTS = logem4)
#############################################################################################

# # TRAIN w/o Neoeurope TEST Neoeurpe
  HoldOutDataset <- Selected_Data %>% 
    as_tibble() %>%  
    filter(rich4!=1)
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
    select(logpgp95,avexpr) %>% 
    mutate(resOLS = logpgp95- (avexpr*ols["avexpr",]),
           resIV = logpgp95- (avexpr*iv["avexpr",]),
           resPULSE = logpgp95- (avexpr*pulse["avexpr","logpgp95"])) %>% 
    summarize(OLS = round(sum(resOLS^2)/n(),3),
              IV = round(sum(resIV^2)/n(),3),
              PULSE = round(sum(resPULSE^2)/n(),3),
              ORDER = case_when(OLS <= IV & IV <= PULSE ~ "OLS<IV<PULSE",
                                OLS < PULSE & PULSE <= IV ~ "OLS<PULSE<IV",
                                OLS == PULSE & PULSE <= IV ~ "OLS=PULSE<IV",
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
              PULSEmIV = PULSE -IV)

# TRAIN w/o africa TEST = africa

  HoldOutDataset <- Selected_Data %>% 
    as_tibble() %>%  
    filter(africa!=1)
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
    select(logpgp95,avexpr) %>% 
    mutate(resOLS = logpgp95- (avexpr*ols["avexpr",]),
           resIV = logpgp95- (avexpr*iv["avexpr",]),
           resPULSE = logpgp95- (avexpr*pulse["avexpr","logpgp95"])) %>% 
    summarize(OLS = round(sum(resOLS^2)/n(),3),
              IV = round(sum(resIV^2)/n(),3),
              PULSE = round(sum(resPULSE^2)/n(),3),
              ORDER = case_when(OLS <= IV & IV <= PULSE ~ "OLS<IV<PULSE",
                                OLS < PULSE & PULSE <= IV ~ "OLS<PULSE<IV",
                                OLS == PULSE & PULSE <= IV ~ "OLS=PULSE<IV",
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
              PULSEmIV = PULSE -IV)

  
MSPE  


# TRAIN w/o africa TEST = africa

HoldOutDataset <- Selected_Data %>% 
  as_tibble() %>%  
  filter(asia!=1)
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
  select(logpgp95,avexpr) %>% 
  mutate(resOLS = logpgp95- (avexpr*ols["avexpr",]),
         resIV = logpgp95- (avexpr*iv["avexpr",]),
         resPULSE = logpgp95- (avexpr*pulse["avexpr","logpgp95"])) %>% 
  summarize(OLS = round(sum(resOLS^2)/n(),3),
            IV = round(sum(resIV^2)/n(),3),
            PULSE = round(sum(resPULSE^2)/n(),3),
            ORDER = case_when(OLS <= IV & IV <= PULSE ~ "OLS<IV<PULSE",
                              OLS < PULSE & PULSE <= IV ~ "OLS<PULSE<IV",
                              OLS == PULSE & PULSE <= IV ~ "OLS=PULSE<IV",
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
            PULSEmIV = PULSE -IV)


MSPE  




#######################################################################
# WITH LAT_ABST


##############################
#Hold out extremes experiment#
##############################

holdoutextremes <- function(Data,noExtremes,instrumentordering){
  
  #  Data <- Selected_Data
  #  noExtremes <- 2 
  #  instrumentordering <- "logem4"
  
  HoldOutDataset <- Data %>% 
    as_tibble() %>%  
    arrange(!!sym(instrumentordering)) %>% 
    head(-noExtremes/2) %>% 
    tail(-noExtremes/2)
  HeldOutSamples <- Data %>% 
    filter(shortnam %notin% HoldOutDataset$shortnam)
  
  n <- nrow(HoldOutDataset)
  #Target
  Y <- HoldOutDataset %>% select(logpgp95) %>%  as.matrix
  #Included Exogenous:
  A_1 <- HoldOutDataset %>% select(lat_abst) %>%  as.matrix
  #All Exogenous
  A <- HoldOutDataset %>% select(logem4,lat_abst) %>% as.matrix
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
    select(logpgp95,avexpr,logem4,lat_abst) %>% 
    mutate(resOLS = logpgp95- (avexpr*ols["avexpr",]+lat_abst*ols["lat_abst",]),
           resIV = logpgp95- (avexpr*iv["avexpr",]+lat_abst*iv["lat_abst",]),
           resPULSE = logpgp95- (avexpr*pulse["avexpr","logpgp95"]+lat_abst*pulse["lat_abst","logpgp95"])) %>% 
    summarize(OLS = round(sum(resOLS^2)/n(),3),
              IV = round(sum(resIV^2)/n(),3),
              PULSE = round(sum(resPULSE^2)/n(),3),
              ORDER = case_when(OLS <= IV & IV <= PULSE ~ "OLS<IV<PULSE",
                                OLS < PULSE & PULSE <= IV ~ "OLS<PULSE<IV",
                                OLS == PULSE & PULSE <= IV ~ "OLS=PULSE<IV",
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
              Coef.ols = ols["avexpr",],
              Coef.iv = iv["avexpr",],
              Coef.pulse = pulse["avexpr","logpgp95"])
  return(MSPE)
}

summarydat <-data.frame(n=seq(2,40,2)) %>% 
  rowwise() %>% 
  mutate(order=list(holdoutextremes(Selected_Data ,n,"logem4")))%>%
  unnest(cols=c(order))

summarydat %>%  print(n=60)


summarydat <-data.frame(n=seq(2,40,2)) %>% 
  rowwise() %>% 
  mutate(order=list(holdoutextremes(Selected_Data,n,"lat_abst")))%>%
  unnest(cols=c(order))

summarydat %>%  print(n=60)



#######################################################################
# WITH LAT_ABST AND CONTINENT


##############################
#Hold out extremes experiment#
##############################

holdoutextremes <- function(Data,noExtremes,instrumentordering){
  
  HoldOutDataset <- Data %>% 
    as_tibble() %>%  
    arrange(!!sym(instrumentordering)) %>% 
    head(-noExtremes/2) %>% 
    tail(-noExtremes/2)
  HeldOutSamples <- Data %>% 
    filter(shortnam %notin% HoldOutDataset$shortnam)
  
  n <- nrow(HoldOutDataset)
  #Target
  Y <- HoldOutDataset %>% select(logpgp95) %>%  as.matrix
  #Included Exogenous:
  A_1 <- HoldOutDataset %>% select(lat_abst,africa,asia) %>%  as.matrix
  #All Exogenous
  A <- HoldOutDataset %>% select(logem4,lat_abst,africa,asia) %>% as.matrix
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
    select(logpgp95,avexpr,logem4,lat_abst,africa,asia) %>% 
    mutate(resOLS = logpgp95- (avexpr*ols["avexpr",]+lat_abst*ols["lat_abst",]+africa*ols["africa",]+asia*ols["asia",]),
           resIV = logpgp95- (avexpr*iv["avexpr",]+lat_abst*iv["lat_abst",]+africa*iv["africa",]+asia*iv["asia",]),
           resPULSE = logpgp95- (avexpr*pulse["avexpr","logpgp95"]+lat_abst*pulse["lat_abst","logpgp95"]+africa*pulse["africa","logpgp95"]+asia*pulse["asia","logpgp95"])) %>% 
    summarize(OLS = round(sum(resOLS^2)/n(),3),
              IV = round(sum(resIV^2)/n(),3),
              PULSE = round(sum(resPULSE^2)/n(),3),
              ORDER = case_when(OLS <= IV & IV <= PULSE ~ "OLS<IV<PULSE",
                                OLS < PULSE & PULSE <= IV ~ "OLS<PULSE<IV",
                                OLS == PULSE & PULSE <= IV ~ "OLS=PULSE<IV",
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
              Coef.ols = ols["avexpr",],
              Coef.iv = iv["avexpr",],
              Coef.pulse = pulse["avexpr","logpgp95"])
  return(MSPE)
}

summarydat <-data.frame(n=seq(2,40,2)) %>% 
  rowwise() %>% 
  mutate(order=list(holdoutextremes(Selected_Data,n,"logem4")))%>%
  unnest(cols=c(order))

summarydat %>%  print(n=60)


summarydat <-data.frame(n=seq(2,40,2)) %>% 
  rowwise() %>% 
  mutate(order=list(holdoutextremes(Selected_Data,n,"lat_abst")))%>%
  unnest(cols=c(order))

summarydat %>%  print(n=60)





#########################################################
# HOW OUT 20 extreme samples and ITTERATIVELY ADD TO HELD OUT SAMPLES TO SEE MSPE 
##########################################################

HoldOutDataset <- Data %>% 
  as_tibble() %>%  
  arrange(logem4) %>% 
  head(-10/2) %>% 
  tail(-10/2)
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


HeldOutSamples %>%  
  select(logpgp95,avexpr,logem4) %>% 
  mutate(resOLS = logpgp95- (avexpr*ols["avexpr",]),
         resIV = logpgp95- (avexpr*iv["avexpr",]),
         resPULSE = logpgp95- (avexpr*pulse["avexpr","logpgp95"])) %>% 
  mutate(n= 1:n()) %>% 
  mutate(OLS = round(cumsum(resOLS^2)/n,3),
         IV = round(cumsum(resIV^2)/n,3),
         PULSE = round(cumsum(resPULSE^2)/n,3))

HeldOutSamples[which(HeldOutSamples$shortnam=="AUS"),]

targetorder <- rev(c("AUS","MLI","NZL","NGA","HKG","GMB","USA","TGO","ZAF","GHA"))

HeldOutSamples[match(targetorder,HeldOutSamples$shortnam),] %>% 
  select(logpgp95,avexpr,logem4) %>% 
  mutate(resOLS = logpgp95- (avexpr*ols["avexpr",]),
         resIV = logpgp95- (avexpr*iv["avexpr",]),
         resPULSE = logpgp95- (avexpr*pulse["avexpr","logpgp95"])) %>% 
  mutate(n= 1:n()) %>% 
  mutate(OLS = round(cumsum(resOLS^2)/n,3),
         IV = round(cumsum(resIV^2)/n,3),
         PULSE = round(cumsum(resPULSE^2)/n,3)) %>% 
  mutate( ORDER = case_when(OLS <= IV & IV <= PULSE ~ "OLS<IV<PULSE",
                            OLS < PULSE & PULSE <= IV ~ "OLS<PULSE<IV",
                            OLS == PULSE & PULSE <= IV ~ "OLS=PULSE<IV",
                            IV <= PULSE & PULSE < OLS ~ "IV<PULSE<OLS",
                            IV <= PULSE & PULSE == OLS ~ "IV<PULSE=OLS",
                            IV <= OLS & OLS <= PULSE ~ "IV<OLS<PULSE",
                            IV <= OLS & OLS == PULSE ~ "IV<OLS=PULSE",
                            PULSE < OLS & OLS <= IV ~ "PULSE<OLS<IV",
                            PULSE == OLS & OLS <= IV ~ "PULSE=OLS<IV",
                            PULSE <= IV & IV <= OLS ~ "PULSE<IV<OLS")) %>% 
  mutate(CumVarLogem4 = cumsum(logem4^2)/n)

data.frame(n=1:10) %>% 
  mutate(nsq=n^2,cumsum(n^2))
MSPE <- HeldOutSamples %>%  
  select(logpgp95,avexpr,logem4) %>% 
  mutate(resOLS = logpgp95- (avexpr*ols["avexpr",]),
         resIV = logpgp95- (avexpr*iv["avexpr",]),
         resPULSE = logpgp95- (avexpr*pulse["avexpr","logpgp95"])) %>% 
  summarize(OLS = round(sum(resOLS^2)/n(),3),
            IV = round(sum(resIV^2)/n(),3),
            PULSE = round(sum(resPULSE^2)/n(),3),
            ORDER = case_when(OLS <= IV & IV <= PULSE ~ "OLS<IV<PULSE",
                              OLS < PULSE & PULSE <= IV ~ "OLS<PULSE<IV",
                              OLS == PULSE & PULSE <= IV ~ "OLS=PULSE<IV",
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
            PULSEmIV = PULSE -IV)



#################################################################################################

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
  A_1 <- HoldOutDataset %>% select(Intercept) %>%  as.matrix
  #All Exogenous
  A <- HoldOutDataset %>% select(logem4,Intercept) %>% as.matrix
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
    mutate(resOLS = logpgp95- (avexpr*ols["avexpr",]+ols["Intercept",]),
           resIV = logpgp95- (avexpr*iv["avexpr",]+iv["Intercept",]),
           resPULSE = logpgp95- (avexpr*pulse["avexpr",]+pulse["Intercept",])) %>% 
    summarize(OLS = round(sum(resOLS^2)/n(),3),
              IV = round(sum(resIV^2)/n(),3),
              PULSE = round(sum(resPULSE^2)/n(),3),
              ORDER = case_when(OLS <= IV & IV <= PULSE ~ "OLS<IV<PULSE",
                                OLS < PULSE & PULSE <= IV ~ "OLS<PULSE<IV",
                                OLS == PULSE & PULSE <= IV ~ "OLS=PULSE<IV",
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
              PULSEmIV = PULSE -IV)
  return(MSPE)
}
summarydat <-data.frame(n=seq(1,1000,1)) %>% 
  rowwise() %>% 
  mutate(order=list(holdoutrandom(Selected_Data,0.9)))

summarydat %>%
  unnest(cols=c(order)) %>% 
  group_by(ORDER) %>% 
  summarise(count= n(),
            OLS = mean(OLS),
            PULSE = mean(PULSE),
            IV = mean(IV)) %>% 
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
                                ORDER == "PULSE<IV<OLS" ~ OLS-IV))


#2 #PULSE < OLS < IV

set.seed(10)
HoldOutDataset <- Selected_Data %>% as_tibble() %>%  sample_frac(0.90)
HeldOutSamples <- Selected_Data %>% filter(shortnam %notin% HoldOutDataset$shortnam)

HoldOutDataset
HeldOutSamples

n <- nrow(HoldOutDataset)
#Target
Y <- HoldOutDataset %>% select(logpgp95) %>%  as.matrix
#Included Exogenous:
A_1 <- HoldOutDataset %>% select(Intercept) %>%  as.matrix
#All Exogenous
A <- HoldOutDataset %>% select(logem4,Intercept) %>% as.matrix
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

ols
iv
pulse

HeldOutSamples %>%  
  select(logpgp95,avexpr,Intercept) %>% 
  mutate(res = logpgp95- (avexpr*ols["avexpr",]+ols["Intercept",]) ) %>% 
  summarize(MSPE = sum(res^2))

HeldOutSamples %>%  
  select(logpgp95,avexpr,Intercept) %>% 
  mutate(res = logpgp95- (avexpr*iv["avexpr",]+iv["Intercept",]) ) %>% 
  summarize(MSPE = sum(res^2))

HeldOutSamples %>%  
  select(logpgp95,avexpr,Intercept) %>% 
  mutate(res = logpgp95- (avexpr*pulse["avexpr",]+pulse["Intercept",]) ) %>% 
  summarize(MSPE = sum(res^2))


#3 #OLS < PULSE < IV

set.seed(100)
HoldOutDataset <- Selected_Data %>% as_tibble() %>%  sample_frac(0.90)
HeldOutSamples <- Selected_Data %>% filter(shortnam %notin% HoldOutDataset$shortnam)

HoldOutDataset
HeldOutSamples

n <- nrow(HoldOutDataset)
#Target
Y <- HoldOutDataset %>% select(logpgp95) %>%  as.matrix
#Included Exogenous:
A_1 <- HoldOutDataset %>% select(Intercept) %>%  as.matrix
#All Exogenous
A <- HoldOutDataset %>% select(logem4,Intercept) %>% as.matrix
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

ols
iv
pulse

HeldOutSamples %>%  
  select(logpgp95,avexpr,Intercept) %>% 
  mutate(res = logpgp95- (avexpr*ols["avexpr",]+ols["Intercept",]) ) %>% 
  summarize(MSPE = sum(res^2))

HeldOutSamples %>%  
  select(logpgp95,avexpr,Intercept) %>% 
  mutate(res = logpgp95- (avexpr*iv["avexpr",]+iv["Intercept",]) ) %>% 
  summarize(MSPE = sum(res^2))

HeldOutSamples %>%  
  select(logpgp95,avexpr,Intercept) %>% 
  mutate(res = logpgp95- (avexpr*pulse["avexpr",]+pulse["Intercept",]) ) %>% 
  summarize(MSPE = sum(res^2))





#### HOLD OUT RANDOM SECOND..









holdoutrandom <- function(Data,percentage){
  HoldOutDataset <- Data %>% 
    as_tibble() %>%  
    sample_frac(percentage)
  HeldOutSamples <- Data %>% 
    filter(ID %notin% HoldOutDataset$ID)
  
  
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
  ols <- K_class(0,A,Z,Y,n)
  #IV
  iv <- K_class(1,A,Z,Y,n)
  #PULSE
  pulse <- PULSE(A,A_1,X,Y,p=0.05,N=10000,n)
  
  
  MSPE <- HeldOutSamples %>%  
    select(LWKLYWGE,AGEQ,AGEQSQ,RACE, MARRIED, SMSA, NEWENG, MIDATL, ENOCENT, WNOCENT, SOATL, ESOCENT, WSOCENT, MT,EDUC) %>% 
    mutate(resOLS = LWKLYWGE - 
             (  EDUC*ols["EDUC",]+
                  AGEQ*ols["AGEQ",]+
                  AGEQSQ*ols["AGEQSQ",]+
                  RACE*ols["RACE",]+
                  MARRIED*ols["MARRIED",]+
                  SMSA*ols["SMSA",]+
                  NEWENG*ols["NEWENG",]+
                  MIDATL*ols["MIDATL",]+
                  ENOCENT*ols["ENOCENT",]+
                  WNOCENT*ols["WNOCENT",]+
                  SOATL*ols["SOATL",]+
                  ESOCENT*ols["ESOCENT",]+
                  WSOCENT*ols["WSOCENT",]+
                  MT*ols["MT",]),
           resIV = LWKLYWGE - 
             (EDUC*iv["EDUC",]+
                AGEQ*iv["AGEQ",]+
                AGEQSQ*iv["AGEQSQ",]+
                RACE*iv["RACE",]+
                MARRIED*iv["MARRIED",]+
                SMSA*iv["SMSA",]+
                NEWENG*iv["NEWENG",]+
                MIDATL*iv["MIDATL",]+
                ENOCENT*iv["ENOCENT",]+
                WNOCENT*iv["WNOCENT",]+
                SOATL*iv["SOATL",]+
                ESOCENT*iv["ESOCENT",]+
                WSOCENT*iv["WSOCENT",]+
                MT*iv["MT",]),
           resPULSE = LWKLYWGE - 
             (EDUC*pulse["EDUC","LWKLYWGE"]+
                AGEQ*pulse["AGEQ","LWKLYWGE"]+
                AGEQSQ*pulse["AGEQSQ","LWKLYWGE"]+
                RACE*pulse["RACE","LWKLYWGE"]+
                MARRIED*pulse["MARRIED","LWKLYWGE"]+
                SMSA*pulse["SMSA","LWKLYWGE"]+
                NEWENG*pulse["NEWENG","LWKLYWGE"]+
                MIDATL*pulse["MIDATL","LWKLYWGE"]+
                ENOCENT*pulse["ENOCENT","LWKLYWGE"]+
                WNOCENT*pulse["WNOCENT","LWKLYWGE"]+
                SOATL*pulse["SOATL","LWKLYWGE"]+
                ESOCENT*pulse["ESOCENT","LWKLYWGE"]+
                WSOCENT*pulse["WSOCENT","LWKLYWGE"]+
                MT*pulse["MT","LWKLYWGE"])) %>% 
    summarize(OLS = round(sum(resOLS^2)/n(),3),
              IV = round(sum(resIV^2)/n(),3),
              PULSE = round(sum(resPULSE^2)/n(),3),
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
              PULSEmIV = PULSE -IV)
  return(MSPE)
}
holdoutrandom(Selected_Data %>% sample_frac(0.005),0.9)


summarydat <-data.frame(n=seq(1,100,1)) %>% 
  rowwise() %>% 
  mutate(order=list(holdoutrandom(Selected_Data %>% sample_frac(0.1),0.9)))

summarydat %>%
  unnest(cols=c(order)) %>% 
  group_by(ORDER) %>% 
  summarise(count= n(),
            OLS = mean(OLS),
            PULSE = mean(PULSE),
            IV = mean(IV)) %>% 
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
                                ORDER == "PULSE<IV<OLS" ~ OLS-IV))





holdoutextremes <- function(Data,noExtremes){
  HoldOutDataset <- Data %>% 
    as_tibble() %>%  
    arrange(AGEQ) %>% 
    head(-noExtremes/2) %>% 
    tail(-noExtremes/2)
  
  HeldOutSamples <- Data %>% 
    filter(ID %notin% HoldOutDataset$ID)
  
  
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
  ols <- K_class(0,A,Z,Y,n)
  #IV
  iv <- K_class(1,A,Z,Y,n)
  #PULSE
  pulse <- PULSE(A,A_1,X,Y,p=0.05,N=10000,n)
  
  
  MSPE <- HeldOutSamples %>%  
    select(LWKLYWGE,AGEQ,AGEQSQ,RACE, MARRIED, SMSA, NEWENG, MIDATL, ENOCENT, WNOCENT, SOATL, ESOCENT, WSOCENT, MT,EDUC) %>% 
    mutate(resOLS = LWKLYWGE - 
             (  EDUC*ols["EDUC",]+
                  AGEQ*ols["AGEQ",]+
                  AGEQSQ*ols["AGEQSQ",]+
                  RACE*ols["RACE",]+
                  MARRIED*ols["MARRIED",]+
                  SMSA*ols["SMSA",]+
                  NEWENG*ols["NEWENG",]+
                  MIDATL*ols["MIDATL",]+
                  ENOCENT*ols["ENOCENT",]+
                  WNOCENT*ols["WNOCENT",]+
                  SOATL*ols["SOATL",]+
                  ESOCENT*ols["ESOCENT",]+
                  WSOCENT*ols["WSOCENT",]+
                  MT*ols["MT",]),
           resIV = LWKLYWGE - 
             (EDUC*iv["EDUC",]+
                AGEQ*iv["AGEQ",]+
                AGEQSQ*iv["AGEQSQ",]+
                RACE*iv["RACE",]+
                MARRIED*iv["MARRIED",]+
                SMSA*iv["SMSA",]+
                NEWENG*iv["NEWENG",]+
                MIDATL*iv["MIDATL",]+
                ENOCENT*iv["ENOCENT",]+
                WNOCENT*iv["WNOCENT",]+
                SOATL*iv["SOATL",]+
                ESOCENT*iv["ESOCENT",]+
                WSOCENT*iv["WSOCENT",]+
                MT*iv["MT",]),
           resPULSE = LWKLYWGE - 
             (EDUC*pulse["EDUC","LWKLYWGE"]+
                AGEQ*pulse["AGEQ","LWKLYWGE"]+
                AGEQSQ*pulse["AGEQSQ","LWKLYWGE"]+
                RACE*pulse["RACE","LWKLYWGE"]+
                MARRIED*pulse["MARRIED","LWKLYWGE"]+
                SMSA*pulse["SMSA","LWKLYWGE"]+
                NEWENG*pulse["NEWENG","LWKLYWGE"]+
                MIDATL*pulse["MIDATL","LWKLYWGE"]+
                ENOCENT*pulse["ENOCENT","LWKLYWGE"]+
                WNOCENT*pulse["WNOCENT","LWKLYWGE"]+
                SOATL*pulse["SOATL","LWKLYWGE"]+
                ESOCENT*pulse["ESOCENT","LWKLYWGE"]+
                WSOCENT*pulse["WSOCENT","LWKLYWGE"]+
                MT*pulse["MT","LWKLYWGE"])) %>% 
    summarize(OLS = round(sum(resOLS^2)/n(),3),
              IV = round(sum(resIV^2)/n(),3),
              PULSE = round(sum(resPULSE^2)/n(),3),
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
              PULSEmIV = PULSE -IV)
  
  return(MSPE)
}

summarydat <-data.frame(n=seq(50000,250000 ,50000)) %>% 
  rowwise() %>% 
  mutate(order=list(holdoutextremes(Selected_Data_Centered,n)))%>%
  unnest(cols=c(order))



