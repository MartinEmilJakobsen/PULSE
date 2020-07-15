library(tidyverse)
library(magrittr)
library(stringr)
library(gridExtra)


#############################################
######### Multivariate Experiments   ########
######### Fixed Confounding Analysis ########
#############################################

# Read data normRho = 0.447,0.949, eta =0.2,0.8
# Data_Location <- "Data/Experiment_Multivariate_FixedConfounding_nSim_5000_nObsPerSim_50_nModel_5000_20200714132856.RDS"
# ID <- "20200714132856"
# dat <- readRDS(file=Data_Location)

# Read data normRho = 0.2,0.8, eta =0.2,0.8
Data_Location <- "Data/Experiment_Multivariate_FixedConfounding_nSim_5000_nObsPerSim_50_nModel_5000_20200715005001.RDS"
ID <- "20200715005001"
dat <- readRDS(file=Data_Location)


Optimal <- dat %>% 
  select(n,nModel,Type,MSE,Cov) %>%
  unique() %>% 
  spread(Type,MSE) %>%  
  rowwise() %>% 
  mutate(TrueSuperior = case_when( min(eigen(PULSE05-Ful4)$values) > 0 ~ "Ful4",
                                   min(eigen(Ful4-PULSE05)$values) > 0 ~ "PULSE",
                                   TRUE ~ "Not comparable")) %>% 
  select(n,nModel,Cov,TrueSuperior) 

Cors <- dat %>% 
  filter(Type=="Ful1") %>%  
  select(nModel,ModelCoefs,Cov) %>% 
  unnest(cols=c(ModelCoefs)) %>% 
  mutate(Coef = rep(c(
    "xi11",
    "xi12",
    "xi21",
    "xi22",
    "phi1",
    "phi2",
    "eta"),20000)  ) %>% 
  unnest(cols=c(ModelCoefs)) %>% 
  spread(Coef,ModelCoefs) %>% 
  mutate(normRho = sqrt((phi1^2+phi2^2-2*eta*phi1*phi2)/(1-eta^2)))


LossData <- dat  %>%
  select(nModel,Cov,nSim,n,Type,MeanGn,Determinant,Trace,Bias) %>% 
  gather(pm, Value, c("Determinant",
                      "Trace",
                      "Bias")) %>% 
  mutate(pm = factor(pm, levels = c("Determinant",
                                    "Trace",
                                    "Bias"))) %>% 
  spread(Type,Value) %>%
  ungroup() %>% 
  mutate("PULSE05 to Fuller4" = pmap_dbl(.l=list(Ful4,PULSE05,pm), .f=function(Ful4,PULSE05,pm){ (Ful4-PULSE05)/PULSE05 }),
         "PULSE05 to Fuller1"  = pmap_dbl(.l=list(Ful1,PULSE05,pm), .f=function(Ful1,PULSE05,pm){ (Ful1-PULSE05)/PULSE05}),
         "PULSE05 to OLS"  = pmap_dbl(.l=list(OLS,PULSE05,pm), .f=function(OLS,PULSE05,pm){ (OLS-PULSE05)/PULSE05}),
         MinEigenMeanGn = map_dbl(.x= MeanGn,.f= function(x) {ifelse(is.double(eigen(x)$values),min(eigen(x)$values),NA)}),
         MaxEigenMeanGn = map_dbl(.x= MeanGn,.f= function(x) {ifelse(is.double(eigen(x)$values),max(eigen(x)$values),NA)})) %>%
  mutate(Superior = ifelse(Ful4<PULSE05,"Ful4","PULSE")) %>% 
  gather(Type,Value,c(-nModel,-Cov,-nSim,-n,-MeanGn,-pm,-Ful1,-Ful4,-OLS,-PULSE05,-MinEigenMeanGn,-MaxEigenMeanGn,-Superior))


Optimal %>% group_by(TrueSuperior) %>%  summarise(count = n())


# Percentage Better
left_join(
  LossData %>% filter(Ful4>= PULSE05) %>% select(Cov,pm,nModel,Type)  %>% group_by(Cov,pm,Type) %>% summarise(Better=n()) %>% select(-Type) %>% unique(),
  LossData %>% filter(PULSE05>= Ful4) %>% select(Cov,pm,nModel,Type)  %>% group_by(Cov,pm,Type) %>% summarise(Worse=n())%>% select(-Type) %>% unique(),
  by = c("Cov","pm")) %>% 
  mutate(sum = Better+Worse, fractionBetter = round(100*Better/5000,2)) %>% arrange(pm) %>% 
  left_join(Cors %>%  filter(nModel==1),by="Cov")
