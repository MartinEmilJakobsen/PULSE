library(tidyverse)
library(magrittr)
library(stringr)
library(gridExtra)


########################################
####### Multivariate Experiments ####### 
#######  Selecting MSE superior  ####### 
########################################

#########################
#### Beta 00 PULSE05 ####
#########################

# Read data
Data_Location <- "Data/Experiment_Multivariate_VaryingConfounding_Beta00_PULSE05_nSim_5000_nObsPerSim_50_nModel_10000_20200714024302.RDS"
ID <- "20200714024302"
dat <- readRDS(file=Data_Location) 

# Calculating MSE superiority percentages
Optimal <- dat %>% 
  select(n,nModel,Type,MSE) %>%
  unique() %>% 
  spread(Type,MSE) %>%  
  rowwise() %>% 
  mutate(TrueSuperior = case_when( min(eigen(PULSE05-Ful4)$values) > 0 ~ "Ful4",
                                   min(eigen(Ful4-PULSE05)$values) > 0 ~ "PULSE05",
                                   TRUE ~ "Not comparable")) %>% 
  select(n,nModel,TrueSuperior) 


# Extracting Model Coefficients from data
Coefs <- dat %>% 
  filter(Type=="Ful1") %>%  
  select(nModel,ModelCoefs) %>% 
  unnest(cols=c(ModelCoefs)) %>% 
  mutate(Coef = rep(c(
    "xi11",
    "xi12",
    "xi21",
    "xi22",
    "delta11",
    "delta12",
    "delta21",
    "delta22",
    "mu11",
    "mu22",
    "VepX1" ,
    "VepX2" ,
    "rhosq"),10000)
    
  ) %>% 
  unnest(cols=c(ModelCoefs)) %>% 
  spread(Coef,ModelCoefs) %>% 
  rename(RhoSq=rhosq)

# Saving Models where PULSE is MSE superior to Fuller(4) for rerun to account for selection bias.

ModelDataForSave <- Coefs %>% filter(nModel %in% {Optimal %>% filter(TrueSuperior=="PULSE05") %$% nModel} )

saveRDS(ModelDataForSave,file=paste0("Data/VaryingConfounding_MSESuperiorModelData_",ID,".RDS"))

