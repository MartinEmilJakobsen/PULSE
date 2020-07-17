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


### Analysis of Initial Experiment:

# Read data
Data_Location <- "Data/Experiment_Multivariate_VaryingConfounding_Beta00_PULSE05_nSim_5000_nObsPerSim_50_nModel_10000_20200716080129.RDS"
ID <- "20200716080129"
dat <- readRDS(file=Data_Location) 

# Calculating MSE superiority percentages
Optimal <- dat %>% 
  select(n,nModel,Type,MSE) %>%
  unique() %>% 
  spread(Type,MSE) %>%  
  rowwise() %>% 
  mutate(TrueSuperior = case_when( min(eigen(PULSE05-Ful4)$values) > 0 ~ "Ful4",
                                   min(eigen(Ful4-PULSE05)$values) > 0 ~ "PULSE05",
                                   TRUE ~ "Not comparable"),
         TrueSuperiorVsFul1 = case_when( min(eigen(PULSE05-Ful1)$values) > 0 ~ "Ful1",
                                         min(eigen(Ful1-PULSE05)$values) > 0 ~ "PULSE05",
                                         TRUE ~ "Not comparable"),
         TrueSuperiorVsOLS = case_when( min(eigen(PULSE05-OLS)$values) > 0 ~ "OLS",
                                        min(eigen(OLS-PULSE05)$values) > 0 ~ "PULSE05",
                                        TRUE ~ "Not comparable")) %>% 
  select(n,nModel,TrueSuperior,TrueSuperiorVsFul1,TrueSuperiorVsOLS) 


Ful4 <- Optimal %>% 
  group_by(TrueSuperior) %>%  
  summarise(count = n()) %>% 
  spread(TrueSuperior,count)
Ful1 <- Optimal %>% 
  group_by(TrueSuperiorVsFul1) %>%  
  summarise(count = n()) %>%  
  spread(TrueSuperiorVsFul1,count)
OLS <- Optimal %>% 
  group_by(TrueSuperiorVsOLS) %>%  
  summarise(count = n()) %>% 
  spread(TrueSuperiorVsOLS,count)

bind_rows(Ful4,Ful1,OLS) %>% 
  mutate_all(.funs=function(x){round(100*x/10000,1)}) %>% 
 # mutate_all(~replace(., is.na(.), 0.0)) %>% 
  mutate(Against = c("Fuller4","Fuller1","OLS"))


### Analysis of Experiment rerun on the optimal models from above:


Data_Location <- "Data/Experiment_Multivariate_VaryingConfounding_SuperiorModels_nSim_25000_nObsPerSim_50_20200716154413.RDS"
ID <- "20200716154413"
dat2 <- readRDS(file=Data_Location) 

Optimal <- dat2 %>% 
  select(n,nModel,Type,MSE) %>%
  unique() %>% 
  spread(Type,MSE) %>%  
  rowwise() %>% 
  mutate(TrueSuperior = case_when( min(eigen(PULSE05-Ful4)$values) > 0 ~ "Ful4",
                                   min(eigen(Ful4-PULSE05)$values) > 0 ~ "PULSE05",
                                   TRUE ~ "Not comparable"),
         TrueSuperiorVsFul1 = case_when( min(eigen(PULSE05-Ful1)$values) > 0 ~ "Ful1",
                                         min(eigen(Ful1-PULSE05)$values) > 0 ~ "PULSE05",
                                         TRUE ~ "Not comparable"),
         TrueSuperiorVsOLS = case_when( min(eigen(PULSE05-OLS)$values) > 0 ~ "OLS",
                                        min(eigen(OLS-PULSE05)$values) > 0 ~ "PULSE05",
                                        TRUE ~ "Not comparable")) %>% 
  select(n,nModel,TrueSuperior,TrueSuperiorVsFul1,TrueSuperiorVsOLS) 

Optimal %>% 
  group_by(TrueSuperior) %>%  
  summarise(count = n())


########################################
########## Count Error Models ########## 
########################################

Data_Location <- "Data/Experiment_Multivariate_VaryingConfounding_CountErrorModels_nSim_5000_nObsPerSim_50_nModel_10000_20200715041723.RDS"
ID <- "20200715041723"
dat <- readRDS(file=Data_Location) 

filter(dat,ErrorCount >0)

#No models encountered error in trycatch.