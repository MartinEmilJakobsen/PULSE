library(tidyverse)
library(dplyr)
library(magrittr)
library(stringr)
library(grid)
library(gridExtra)


########################################
####### Multivariate Experiments ####### 
#######          Plots           ####### 
########################################

#############################
#### Varying confounding ####
#############################

#########################
#### Beta 00 PULSE05 ####
#########################

# Read data
Data_Location <- "Data/Experiment_Multivariate_VaryingConfounding_Beta00_PULSE05_nSim_5000_nObsPerSim_50_nModel_10000_20200716080129.RDS"
ID <- "20200716080129"
dat <- readRDS(file=Data_Location) 


# Finding models where PULSE is MSE superior to Fuller(4)
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

# Computing relative change in performance measures
LossData <- dat  %>%
  select(nModel,nSim,n,Type,MeanGn,Determinant,Trace,Bias) %>% 
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
  gather(Type,Value,c(-nModel,-nSim,-n,-MeanGn,-pm,-Ful1,-Ful4,-OLS,-PULSE05,-MinEigenMeanGn,-MaxEigenMeanGn,-Superior))


# Merging Data
PlotData <- left_join(left_join(LossData,Optimal,by=c("n","nModel")) ,Coefs,by=c("nModel"))


# Plot for relative change in performance measures
ggplot(data=PlotData) +
  geom_hline(yintercept =0,color="black",linetype="solid") +
  geom_vline(xintercept =log(15.5),color="black",linetype="dotted") +
  geom_point(aes(x=log(MinEigenMeanGn),y=Value,color=sqrt(RhoSq)),alpha=0.5,size=0.8)+
  scale_color_gradient(low = "blue", high = "orange")+
  facet_wrap(    Type ~ pm ,scales="free_y",ncol=3, labeller = label_both)+
  xlab(expression(log(lambda[min](hat(E)[N](G[n])))))+
  labs(colour=expression("||"*rho*"||"[2]))+
  ylab(expression(paste("Relative change of performance measure")))+
  theme(plot.margin = unit(c(0,0.8,0,0), "cm"))+ 
  theme(legend.position="bottom")

ggsave(paste0("Plots/Multivariate_VaryingConfounding_Beta00_PULSE05_",ID,".png"), plot = last_plot(), device = NULL, path = NULL,
       scale = 1, width = 12, height = 9, units = c("in"),
       dpi = 200, limitsize = FALSE)


########################
### pm Superior plot ###
########################


ggplot(data=PlotData) +
  geom_hline(yintercept =0,color="black",linetype="solid") +
  geom_vline(xintercept =log(15.5),color="black",linetype="dotted") +
  geom_point(aes(x=log(MinEigenMeanGn),y=Value,color=Superior),alpha=0.5,size=0.8)+
  facet_wrap(    Type ~ pm ,scales="free_y",ncol=3, labeller = label_both)+
  xlab(expression(log(lambda[min](hat(E)[N](G[n])))))+
  labs(colour="Performance Measure Superiorty - Fuller(4) vs PULSE(05)")+
  ylab(expression(paste("Relative change of performance measure")))+
  theme(plot.margin = unit(c(0,0.8,0,0), "cm"))+
  theme(legend.position="bottom")

ggsave(paste0("Plots/Multivariate_VaryingConfounding_Beta00_PULSE05_pmSuperior_",ID,".png"), plot = last_plot(), device = NULL, path = NULL,
       scale = 1, width = 12, height = 9, units = c("in"),
       dpi = 200, limitsize = FALSE)


##########################
### MSE Superior Plots ###
##########################





p1 <- ggplot(data=PlotData %>% arrange(TrueSuperior)) +
  geom_hline(yintercept =0,color="black",linetype="solid") +
  geom_vline(xintercept =log(15.5),color="black",linetype="dotted") +
  geom_point(aes(x=log(MinEigenMeanGn),y=Value,color=TrueSuperior),alpha=0.6,size=0.8)+
  facet_wrap(    Type ~ pm ,scales="free_y",ncol=3, labeller = label_both)+
  xlab(expression(log(lambda[min](hat(E)[N](G[n])))))+
  labs(colour="MSE Superiority - Fuller(4) vs PULSE(05)")+
  ylab(expression(paste("Relative change of performance measure")))+
  theme(plot.margin = unit(c(0,0.8,0,0), "cm"))+
  theme(legend.position="bottom")+ 
  theme(axis.title.y = element_blank())+ 
  theme(axis.title.x = element_blank())+
  scale_color_manual(values=c("#999999", "#FF0000"))


# Read Superior Models Rerun data

Data_Location <- "Data/Experiment_Multivariate_VaryingConfounding_SuperiorModels_nSim_25000_nObsPerSim_50_20200716154413.RDS"
ID2 <- "20200716154413"
dat <- readRDS(file=Data_Location) 


# Finding optimal
Optimal <- dat %>% 
  select(n,nModel,Type,MSE) %>%
  unique() %>% 
  spread(Type,MSE) %>%  
  rowwise() %>% 
  mutate(TrueSuperior = case_when( min(eigen(PULSE05-Ful4)$values) > 0 ~ "Ful4",
                                   min(eigen(Ful4-PULSE05)$values) > 0 ~ "PULSE05",
                                   TRUE ~ "Not comparable")) %>% 
  select(n,nModel,TrueSuperior) 

Optimal %>% group_by(TrueSuperior) %>%  summarise(count = n())

CoefsSup <- dat %>% 
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
    "rhosq"),461)
    
  ) %>% 
  unnest(cols=c(ModelCoefs)) %>% 
  spread(Coef,ModelCoefs) %>% 
  rename(RhoSq=rhosq)

LossDataSup <- dat  %>%
  select(nModel,nSim,n,Type,MeanGn,Determinant,Trace,Bias) %>% 
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
  gather(Type,Value,c(-nModel,-nSim,-n,-MeanGn,-pm,-Ful1,-Ful4,-OLS,-PULSE05,-MinEigenMeanGn,-MaxEigenMeanGn,-Superior))


PlotDataSup <- left_join(left_join(LossDataSup,Optimal,by=c("n","nModel")) ,CoefsSup,by=c("nModel"))


p2 <- ggplot(data=PlotDataSup %>% filter(nModel %in% Optimal$nModel)) +
  geom_hline(yintercept =0,color="black",linetype="solid") +
  geom_vline(xintercept =log(15.5),color="black",linetype="dotted") +
  geom_point(aes(x=log(MinEigenMeanGn),y=Value,color=sqrt(RhoSq)),alpha=0.6,size=0.8)+
  scale_color_gradient(low = "blue", high = "orange")+
  facet_wrap(    Type ~ pm ,scales="free_y",ncol=3, labeller = label_both)+
  xlab(expression(log(lambda[min](hat(E)[N](G[n])))))+
  labs(colour=expression("||"*rho*"||"[2]))+
  ylab(expression(paste("Relative change of performance measure")))+
  theme(plot.margin = unit(c(0,0.8,0,0), "cm"))+
  theme(legend.position="bottom")+ 
  theme(axis.title.y = element_blank())



plot <- arrangeGrob(p1,p2, ncol=1,widths=c(1),left = textGrob("Relative Change in Performance Measure", rot = 90, vjust = 1))


ggsave(paste0("Plots/Multivariate_VaryingConfounding_Beta00_PULSE05_MSE_Superior",ID,"_",ID2,".png"), plot = plot, device = NULL, path = NULL,
       scale = 1, width = 12, height = 18, units = c("in"),
       dpi = 200, limitsize = FALSE)


#########################
#### Beta 00 PULSE10 ####
#########################

Data_Location <- "Data/Experiment_Multivariate_VaryingConfounding_Beta00_PULSE10_nSim_5000_nObsPerSim_50_nModel_10000_20200716122010.RDS"
ID <- "20200716122010"
dat <- readRDS(file=Data_Location)


Optimal <- dat %>% 
  select(n,nModel,Type,MSE) %>%
  unique() %>% 
  spread(Type,MSE) %>%  
  rowwise() %>% 
  mutate(TrueSuperior = case_when( min(eigen(PULSE10-Ful4)$values) > 0 ~ "Ful4",
                                   min(eigen(Ful4-PULSE10)$values) > 0 ~ "PULSE10",
                                   TRUE ~ "Not comparable")) %>% 
  select(n,nModel,TrueSuperior) 

Rhosq <- dat %>% 
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
  filter(Coef=="rhosq") %>% 
  unnest(cols=c(ModelCoefs))
names(Rhosq) <- c("nModel","RhoSq","t")
Rhosq <- Rhosq %>%  mutate(RhoSq= as.numeric(RhoSq)) %>% select(-t)

LossData <- dat  %>%
  select(nModel,nSim,n,Type,MeanGn,Determinant,Trace,Bias) %>% 
  gather(pm, Value, c("Determinant",
                      "Trace",
                      "Bias")) %>% 
  mutate(pm = factor(pm, levels = c("Determinant",
                                    "Trace",
                                    "Bias"))) %>% 
  spread(Type,Value) %>%
  ungroup() %>% 
  mutate("PULSE10 to Fuller4" = pmap_dbl(.l=list(Ful4,PULSE10,pm), .f=function(Ful4,PULSE10,pm){ (Ful4-PULSE10)/PULSE10 }),
         "PULSE10 to Fuller1"  = pmap_dbl(.l=list(Ful1,PULSE10,pm), .f=function(Ful1,PULSE10,pm){ (Ful1-PULSE10)/PULSE10}),
         "PULSE10 to OLS"  = pmap_dbl(.l=list(OLS,PULSE10,pm), .f=function(OLS,PULSE10,pm){ (OLS-PULSE10)/PULSE10}),
         MinEigenMeanGn = map_dbl(.x= MeanGn,.f= function(x) {ifelse(is.double(eigen(x)$values),min(eigen(x)$values),NA)}),
         MaxEigenMeanGn = map_dbl(.x= MeanGn,.f= function(x) {ifelse(is.double(eigen(x)$values),max(eigen(x)$values),NA)})) %>%
  mutate(Superior = ifelse(Ful4<PULSE10,"Ful4","PULSE")) %>% 
  gather(Type,Value,c(-nModel,-nSim,-n,-MeanGn,-pm,-Ful1,-Ful4,-OLS,-PULSE10,-MinEigenMeanGn,-MaxEigenMeanGn,-Superior))

PlotData <- left_join(left_join(LossData,Optimal,by=c("n","nModel")) ,Rhosq,by=c("nModel"))

ggplot(data=PlotData) +
  geom_hline(yintercept =0,color="black",linetype="solid") +
  geom_vline(xintercept =log(15.5),color="black",linetype="dotted") +
  geom_point(aes(x=log(MinEigenMeanGn),y=Value,color=sqrt(RhoSq)),alpha=0.5,size=0.8)+
  scale_color_gradient(low = "blue", high = "orange")+
  facet_wrap(    Type ~ pm ,scales="free_y",ncol=3, labeller = label_both)+
  xlab(expression(log(lambda[min](hat(E)[N](G[n])))))+
  labs(colour=expression("||"*rho*"||"[2]))+
  ylab(expression(paste("Relative change of performance measure")))+
  theme(plot.margin = unit(c(0,0.8,0,0), "cm"))+ 
  theme(legend.position="bottom")

ggsave(paste0("Plots/Multivariate_VaryingConfounding_Beta00_PULSE10_",ID,".png"), plot = last_plot(), device = NULL, path = NULL,
       scale = 1, width = 12, height = 9, units = c("in"),
       dpi = 200, limitsize = FALSE)



#########################
#### Beta 11 PULSE05 ####
#########################

# Read data
Data_Location <- "Data/Experiment_Multivariate_VaryingConfounding_Beta11_PULSE05_nSim_5000_nObsPerSim_50_nModel_10000_20200716180808.RDS"
ID <- "20200716180808"
dat <- readRDS(file=Data_Location)

Optimal <- dat %>% 
  select(n,nModel,Type,MSE) %>%
  unique() %>% 
  spread(Type,MSE) %>%  
  rowwise() %>% 
  mutate(TrueSuperior = case_when( min(eigen(PULSE05-Ful4)$values) > 0 ~ "Ful4",
                                   min(eigen(Ful4-PULSE05)$values) > 0 ~ "PULSE05",
                                   TRUE ~ "Not comparable")) %>% 
  select(n,nModel,TrueSuperior) 

Rhosq <- dat %>% 
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
  filter(Coef=="rhosq") %>% 
  unnest(cols=c(ModelCoefs))
names(Rhosq) <- c("nModel","RhoSq","t")
Rhosq <- Rhosq %>%  mutate(RhoSq= as.numeric(RhoSq)) %>% select(-t)

LossData <- dat  %>%
  select(nModel,nSim,n,Type,MeanGn,Determinant,Trace,Bias) %>% 
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
  gather(Type,Value,c(-nModel,-nSim,-n,-MeanGn,-pm,-Ful1,-Ful4,-OLS,-PULSE05,-MinEigenMeanGn,-MaxEigenMeanGn,-Superior))

Optimal %>% group_by(TrueSuperior) %>%  summarise(count = n())

PlotData <- left_join(left_join(LossData,Optimal,by=c("n","nModel")) ,Rhosq,by=c("nModel"))

ggplot(data=PlotData) +
  geom_hline(yintercept =0,color="black",linetype="solid") +
  geom_vline(xintercept =log(15.5),color="black",linetype="dotted") +
  geom_point(aes(x=log(MinEigenMeanGn),y=Value,color=sqrt(RhoSq)),alpha=0.5,size=0.8)+
  scale_color_gradient(low = "blue", high = "orange")+
  facet_wrap(    Type ~ pm ,scales="free_y",ncol=3, labeller = label_both)+
  xlab(expression(log(lambda[min](hat(E)[N](G[n])))))+
  labs(colour=expression("||"*rho*"||"[2]))+
  ylab(expression(paste("Relative change of performance measure")))+
  theme(plot.margin = unit(c(0,0.8,0,0), "cm"))+ 
  theme(legend.position="bottom")

ggsave(paste0("Plots/Multivariate_VaryingConfounding_Beta11_PULSE05_",ID,".png"), plot = last_plot(), device = NULL, path = NULL,
       scale = 1, width = 12, height = 9, units = c("in"),
       dpi = 200, limitsize = FALSE)

##########################
#### Beta -11 PULSE05 ####
##########################

# Read data
Data_Location <- "Data/Experiment_Multivariate_VaryingConfounding_Beta-11_PULSE05_nSim_5000_nObsPerSim_50_nModel_10000_20200716033957.RDS"
ID <- "20200716033957"
dat <- readRDS(file=Data_Location)


# Find Optimal Models
Optimal <- dat %>% 
  select(n,nModel,Type,MSE) %>%
  unique() %>% 
  spread(Type,MSE) %>%  
  rowwise() %>% 
  mutate(TrueSuperior = case_when( min(eigen(PULSE05-Ful4)$values) > 0 ~ "Ful4",
                                   min(eigen(Ful4-PULSE05)$values) > 0 ~ "PULSE05",
                                   TRUE ~ "Not comparable")) %>% 
  select(n,nModel,TrueSuperior) 


# Extract coefficients
Rhosq <- dat %>% 
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
  filter(Coef=="rhosq") %>% 
  unnest(cols=c(ModelCoefs))
names(Rhosq) <- c("nModel","RhoSq","t")
Rhosq <- Rhosq %>%  mutate(RhoSq= as.numeric(RhoSq)) %>% select(-t)

LossData <- dat  %>%
  select(nModel,nSim,n,Type,MeanGn,Determinant,Trace,Bias) %>% 
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
  gather(Type,Value,c(-nModel,-nSim,-n,-MeanGn,-pm,-Ful1,-Ful4,-OLS,-PULSE05,-MinEigenMeanGn,-MaxEigenMeanGn,-Superior))

PlotData <- left_join(left_join(LossData,Optimal,by=c("n","nModel")) ,Rhosq,by=c("nModel"))

ggplot(data=PlotData) +
  geom_hline(yintercept =0,color="black",linetype="solid") +
  geom_vline(xintercept =log(15.5),color="black",linetype="dotted") +
  geom_point(aes(x=log(MinEigenMeanGn),y=Value,color=sqrt(RhoSq)),alpha=0.5,size=0.8)+
  scale_color_gradient(low = "blue", high = "orange")+
  facet_wrap(    Type ~ pm ,scales="free_y",ncol=3, labeller = label_both)+
  xlab(expression(log(lambda[min](hat(E)[N](G[n])))))+
  labs(colour=expression("||"*rho*"||"[2]))+
  ylab(expression(paste("Relative change of performance measure")))+
  theme(plot.margin = unit(c(0,0.8,0,0), "cm"))+ 
  theme(legend.position="bottom")


ggsave(paste0("Plots/Multivariate_VaryingConfounding_Beta-11_PULSE05_",ID,".png"), plot = last_plot(), device = NULL, path = NULL,
       scale = 1, width = 12, height = 9, units = c("in"),
       dpi = 200, limitsize = FALSE)

###########################
#### Fixed confounding ####
###########################

# Read data normRho = 0.2,0.8, eta =0.2,0.8
Data_Location <- "Data/Experiment_Multivariate_FixedConfounding_nSim_5000_nObsPerSim_50_nModel_5000_20200716060820.RDS"
ID <- "20200716060820"
dat <- readRDS(file=Data_Location)

#Finding MSE superior models
Optimal <- dat %>% 
  select(n,nModel,Type,MSE,Cov) %>%
  unique() %>% 
  spread(Type,MSE) %>%  
  rowwise() %>% 
  mutate(TrueSuperior = case_when( min(eigen(PULSE05-Ful4)$values) > 0 ~ "Ful4",
                                   min(eigen(Ful4-PULSE05)$values) > 0 ~ "PULSE",
                                   TRUE ~ "Not comparable")) %>% 
  select(n,nModel,Cov,TrueSuperior) 

#Extracting models coefficients
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

# Calculating relative change in performance measures
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


# Gathering Data for plot
PlotData <- left_join(left_join(LossData,Optimal,by=c("n","nModel","Cov")) ,Cors,by=c("nModel","Cov")) %>%  
  arrange(normRho,-eta) %>% 
  mutate(Label = paste0("'||'*rho*'||'[2]:'",sprintf("%.3f",normRho),"'~~eta:'",sprintf("%.3f",eta),"'~~phi[1]:'",sprintf("%.3f",phi1),"'~~phi[2]:'",sprintf("%.3f",phi2),"'")) 

scaleFUN <- function(x) sprintf("%.2f", x)

p1 <- ggplot(data=PlotData %>% filter(Type=="PULSE05 to Fuller4",Cov %in% c(1))) +
  geom_hline(yintercept =0,color="black",linetype="solid") +
  geom_vline(xintercept =log(15.5),color="black",linetype="dotted") +
  geom_point(aes(x=log(MinEigenMeanGn),y=Value),alpha=0.1,size=1)+
  facet_wrap( Label ~ pm ,scales="free_y",ncol=3, labeller = label_parsed)+
  xlab(expression(log(lambda[min](hat(E)[N](G[n])))))+
  ylab(expression(paste("Relative change of performance measure")))+
  theme(plot.margin = unit(c(0,0.8,0,0), "cm"))+ 
  theme(axis.title.y = element_blank())+
  theme(axis.title.x = element_blank())+ scale_y_continuous(labels=scaleFUN)


p2 <- ggplot(data=PlotData %>% filter(Type=="PULSE05 to Fuller4",Cov %in% c(2))) +
  geom_hline(yintercept =0,color="black",linetype="solid") +
  geom_vline(xintercept =log(15.5),color="black",linetype="dotted") +
  geom_point(aes(x=log(MinEigenMeanGn),y=Value),alpha=0.1,size=1)+
  facet_wrap(  Label ~ pm ,scales="free_y",ncol=3, labeller = label_parsed)+
  xlab(expression(log(lambda[min](hat(E)[N](G[n])))))+
  ylab(expression(paste("Relative change of performance measure")))+
  theme(plot.margin = unit(c(0,0.8,0,0), "cm"))+ 
  theme(axis.title.y = element_blank())+ 
  theme(axis.title.x = element_blank())+ scale_y_continuous(labels=scaleFUN)


p3 <- ggplot(data=PlotData %>% filter(Type=="PULSE05 to Fuller4",Cov %in% c(3))) +
  geom_hline(yintercept =0,color="black",linetype="solid") +
  geom_vline(xintercept =log(15.5),color="black",linetype="dotted") +
  geom_point(aes(x=log(MinEigenMeanGn),y=Value),alpha=0.1,size=1)+
  facet_wrap(  Label ~ pm ,scales="free_y",ncol=3, labeller = label_parsed)+
  xlab(expression(log(lambda[min](hat(E)[N](G[n])))))+
  ylab(expression(paste("Relative change of performance measure")))+
  theme(plot.margin = unit(c(0,0.8,0,0), "cm"))+ 
  theme(axis.title.y = element_blank())+ 
  theme(axis.title.x = element_blank())+ scale_y_continuous(labels=scaleFUN)

p4 <- ggplot(data=PlotData %>% filter(Type=="PULSE05 to Fuller4",Cov %in% c(4))) +
  geom_hline(yintercept =0,color="black",linetype="solid") +
  geom_vline(xintercept =log(15.5),color="black",linetype="dotted") +
  geom_point(aes(x=log(MinEigenMeanGn),y=Value),alpha=0.1,size=1)+
  facet_wrap(  Label ~ pm ,scales="free_y",ncol=3, labeller = label_parsed)+
  xlab(expression(log(lambda[min](hat(E)[N](G[n])))))+
  ylab(expression(paste("Relative change of performance measure")))+
  theme(plot.margin = unit(c(0,0.8,0,0), "cm"))+ theme(axis.title.y = element_blank())+ scale_y_continuous(labels=scaleFUN)


plot <- arrangeGrob(p1,p2,p3,p4, ncol=1,left = textGrob("Relative Change in Performance Measure", rot = 90, vjust = 0))

ggsave(paste0("Plots/Multivariate_FixedConfounding_",ID,".png"), plot =plot, device = NULL, path = NULL,
       scale = 1, width = 12, height = 10, units = c("in"),
       dpi = 200, limitsize = FALSE)
