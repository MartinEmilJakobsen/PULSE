library(tidyverse)
library(magrittr)
library(stringr)
library(gridExtra)


########################################
####### Multivariate Experiments ####### 
#######          Plots           ####### 
########################################

#########################
#### Beta 00 PULSE10 ####
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
  mutate_all(~replace(., is.na(.), 0.0))


# Saving only the MSE superiority versus Fuller(4)
Optimal <- Optimal %>% 
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

ModelDataForSave <- Coefs %>% filter(nModel %in% {Optimal %>% filter(TrueSuperior=="PULSE05") %$% nModel} )

saveRDS(ModelDataForSave,file="AllRandomCoefs_SuperiorModelData.RDS")

LossData <- dat  %>%
  select(nModel,nSim,n,Type,MeanGn,Determinant,Trace,BiasTwoNorm) %>% 
  gather(pm, Value, c("Determinant",
                      "Trace",
                      "BiasTwoNorm")) %>% 
  mutate(pm = factor(pm, levels = c("Determinant",
                                    "Trace",
                                    "BiasTwoNorm"))) %>% 
  spread(Type,Value) %>%
  ungroup() %>% 
  mutate("PULSE05 to Fuller4" = pmap_dbl(.l=list(Ful4,PULSE05,pm), .f=function(Ful4,PULSE05,pm){ (Ful4-PULSE05)/PULSE05 }),
         "PULSE05 to Fuller1"  = pmap_dbl(.l=list(Ful1,PULSE05,pm), .f=function(Ful1,PULSE05,pm){ (Ful1-PULSE05)/PULSE05}),
         "PULSE05 to OLS"  = pmap_dbl(.l=list(OLS,PULSE05,pm), .f=function(OLS,PULSE05,pm){ (OLS-PULSE05)/PULSE05}),
         MinEigenMeanGn = map_dbl(.x= MeanGn,.f= function(x) {ifelse(is.double(eigen(x)$values),min(eigen(x)$values),NA)}),
         MaxEigenMeanGn = map_dbl(.x= MeanGn,.f= function(x) {ifelse(is.double(eigen(x)$values),max(eigen(x)$values),NA)})) %>%
  mutate(Superior = ifelse(Ful4<PULSE05,"Ful4","PULSE")) %>% 
  gather(Type,Value,c(-nModel,-nSim,-n,-MeanGn,-pm,-Ful1,-Ful4,-OLS,-PULSE05,-MinEigenMeanGn,-MaxEigenMeanGn,-Superior))


PlotData <- left_join(left_join(LossData,Optimal,by=c("n","nModel")) ,Coefs,by=c("nModel"))

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

ggsave("Plots/AllRandom_Beta00_20200627014219.png", plot = last_plot(), device = NULL, path = NULL,
       scale = 1, width = 12, height = 9, units = c("in"),
       dpi = 200, limitsize = FALSE)

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


dat <- readRDS(file="Data/FinalSim_MSE_IV2d_AllRandomCoefs_SuperiorModels_nSim_25000_nObsPerSim_50_20200702134109.RDS") %>% filter(Type != "PULSE10") #SuperiorModelswith 25k sims

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
    "rhosq"),456)
    
  ) %>% 
  unnest(cols=c(ModelCoefs)) %>% 
  spread(Coef,ModelCoefs) %>% 
  rename(RhoSq=rhosq)

LossDataSup <- dat  %>%
  select(nModel,nSim,n,Type,MeanGn,Determinant,Trace,BiasTwoNorm) %>% 
  gather(pm, Value, c("Determinant",
                      "Trace",
                      "BiasTwoNorm")) %>% 
  mutate(pm = factor(pm, levels = c("Determinant",
                                    "Trace",
                                    "BiasTwoNorm"))) %>% 
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


ggsave(file="Plots/AllRandom_Beta00_TrueSuperior_20200627014219.png", plot = plot, device = NULL, path = NULL,
       scale = 1, width = 12, height = 15, units = c("in"),
       dpi = 300, limitsize = FALSE)


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

ggsave("Plots/AllRandom_Beta00_PmSuperior_20200627014219.png", plot = last_plot(), device = NULL, path = NULL,
       scale = 1, width = 12, height = 9, units = c("in"),
       dpi = 200, limitsize = FALSE)


############# PULSE10 BETA 0 0 ##################


dat <- readRDS(file="Data/FinalSim_MSE_IV2d_AllRandomCoefs_nSim_5000_nObsPerSim_50_nModel_10000_20200701132236.RDS") #RandomVar+beta00

dat %<>% mutate(Type = ifelse(Type=="PULSE05","PULSE10",Type))

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
  select(nModel,nSim,n,Type,MeanGn,Determinant,Trace,BiasTwoNorm) %>% 
  gather(pm, Value, c("Determinant",
                      "Trace",
                      "BiasTwoNorm")) %>% 
  mutate(pm = factor(pm, levels = c("Determinant",
                                    "Trace",
                                    "BiasTwoNorm"))) %>% 
  spread(Type,Value) %>%
  ungroup() %>% 
  mutate("PULSE10 to Fuller4" = pmap_dbl(.l=list(Ful4,PULSE10,pm), .f=function(Ful4,PULSE10,pm){ (Ful4-PULSE10)/PULSE10 }),
         "PULSE10 to Fuller1"  = pmap_dbl(.l=list(Ful1,PULSE10,pm), .f=function(Ful1,PULSE10,pm){ (Ful1-PULSE10)/PULSE10}),
         "PULSE10 to OLS"  = pmap_dbl(.l=list(OLS,PULSE10,pm), .f=function(OLS,PULSE10,pm){ (OLS-PULSE10)/PULSE10}),
         MinEigenMeanGn = map_dbl(.x= MeanGn,.f= function(x) {ifelse(is.double(eigen(x)$values),min(eigen(x)$values),NA)}),
         MaxEigenMeanGn = map_dbl(.x= MeanGn,.f= function(x) {ifelse(is.double(eigen(x)$values),max(eigen(x)$values),NA)})) %>%
  mutate(Superior = ifelse(Ful4<PULSE10,"Ful4","PULSE")) %>% 
  gather(Type,Value,c(-nModel,-nSim,-n,-MeanGn,-pm,-Ful1,-Ful4,-OLS,-PULSE10,-MinEigenMeanGn,-MaxEigenMeanGn,-Superior))

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

ggsave("Plots/AllRandom_Beta00_PULSE10_20200701132236.png", plot = last_plot(), device = NULL, path = NULL,
       scale = 1, width = 12, height = 9, units = c("in"),
       dpi = 200, limitsize = FALSE)




############# BETA 1 1 ##################

dat <- readRDS(file="Data/FinalSim_MSE_IV2d_AllRandomCoefs_nSim_5000_nObsPerSim_50_nModel_10000_20200626125859.RDS") #RandomVar + beta11


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
  select(nModel,nSim,n,Type,MeanGn,Determinant,Trace,BiasTwoNorm) %>% 
  gather(pm, Value, c("Determinant",
                      "Trace",
                      "BiasTwoNorm")) %>% 
  mutate(pm = factor(pm, levels = c("Determinant",
                                    "Trace",
                                    "BiasTwoNorm"))) %>% 
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

ggsave("Plots/AllRandom_Beta11_20200626125859.png", plot = last_plot(), device = NULL, path = NULL,
       scale = 1, width = 12, height = 9, units = c("in"),
       dpi = 200, limitsize = FALSE)

############# BETA -1 1 ##################

dat <- readRDS(file="Data/FinalSim_MSE_IV2d_AllRandomCoefs_nSim_5000_nObsPerSim_50_nModel_10000_20200702230529.RDS") #RandomVar + beta-11


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
  select(nModel,nSim,n,Type,MeanGn,Determinant,Trace,BiasTwoNorm) %>% 
  gather(pm, Value, c("Determinant",
                      "Trace",
                      "BiasTwoNorm")) %>% 
  mutate(pm = factor(pm, levels = c("Determinant",
                                    "Trace",
                                    "BiasTwoNorm"))) %>% 
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

ggsave("Plots/AllRandom_Beta-11_20200626125859.png", plot = last_plot(), device = NULL, path = NULL,
       scale = 1, width = 12, height = 9, units = c("in"),
       dpi = 200, limitsize = FALSE)


############################################################
################### Experiment 2    ########################
############# 2d All random - Fixed Cor ####################
############################################################

#dat <- readRDS(file="Data/FinalSim_MSE_IV2d_AllRandomCoefsFixedNoiseCorr_nSim_3000_nObsPerSim_50_nModel_2500_20200627020030.RDS") #RandomVar+fixedcorr+alot

#dat <- readRDS(file="Data/FinalSim_MSE_IV2d_AllRandomCoefsFixedNoiseCorr_nSim_5000_nObsPerSim_50_nModel_5000_20200629104506.RDS") #RandomVar+fixedcorr+selected
#dat1 <- readRDS(file="Data/FinalSim_MSE_IV2d_AllRandomCoefsFixedNoiseCorr_nSim_5000_nObsPerSim_50_nModel_5000_20200629173217.RDS") #RandomVar+fixedcorr+12
#dat2 <- readRDS(file="Data/FinalSim_MSE_IV2d_AllRandomCoefsFixedNoiseCorr_nSim_5000_nObsPerSim_50_nModel_5000_20200629182444.RDS") #RandomVar+fixedcorr+13
#dat3 <- readRDS(file="Data/FinalSim_MSE_IV2d_AllRandomCoefsFixedNoiseCorr_nSim_5000_nObsPerSim_50_nModel_5000_20200629212227.RDS") #RandomVar+fixedcorr+14
#dat4 <- readRDS(file="Data/FinalSim_MSE_IV2d_AllRandomCoefsFixedNoiseCorr_nSim_5000_nObsPerSim_50_nModel_5000_20200629214756.RDS") #RandomVar+fixedcorr+15
#dat5 <- readRDS(file="Data/FinalSim_MSE_IV2d_AllRandomCoefsFixedNoiseCorr_nSim_5000_nObsPerSim_50_nModel_5000_20200629230200.RDS") #RandomVar+fixedcorr+16
#dat <- bind_rows(dat,dat1,dat2,dat3,dat4,dat5)

dat <- readRDS(file="Data/FinalSim_MSE_IV2d_AllRandomCoefsFixedNoiseCorr_nSim_5000_nObsPerSim_50_nModel_5000_20200701031153.RDS") #RandomVar+fixedcorr+selectedfinal



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
    "eta"),25000)  ) %>% 
  # filter(Coef %in% c("phi1","phi2","eta")) %>% 
  unnest(cols=c(ModelCoefs)) %>% 
  spread(Coef,ModelCoefs) %>% 
  mutate(normRho = (phi1^2+phi2^2-2*eta*phi1*phi2)/(1-eta^2))


LossData <- dat  %>%
  select(nModel,Cov,nSim,n,Type,MeanGn,Determinant,Trace,BiasTwoNorm) %>% 
  gather(pm, Value, c("Determinant",
                      "Trace",
                      "BiasTwoNorm")) %>% 
  mutate(pm = factor(pm, levels = c("Determinant",
                                    "Trace",
                                    "BiasTwoNorm"))) %>% 
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
  mutate(sum = Better+Worse, fractionBetter = round(100*Better/5000,1)) %>% arrange(pm)


PlotData <- left_join(left_join(LossData,Optimal,by=c("n","nModel","Cov")) ,Cors,by=c("nModel","Cov")) %>%  
  mutate(normRho = round(normRho,1)) %>% 
  arrange(normRho,-eta) %>% 
  mutate(Label = paste0("'||'*rho*'||'[2]:'",sprintf("%.3f",normRho),"'~~eta:'",sprintf("%.3f",eta),"'~~phi[1]:'",sprintf("%.3f",phi1),"'~~phi[2]:'",sprintf("%.3f",phi2),"'")) 

scaleFUN <- function(x) sprintf("%.2f", x)

p1 <- ggplot(data=PlotData %>% filter(Type=="PULSE05 to Fuller4",Cov %in% c(4))) +
  geom_hline(yintercept =0,color="black",linetype="solid") +
  geom_vline(xintercept =log(15.5),color="black",linetype="dotted") +
  #geom_point(aes(x=log(MinEigenMeanGn),y=Value,color=sqrt(MaxEigenMeanGn)),alpha=0.5,size=1)+
  geom_point(aes(x=log(MinEigenMeanGn),y=Value),alpha=0.1,size=1)+
  # scale_color_gradient(low = "blue", high = "orange")+
  facet_wrap( Label ~ pm ,scales="free_y",ncol=3, labeller = label_parsed)+
  xlab(expression(log(lambda[min](hat(E)[N](G[n])))))+
  # labs(colour="")+
  ylab(expression(paste("Relative change of performance measure")))+
  theme(plot.margin = unit(c(0,0.8,0,0), "cm"))+ 
  theme(axis.title.y = element_blank())+
  theme(axis.title.x = element_blank())+ scale_y_continuous(labels=scaleFUN)


p2 <- ggplot(data=PlotData %>% filter(Type=="PULSE05 to Fuller4",Cov %in% c(3))) +
  geom_hline(yintercept =0,color="black",linetype="solid") +
  geom_vline(xintercept =log(15.5),color="black",linetype="dotted") +
  #geom_point(aes(x=log(MinEigenMeanGn),y=Value,color=sqrt(MaxEigenMeanGn)),alpha=0.5,size=1)+
  geom_point(aes(x=log(MinEigenMeanGn),y=Value),alpha=0.1,size=1)+
  #scale_color_gradient(low = "blue", high = "orange")+
  facet_wrap(  Label ~ pm ,scales="free_y",ncol=3, labeller = label_parsed)+
  xlab(expression(log(lambda[min](hat(E)[N](G[n])))))+
  #labs(colour="")+
  ylab(expression(paste("Relative change of performance measure")))+
  theme(plot.margin = unit(c(0,0.8,0,0), "cm"))+ 
  theme(axis.title.y = element_blank())+ 
  theme(axis.title.x = element_blank())+ scale_y_continuous(labels=scaleFUN)


p3 <- ggplot(data=PlotData %>% filter(Type=="PULSE05 to Fuller4",Cov %in% c(2))) +
  geom_hline(yintercept =0,color="black",linetype="solid") +
  geom_vline(xintercept =log(15.5),color="black",linetype="dotted") +
  #geom_point(aes(x=log(MinEigenMeanGn),y=Value,color=sqrt(MaxEigenMeanGn)),alpha=0.5,size=1)+
  geom_point(aes(x=log(MinEigenMeanGn),y=Value),alpha=0.1,size=1)+
  #scale_color_gradient(low = "blue", high = "orange")+
  facet_wrap(  Label ~ pm ,scales="free_y",ncol=3, labeller = label_parsed)+
  xlab(expression(log(lambda[min](hat(E)[N](G[n])))))+
  #labs(colour="")+
  ylab(expression(paste("Relative change of performance measure")))+
  theme(plot.margin = unit(c(0,0.8,0,0), "cm"))+ 
  theme(axis.title.y = element_blank())+ 
  theme(axis.title.x = element_blank())+ scale_y_continuous(labels=scaleFUN)

p4 <- ggplot(data=PlotData %>% filter(Type=="PULSE05 to Fuller4",Cov %in% c(5))) +
  geom_hline(yintercept =0,color="black",linetype="solid") +
  geom_vline(xintercept =log(15.5),color="black",linetype="dotted") +
  #geom_point(aes(x=log(MinEigenMeanGn),y=Value,color=sqrt(MaxEigenMeanGn)),alpha=0.5,size=1)+
  geom_point(aes(x=log(MinEigenMeanGn),y=Value),alpha=0.1,size=1)+
  #scale_color_gradient(low = "blue", high = "orange")+
  facet_wrap(  Label ~ pm ,scales="free_y",ncol=3, labeller = label_parsed)+
  xlab(expression(log(lambda[min](hat(E)[N](G[n])))))+
  #labs(colour="")+
  ylab(expression(paste("Relative change of performance measure")))+
  theme(plot.margin = unit(c(0,0.8,0,0), "cm"))+ theme(axis.title.y = element_blank())+ scale_y_continuous(labels=scaleFUN)

grid.arrange(p1,p2,p3,p4,ncol=1)

plot <- arrangeGrob(p1,p2,p3,p4, ncol=1,left = textGrob("Relative Change in Performance Measure", rot = 90, vjust = 0))
#,right = textGrob(expression(sqrt(lambda[max](hat(E)[N](G[n])))), rot = 270, vjust = 2))



ggsave(file="Plots/AllRandom_fixedCorr_SMallRhoLargeEtaDecrease.png", plot = plot, device = NULL, path = NULL,
       scale = 1, width = 12, height = 12, units = c("in"),
       dpi = 300, limitsize = FALSE)


p1 <- ggplot(data=PlotData %>% filter(Type=="PULSE05 to Fuller4",Cov %in% c(4))) +
  geom_hline(yintercept =0,color="black",linetype="solid") +
  geom_vline(xintercept =log(15.5),color="black",linetype="dotted") +
  geom_point(aes(x=log(MinEigenMeanGn),y=Value,color=sqrt(MaxEigenMeanGn)),alpha=0.5,size=0.8)+
  scale_color_gradient(low = "blue", high = "orange")+
  facet_wrap( Label ~ pm ,scales="free_y",ncol=3, labeller = label_parsed)+
  xlab(expression(log(lambda[min](E[N](G[n])))))+
  labs(colour="")+
  ylab(expression(paste("Relative change of performance measure")))+
  theme(plot.margin = unit(c(0,0.8,0,0), "cm"))+ 
  theme(axis.title.y = element_blank())+ 
  theme(axis.title.x = element_blank())+ scale_y_continuous(labels=scaleFUN)

p2 <- ggplot(data=PlotData %>% filter(Type=="PULSE05 to Fuller4",Cov %in% c(1))) +
  geom_hline(yintercept =0,color="black",linetype="solid") +
  geom_vline(xintercept =log(15.5),color="black",linetype="dotted") +
  geom_point(aes(x=log(MinEigenMeanGn),y=Value,color=sqrt(MaxEigenMeanGn)),alpha=0.5,size=0.8)+
  scale_color_gradient(low = "blue", high = "orange")+
  facet_wrap(  Label ~ pm ,scales="free_y",ncol=3, labeller = label_parsed)+
  xlab(expression(log(lambda[min](E[N](G[n])))))+
  labs(colour="")+
  ylab(expression(paste("Relative change of performance measure")))+
  theme(plot.margin = unit(c(0,0.8,0,0), "cm"))+ theme(axis.title.y = element_blank())+ scale_y_continuous(labels=scaleFUN)


plot <- arrangeGrob(p1,p2, ncol=1,widths=c(1),left = textGrob("Relative Change in Performance Measure", rot = 90, vjust = 0),
                    right = textGrob(expression(sqrt(lambda[max](hat(E)[N](G[n])))), rot = 270, vjust = 2))

ggsave(file="Plots/AllRandom_fixedCorr_HighEtaSymmetricToNonSymmetricPhi.png", plot = plot, device = NULL, path = NULL,
       scale = 1, width = 12, height = 6, units = c("in"),
       dpi = 300, limitsize = FALSE)

