library(tidyverse)
library(magrittr)
library(stringr)
library(gridExtra)
library(dplyr)

#####################################
####### Univariate Experiment ####### 
#######         Plots         ####### 
#####################################

# Read data
Data_Location <- "Data/Experiment_Univariate_nSim_15000_nObsPerSim_50_100_150_20200716034113.RDS"
ID <- "20200716034113"
dat <- readRDS(file=Data_Location) 

# Select relevant data
dat2 <- dat %>% 
  rowwise() %>% 
  mutate(Conc = n*Rsq/(1-Rsq), Fstat = Conc/nInst+1) %>%   
  select(rho,n,nInst,Rsq,Type,Conc,Fstat,
         MeanGn,MedianBias,MeanEstimate,
         RMSE,VarEstimate,IQRestimate) %>% 
  arrange(rho,Rsq,nInst,n) 

#F test lower rejection thresholds for relevancy of instruments

n <- c(50,100,150)
q <- c(1,2,3,4,5,10,20,30)

expand_grid(n= n,q =q ) %>% 
  mutate(quantile= qf(c(0.95),df1=q,df2=n-q)) %>% 
  summarise(min=min(quantile),max=max(quantile)) %>% 
  mutate(logmin=log(min),logmax=log(max))

####################
#### RMSE Plots ####
####################

RMSE_plotdata <- dat2 %>% 
  select(-MedianBias,-IQRestimate,-MeanEstimate,-VarEstimate) %>%  
  spread(Type,RMSE) %>% 
  arrange(rho,nInst,Rsq,n,nInst) %>%  
  mutate("PULSE05 to Fuller1" = (Ful1-PULSE05)/PULSE05,
         "PULSE05 to Fuller4" = (Ful4-PULSE05)/PULSE05,
         "PULSE05 to TSLS" = (TSLS-PULSE05)/PULSE05,
         "PULSE05 to OLS" =(OLS-PULSE05)/PULSE05 
  )  %>% 
  select(-Ful1,-Ful4,-OLS,-PULSE05,-TSLS,-LIML)  %>% 
  gather(Type,Value,c(-Fstat,-rho,-nInst,-Conc,-Rsq,-n,-MeanGn )) %>% 
  mutate(rho=factor(rho), n = factor(n)) %>% 
  filter(!(Type=="PULSE05 to TSLS" & nInst <3 )) 

dummy <- data.frame(MeanGn=2,Value=-1)

ggplot(data = RMSE_plotdata %>% filter(!(rho %in% c("0.01","0.05")) ) )+
  geom_point(aes(x=log(MeanGn),y=Value,color=rho,shape=n)) +
  geom_vline(xintercept =log(10),color="black",linetype="dotted") +
  geom_vline(xintercept =log(1.55),color="black",linetype="dashed") +
 # geom_blank(data=dummy,aes(x=MeanGn,y=Value))+
  facet_wrap(~Type,scales= "free")+
  xlab(expression(log(hat(E)[N](G[n]))))+
  labs(colour=expression(rho))+
  ylab(expression(paste("Relative Change in RMSE")))+
  theme(plot.margin = unit(c(0,0.8,0,0), "cm"))

ggsave(paste0("Plots/Univariate_RMSE_",ID,".png"), plot = last_plot(), device = NULL, path = NULL,
       scale = 1, width = 12, height = 6, units = c("in"),
       dpi = 200, limitsize = FALSE)

########################
#### MeanBias Plots ####
########################

MeanBias_plotdata <- dat2 %>% 
  mutate(Bias = abs(MeanEstimate)) %>% 
  select(-MedianBias,-IQRestimate,-RMSE,-VarEstimate,-MeanEstimate) %>%
  spread(Type,Bias) %>% 
  arrange(rho,nInst,Rsq,n,nInst) %>%  
  mutate("PULSE05 to Fuller1" = (Ful1-PULSE05)/PULSE05,
         "PULSE05 to Fuller4" = (Ful4-PULSE05)/PULSE05,
         "PULSE05 to TSLS" = (TSLS-PULSE05)/PULSE05,
         "PULSE05 to OLS" =(OLS-PULSE05)/PULSE05 
  )  %>% 
  select(-Ful1,-Ful4,-OLS,-PULSE05,-TSLS,-LIML)  %>% 
  gather(Type,Value,c(-Fstat,-rho,-nInst,-Conc,-Rsq,-n,-MeanGn )) %>% 
  mutate(rho=factor(rho), n = factor(n)) %>% 
  filter(!(Type=="PULSE05 to TSLS" & nInst <2 ))


ggplot(data = MeanBias_plotdata)+
  geom_point(aes(x=log(MeanGn),y=Value,color=rho,shape=n)) +
  geom_vline(xintercept =log(10),color="black",linetype="dotted") +
  geom_vline(xintercept =log(1.55),color="black",linetype="dashed") +
  facet_wrap(~Type,scales= "free")+
  xlab(expression(log(hat(E)[N](G[n]))))+
  labs(colour=expression(rho))+
  ylab(expression(paste("Relative Change in Mean Bias")))+
  theme(plot.margin = unit(c(0,0.8,0,0), "cm"))


ggsave(paste0("Plots/Univariate_MeanBias_",ID,".png"), plot = last_plot(), device = NULL, path = NULL,
       scale = 1, width = 12, height = 6, units = c("in"),
       dpi = 200, limitsize = FALSE)

##########################
#### MedianBias Plots ####
##########################

MedianBias_plotdata <- dat2 %>% 
  mutate(MedianBias = abs(MedianBias)) %>% 
  select(-MeanEstimate,-IQRestimate,-RMSE,-VarEstimate,-MeanEstimate) %>%
  spread(Type,MedianBias) %>% 
  arrange(rho,nInst,Rsq,n,nInst) %>%  
  mutate("PULSE05 to Fuller1" = (Ful1-PULSE05)/PULSE05,
         "PULSE05 to Fuller4" = (Ful4-PULSE05)/PULSE05,
         "PULSE05 to TSLS" = (TSLS-PULSE05)/PULSE05,
         "PULSE05 to OLS" =(OLS-PULSE05)/PULSE05 
  )  %>% 
  select(-Ful1,-Ful4,-OLS,-PULSE05,-TSLS,-LIML)  %>% 
  gather(Type,Value,c(-Fstat,-rho,-nInst,-Conc,-Rsq,-n,-MeanGn )) %>% 
  mutate(rho=factor(rho), n = factor(n)) %>% 
  filter(!(Type=="PULSE05 to TSLS" & nInst <2 ))


ggplot(data = MedianBias_plotdata)+
  geom_point(aes(x=log(MeanGn),y=Value,color=rho,shape=n)) +
  geom_vline(xintercept =log(10),color="black",linetype="dotted") +
  geom_vline(xintercept =log(1.55),color="black",linetype="dashed") +
  facet_wrap(~Type,scales= "free")+
  xlab(expression(log(hat(E)[N](G[n]))))+
  labs(colour=expression(rho))+
  ylab(expression(paste("Relative Change in Median Bias")))+
  theme(plot.margin = unit(c(0,0.8,0,0), "cm"))


ggsave(paste0("Plots/Univariate_MedianBias_",ID,".png"), plot = last_plot(), device = NULL, path = NULL,
       scale = 1, width = 12, height = 6, units = c("in"),
       dpi = 200, limitsize = FALSE)

########################
#### Variance Plots ####
########################

Variance_plotdata <- dat2 %>% 
  select(-MedianBias,-IQRestimate,-RMSE,-MeanEstimate) %>%  
  spread(Type,VarEstimate) %>% 
  arrange(rho,nInst,Rsq,n,nInst) %>%  
  mutate("PULSE05 to Fuller1" = (Ful1-PULSE05)/PULSE05,
         "PULSE05 to Fuller4" = (Ful4-PULSE05)/PULSE05,
         "PULSE05 to TSLS" = (TSLS-PULSE05)/PULSE05,
         "PULSE05 to OLS" =(OLS-PULSE05)/PULSE05 
  )  %>% 
  select(-Ful1,-Ful4,-OLS,-PULSE05,-TSLS,-LIML)  %>% 
  gather(Type,Value,c(-Fstat,-rho,-nInst,-Conc,-Rsq,-n,-MeanGn )) %>% 
  mutate(rho=factor(rho), n = factor(n)) %>% 
  filter(!(Type=="PULSE05 to TSLS" & nInst <3) & !(Value>=100))

ggplot(data = Variance_plotdata)+
  geom_point(aes(x=log(MeanGn),y=Value,color=rho,shape=n)) +
  geom_vline(xintercept =log(10),color="black",linetype="dotted") +
  geom_vline(xintercept =log(1.55),color="black",linetype="dashed") +
  facet_wrap(~Type,scales= "free")+
  xlab(expression(log(hat(E)[N](G[n]))))+
  labs(colour=expression(rho))+
  ylab(expression(paste("Relative Change in Variance")))+
  theme(plot.margin = unit(c(0,0.8,0,0), "cm"))

ggsave(paste0("Plots/Univariate_Variance_",ID,".png"), plot = last_plot(), device = NULL, path = NULL,
       scale = 1, width = 12, height = 6, units = c("in"),
       dpi = 200, limitsize = FALSE)


###################
#### IQR Plots ####
###################

IQR_plotdata <- dat2 %>% 
  select(-MedianBias,-VarEstimate,-RMSE,-MeanEstimate) %>%  
  spread(Type,IQRestimate) %>% 
  arrange(rho,nInst,Rsq,n,nInst) %>%  
  mutate("PULSE05 to Fuller1" = (Ful1-PULSE05)/PULSE05,
         "PULSE05 to Fuller4" = (Ful4-PULSE05)/PULSE05,
         "PULSE05 to TSLS" = (TSLS-PULSE05)/PULSE05,
         "PULSE05 to OLS" =(OLS-PULSE05)/PULSE05 
  )  %>% 
  select(-Ful1,-Ful4,-OLS,-PULSE05,-TSLS,-LIML)  %>% 
  gather(Type,Value,c(-Fstat,-rho,-nInst,-Conc,-Rsq,-n,-MeanGn )) %>% 
  mutate(rho=factor(rho), n = factor(n)) 


ggplot(data = IQR_plotdata)+
  geom_point(aes(x=log(MeanGn),y=Value,color=rho,shape=n)) +
  geom_vline(xintercept =log(10),color="black",linetype="dotted") +
  geom_vline(xintercept =log(1.55),color="black",linetype="dashed") +
  facet_wrap(~Type,scales= "free")+
  xlab(expression(log(hat(E)[N](G[n]))))+
  labs(colour=expression(rho))+
  ylab(expression(paste("Relative Change in IQR")))+
  theme(plot.margin = unit(c(0,0.8,0,0), "cm"))


ggsave(paste0("Plots/Univariate_IQR_",ID,".png"), plot = last_plot(), device = NULL, path = NULL,
       scale = 1, width = 12, height = 6, units = c("in"),
       dpi = 200, limitsize = FALSE)


#####################################
####### Univariate Experiment ####### 
#######         Plots         ####### 
#######  Alternative Settings #######
#####################################

# Read data
Data_Location <- "Data/Experiment_Univariate_AlternativeSetups_nSim_5000_nObsPerSim_50_100_20200716105326.RDS"
ID <- "20200716105326"
dat <- readRDS(file=Data_Location) 

# Select relevant data
dat2 <- dat %>% 
  rowwise() %>% 
  mutate(Conc = n*Rsq/(1-Rsq), Fstat = Conc/nInst+1) %>%   
  select(truealpha,xis,rho,n,nInst,Rsq,Type,Conc,Fstat,
         MeanGn,MedianBias,MeanEstimate,
         RMSE,VarEstimate,IQRestimate) %>% 
  arrange(rho,Rsq,nInst,n) 

####################
#### RMSE Plots ####
####################

RMSE_plotdata <- dat2 %>% 
  select(-MedianBias,-IQRestimate,-MeanEstimate,-VarEstimate) %>%  
  spread(Type,RMSE) %>% 
  arrange(rho,nInst,Rsq,n,nInst) %>%  
  mutate("PULSE05 to Fuller1" = (Ful1-PULSE05)/PULSE05,
         "PULSE05 to Fuller4" = (Ful4-PULSE05)/PULSE05,
         "PULSE05 to TSLS" = (TSLS-PULSE05)/PULSE05,
         "PULSE05 to OLS" =(OLS-PULSE05)/PULSE05 
  )  %>% 
  select(-Ful1,-Ful4,-OLS,-PULSE05,-TSLS,-LIML)  %>% 
  gather(Type,Value,c(-truealpha,-xis,-Fstat,-rho,-nInst,-Conc,-Rsq,-n,-MeanGn )) %>% 
  mutate(rho=factor(rho), n = factor(n)) %>% 
  filter(!(Type=="PULSE05 to TSLS" & nInst <3 )) 

dummy <- data.frame(MeanGn=2,Value=-1)

ggplot(data = RMSE_plotdata %>% filter(!(rho %in% c("0.01","0.05")) ) )+
  geom_point(aes(x=log(MeanGn),y=Value,color=rho,shape=n)) +
  geom_vline(xintercept =log(10),color="black",linetype="dotted") +
  geom_vline(xintercept =log(1.55),color="black",linetype="dashed") +
  # geom_blank(data=dummy,aes(x=MeanGn,y=Value))+
  facet_wrap(truealpha+xis~Type,scales= "free",labeller = "label_both",ncol=4)+
  xlab(expression(log(hat(E)[N](G[n]))))+
  labs(colour=expression(rho))+
  ylab(expression(paste("Relative Change in RMSE")))+
  theme(plot.margin = unit(c(0,0.8,0,0), "cm"))

ggsave(paste0("Plots/Univariate_AlternativeSetups_RMSE_",ID,".png"), plot = last_plot(), device = NULL, path = NULL,
       scale = 1, width = 12, height = 18, units = c("in"),
       dpi = 200, limitsize = FALSE)

########################
#### MeanBias Plots ####
########################

MeanBias_plotdata <- dat2 %>% 
  mutate(Bias = abs(MeanEstimate)) %>% 
  select(-MedianBias,-IQRestimate,-RMSE,-VarEstimate,-MeanEstimate) %>%
  spread(Type,Bias) %>% 
  arrange(rho,nInst,Rsq,n,nInst) %>%  
  mutate("PULSE05 to Fuller1" = (Ful1-PULSE05)/PULSE05,
         "PULSE05 to Fuller4" = (Ful4-PULSE05)/PULSE05,
         "PULSE05 to TSLS" = (TSLS-PULSE05)/PULSE05,
         "PULSE05 to OLS" =(OLS-PULSE05)/PULSE05 
  )  %>% 
  select(-Ful1,-Ful4,-OLS,-PULSE05,-TSLS,-LIML)  %>% 
  gather(Type,Value,c(-truealpha,-xis,-Fstat,-rho,-nInst,-Conc,-Rsq,-n,-MeanGn )) %>% 
  mutate(rho=factor(rho), n = factor(n)) %>% 
  filter(!(Type=="PULSE05 to TSLS" & nInst <2 ))


ggplot(data = MeanBias_plotdata)+
  geom_point(aes(x=log(MeanGn),y=Value,color=rho,shape=n)) +
  geom_vline(xintercept =log(10),color="black",linetype="dotted") +
  geom_vline(xintercept =log(1.55),color="black",linetype="dashed") +
  facet_wrap(truealpha+xis~Type,scales= "free",labeller = "label_both",ncol=4)+
  xlab(expression(log(hat(E)[N](G[n]))))+
  labs(colour=expression(rho))+
  ylab(expression(paste("Relative Change in Mean Bias")))+
  theme(plot.margin = unit(c(0,0.8,0,0), "cm"))


ggsave(paste0("Plots/Univariate_AlternativeSetups_MeanBias_",ID,".png"), plot = last_plot(), device = NULL, path = NULL,
       scale = 1, width = 12, height = 18, units = c("in"),
       dpi = 200, limitsize = FALSE)

##########################
#### MedianBias Plots ####
##########################

MedianBias_plotdata <- dat2 %>% 
  mutate(MedianBias = abs(MedianBias)) %>% 
  select(-MeanEstimate,-IQRestimate,-RMSE,-VarEstimate,-MeanEstimate) %>%
  spread(Type,MedianBias) %>% 
  arrange(rho,nInst,Rsq,n,nInst) %>%  
  mutate("PULSE05 to Fuller1" = (Ful1-PULSE05)/PULSE05,
         "PULSE05 to Fuller4" = (Ful4-PULSE05)/PULSE05,
         "PULSE05 to TSLS" = (TSLS-PULSE05)/PULSE05,
         "PULSE05 to OLS" =(OLS-PULSE05)/PULSE05 
  )  %>% 
  select(-Ful1,-Ful4,-OLS,-PULSE05,-TSLS,-LIML)  %>% 
  gather(Type,Value,c(-truealpha,-xis,-Fstat,-rho,-nInst,-Conc,-Rsq,-n,-MeanGn )) %>% 
  mutate(rho=factor(rho), n = factor(n)) %>% 
  filter(!(Type=="PULSE05 to TSLS" & nInst <2 ))


ggplot(data = MedianBias_plotdata)+
  geom_point(aes(x=log(MeanGn),y=Value,color=rho,shape=n)) +
  geom_vline(xintercept =log(10),color="black",linetype="dotted") +
  geom_vline(xintercept =log(1.55),color="black",linetype="dashed") +
  facet_wrap(truealpha+xis~Type,scales= "free",labeller = "label_both",ncol=4)+
  xlab(expression(log(hat(E)[N](G[n]))))+
  labs(colour=expression(rho))+
  ylab(expression(paste("Relative Change in Median Bias")))+
  theme(plot.margin = unit(c(0,0.8,0,0), "cm"))


ggsave(paste0("Plots/Univariate_AlternativeSetups_MedianBias_",ID,".png"), plot = last_plot(), device = NULL, path = NULL,
       scale = 1, width = 12, height = 18, units = c("in"),
       dpi = 200, limitsize = FALSE)

########################
#### Variance Plots ####
########################

Variance_plotdata <- dat2 %>% 
  select(-MedianBias,-IQRestimate,-RMSE,-MeanEstimate) %>%  
  spread(Type,VarEstimate) %>% 
  arrange(rho,nInst,Rsq,n,nInst) %>%  
  mutate("PULSE05 to Fuller1" = (Ful1-PULSE05)/PULSE05,
         "PULSE05 to Fuller4" = (Ful4-PULSE05)/PULSE05,
         "PULSE05 to TSLS" = (TSLS-PULSE05)/PULSE05,
         "PULSE05 to OLS" =(OLS-PULSE05)/PULSE05 
  )  %>% 
  select(-Ful1,-Ful4,-OLS,-PULSE05,-TSLS,-LIML)  %>% 
  gather(Type,Value,c(-truealpha,-xis,-Fstat,-rho,-nInst,-Conc,-Rsq,-n,-MeanGn )) %>% 
  mutate(rho=factor(rho), n = factor(n)) %>% 
  filter(!(Type=="PULSE05 to TSLS" & nInst <3) & !(Value>=100))

ggplot(data = Variance_plotdata)+
  geom_point(aes(x=log(MeanGn),y=Value,color=rho,shape=n)) +
  geom_vline(xintercept =log(10),color="black",linetype="dotted") +
  geom_vline(xintercept =log(1.55),color="black",linetype="dashed") +
  facet_wrap(truealpha+xis~Type,scales= "free",labeller = "label_both",ncol=4)+
  xlab(expression(log(hat(E)[N](G[n]))))+
  labs(colour=expression(rho))+
  ylab(expression(paste("Relative Change in Variance")))+
  theme(plot.margin = unit(c(0,0.8,0,0), "cm"))

ggsave(paste0("Plots/Univariate_AlternativeSetups_Variance_",ID,".png"), plot = last_plot(), device = NULL, path = NULL,
       scale = 1, width = 12, height = 18, units = c("in"),
       dpi = 200, limitsize = FALSE)


###################
#### IQR Plots ####
###################

IQR_plotdata <- dat2 %>% 
  select(-MedianBias,-VarEstimate,-RMSE,-MeanEstimate) %>%  
  spread(Type,IQRestimate) %>% 
  arrange(rho,nInst,Rsq,n,nInst) %>%  
  mutate("PULSE05 to Fuller1" = (Ful1-PULSE05)/PULSE05,
         "PULSE05 to Fuller4" = (Ful4-PULSE05)/PULSE05,
         "PULSE05 to TSLS" = (TSLS-PULSE05)/PULSE05,
         "PULSE05 to OLS" =(OLS-PULSE05)/PULSE05 
  )  %>% 
  select(-Ful1,-Ful4,-OLS,-PULSE05,-TSLS,-LIML)  %>% 
  gather(Type,Value,c(-truealpha,-xis,-Fstat,-rho,-nInst,-Conc,-Rsq,-n,-MeanGn )) %>% 
  mutate(rho=factor(rho), n = factor(n)) 


ggplot(data = IQR_plotdata)+
  geom_point(aes(x=log(MeanGn),y=Value,color=rho,shape=n)) +
  geom_vline(xintercept =log(10),color="black",linetype="dotted") +
  geom_vline(xintercept =log(1.55),color="black",linetype="dashed") +
  facet_wrap(truealpha+xis~Type,scales= "free",labeller = "label_both",ncol=4)+
  xlab(expression(log(hat(E)[N](G[n]))))+
  labs(colour=expression(rho))+
  ylab(expression(paste("Relative Change in IQR")))+
  theme(plot.margin = unit(c(0,0.8,0,0), "cm"))


ggsave(paste0("Plots/Univariate_AlternativeSetups_IQR_",ID,".png"), plot = last_plot(), device = NULL, path = NULL,
       scale = 1, width = 12, height = 18, units = c("in"),
       dpi = 200, limitsize = FALSE)
