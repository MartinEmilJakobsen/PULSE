library(MASS) # mvrnorm
library(matrixcalc)
library(stats)
library(lqmm) #makepositive.
library(expm) #sqrtm A=SS
library(furrr)
library(tidyverse)
library(metR)
library(gridExtra)
library(dplyr)

source("Estimators_Fast.R")

i <- 0.1  
a <- i*matrix(c(5,2,2,3),ncol=2)
delta <- matrix(c(1,3,2,10),ncol=2)
beta <- matrix(c(1,1),ncol=1)
mu <- matrix(c(15,4),ncol=1)

set.seed(5)
n <- 500
dA <- 2
dX <- 2
dZ <- 2
A <- mvrnorm(n,mu=c(0,0),Sigma=matrix(c(1,0,0,1),ncol=2,byrow=TRUE))
H <- mvrnorm(n,mu=c(0,0),Sigma=matrix(c(1,0,0,1),ncol=2,byrow=TRUE))
X <- A%*%a+H%*%delta + mvrnorm(n,mu=c(0,0),Sigma=matrix(c(1,0,0,1),ncol=2,byrow=TRUE))
Y <- X%*%beta+H%*%mu + rnorm(n)
Z <- cbind(X)


P_A <- A%*%solve(t(A)%*%A)%*%t(A)
YtP_AY <- t(Y)%*%P_A%*%Y
ZtP_AZ <- t(Z)%*%P_A%*%Z
YtP_AZ <- t(Y)%*%P_A%*%Z
YtY <- t(Y)%*%Y
ZtZ <- t(Z)%*%Z
YtZ <- t(Y)%*%Z


alphaPval <- BinarySearch(Z,Y,dZ,dA,p=0.05,N=1000000,n,YtP_AY,ZtP_AZ,YtP_AZ,YtY,ZtZ,YtZ,P_A,A_1,X)


lOLS <- function(a1,a2) {
  return(n^(-1)*t(Y-Z%*%c(a1,a2))%*%(Y-Z%*%c(a1,a2)))
}
lIV <- function(a1,a2) {
  return(n^(-1)*t(Y-Z%*%c(a1,a2))%*%P_A%*%(Y-Z%*%c(a1,a2)))
}
T <- function(a1,a2) { n*lIV(a1,a2)/lOLS(a1,a2)}

dat1 <- expand_grid(alpha1=seq(-4,7,0.01),alpha2=seq(-2,5,0.01)) %>% 
  rowwise() %>% 
  mutate(lOLS = lOLS(alpha1,alpha2),
         lIV = lIV(alpha1,alpha2),
         Test = n*lIV/lOLS) %>% 
  gather(Type,Value,c(lOLS,lIV,Test))

alphaOLS = solve(t(Z)%*%Z)%*%t(Z)%*%Y
datOLS <- data.frame(alpha1=alphaOLS[1],alpha2=alphaOLS[2])

alpha2SLS <- K_class(1,Z,Y,n,P_A)
dat2SLS <- data.frame(alpha1=alpha2SLS[1],alpha2=alpha2SLS[2])

datPval <- data.frame(alpha1=alphaPval[1],alpha2=alphaPval[2])

lOLSinPval <- lOLS(alphaPval[1],alphaPval[2])
lIVinPval <- lIV(alphaPval[1],alphaPval[2])

datTrue <-  data.frame(alpha1=1,alpha2=1)


############################
#### All plots together ####
############################


cols <- c("OLS" = "red","Test" = "blue","TSLS" = "green4", "AR" = "dodgerblue1", "Path" = "black")
q <- qchisq(0.95,df=2)
datSublevel <- dat1 %>% filter(Type=="Test") %>% filter(Value <=5.991465)


K_class_fixed <- function(k){
  W <- t(Z) %*% ((1-k)*diag(n)+k* P_A )  %*%Z
  if(is.singular.matrix(W, tol = 1e-20))
  {a <- "Inversion of Singular Matrix"}
  else{
    a <- solve(W) %*% t(Z) %*% ( (1-k)*diag(n) + k*P_A ) %*% Y 
  }
  return(a)
}


K_class_path <- data.frame(k = seq(0,1,0.001)) %>% 
  rowwise() %>% 
  mutate(a = list(K_class_fixed(k))) %>% 
  ungroup() %>% 
  mutate(ID = seq(1,nrow(.),1)) %>% 
  unnest(cols=c(a)) %>% 
  group_by(ID) %>% 
  mutate(ID2 = c("alpha1","alpha2")) %>% 
  ungroup %>% 
  spread(ID2,a)


#t text
K_class_path_text <-data.frame(k = c(0.9366,0.9921)) %>% 
  rowwise() %>% 
  mutate(a = list(K_class_fixed(k))) %>% 
  ungroup() %>% 
  mutate(ID = seq(1,nrow(.),1)) %>% 
  unnest(cols=c(a)) %>% 
  group_by(ID) %>% 
  mutate(ID2 = c("alpha1","alpha2")) %>% 
  ungroup %>% 
  spread(ID2,a) %>% 
  rowwise() %>%  
  mutate(t = as.character(round(lIV(alpha1,alpha2),2))) 

K_class_path_text <- bind_rows(K_class_path_text, 
                               data.frame(k=NA,ID=NA,
                                          alpha1 =alphaPval[1] , 
                                          alpha2=alphaPval[2] , 
                                          t="t*(p)"))

p1<- ggplot(data=dat1,aes(alpha1,alpha2))+
  geom_raster(data=datSublevel, aes(x=alpha1, y=alpha2,fill="AR"),interpolate = TRUE,alpha=0.1)+
  scale_fill_manual(values = cols,
                    name="Acceptance Region:",
                    labels=c(""))+
  geom_contour(data=dat1 %>% filter(Type=="Test"),aes(z=Value,colour="Test"),breaks=c(2,4,5.99, 8,10,12,14,16,18,20))+
  geom_contour(data=dat1 %>% filter(Type=="lIV"),aes(z=Value,colour="TSLS"),breaks=c(1,lIVinPval,3))+
  geom_text_contour(data=dat1 %>% filter(Type=="Test") %>% mutate(Value=round(Value,digits=1)),aes(z = Value,group=Type),color="black")+
  geom_contour(data=dat1 %>% filter(Type=="lOLS"),aes(z=Value,colour="OLS"),breaks=c(180,200,250,500))+
  geom_path(data=K_class_path,aes(x=alpha1,y=alpha2,colour="Path"))+
  geom_point(data=K_class_path_text ,aes(x=alpha1,y=alpha2),color="black",size=1)+
  geom_text(data=K_class_path_text ,aes(x=alpha1,y=alpha2,label=t),hjust=-0.1, vjust=-0.2)+
  geom_point(data=datOLS,color="red",size=1)+
  geom_point(data=dat2SLS, color ="green4",size=1)+
  geom_point(data=datPval, color ="black",size=1)+
  scale_colour_manual(values = cols,
                      name = 'Levelsets and primal path:',
                      breaks = c("OLS","Test","TSLS","Path"),
                      labels=c(expression(~l[OLS]^n*(alpha)),
                               expression(~T[n](alpha)),
                               expression(~l[IV]^n*(alpha)),
                               expression(  "      "~bgroup("{",hat(alpha)[Pr]^n*(t):t%in%D[Pr],"}"))))+
  xlab(expression(alpha[1]))+
  ylab(expression(alpha[2]))+ 
  theme_bw()+
  theme(legend.position="bottom")+
  scale_x_continuous(expand = c(0,0))+
  scale_y_continuous(expand = c(0,0))



p2 <-ggplot(data=dat1,aes(alpha1,alpha2))+
  geom_raster(data=datSublevel, aes(x=alpha1, y=alpha2,fill="AR"),interpolate = TRUE,alpha=0.1)+
  scale_fill_manual(values = cols,
                    name="Acceptance Region:",
                    labels=c(""))+
  geom_contour(data=dat1 %>% filter(Type=="Test"),aes(z=Value,colour="Test"),breaks=c(2,4,5.99, 8,10,12,14,16,18,20))+
  geom_contour(data=dat1 %>% filter(Type=="lIV"),aes(z=Value,colour="TSLS"),breaks=c(1,lIVinPval,3))+
  geom_text_contour(data=dat1 %>% filter(Type=="Test") %>% mutate(Value=round(Value,digits=1)),aes(z = Value,group=Type),color="black")+
  geom_contour(data=dat1 %>% filter(Type=="lOLS"),aes(z=Value,colour="OLS"),breaks=c(lOLSinPval,180,200,250,500))+
  geom_path(data=K_class_path,aes(x=alpha1,y=alpha2,colour="Path"))+
  geom_point(data=K_class_path_text ,aes(x=alpha1,y=alpha2),color="black",size=1)+
  geom_text(data=K_class_path_text ,aes(x=alpha1,y=alpha2,label=t),hjust=-0.1, vjust=-0.2)+
  geom_point(data=datOLS,color="red",size=1)+
  geom_point(data=dat2SLS, color ="green4",size=1)+
  geom_point(data=datPval, color ="black",size=1)+
  scale_colour_manual(values = cols,
                      name = 'Levelsets and primal path:',
                      breaks = c("OLS","Test","TSLS","Path"),
                      labels=c(expression(~l[OLS]^n*(alpha)),
                               expression(~T[n](alpha)),
                               expression(~l[IV]^n*(alpha)),
                               expression(  "     "~bgroup("{",hat(alpha)[Pr]^n*(t):t%in%D[Pr],"}"))))+
  xlab(expression(alpha[1]))+
  ylab(expression(alpha[2]))+ 
  theme_bw()+
  theme(legend.position="bottom")+
  scale_x_continuous(expand = c(0,0),limits=c(3,5))+
  scale_y_continuous(expand = c(0,0),limits= c(0,1.5))

g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

mylegend<-g_legend(p1)


p3 <- grid.arrange(arrangeGrob(p1 + theme(legend.position="none"),
                               p2 + theme(legend.position="none"),
                               nrow=1),
                   mylegend, nrow=2,heights=c(10, 0.5))


ggsave(filename="Plots/Levelsets_Test_OLS_IV_Combined.pdf", plot = p3, device = NULL, path = NULL,
       scale = 1, width = 12, height = 6, units = c("in"),
       dpi = 200)
ggsave(filename="Plots/Levelsets_Test_OLS_IV_Combined.eps", plot = p3, device = cairo_ps, path = NULL,
       scale = 1, width = 12, height = 6, units = c("in"),
       dpi = 200)
