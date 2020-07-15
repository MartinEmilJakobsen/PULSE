CovStructure <- data.frame(
  phi1 = c(0.1897366597,  0.1549193338, 0.7589466387,  0.6196773353),
  phi2 = c(0.1897366597,  0.1549193338, 0.7589466387,  0.6196773353),
  eta =  c(0.8,           0.2,          0.8,           0.2)
) %>%  mutate(RhoNorm = sqrt((phi1^2+phi2^2-2*eta*phi1*phi2)/(1-eta^2)))




phi1 <- CovStructure[1,1]
phi2 <- CovStructure[1,2]
eta <-  CovStructure[1,3]


xi11 <-1
xi12 <-0
xi21 <-0
xi22 <-1

xi = matrix(c(xi11,xi12,xi21,xi22),ncol=2,byrow=TRUE)

n<- 1500
truealpha <- c(1,1)

  U_XY = mvrnorm(n,mu=c(0,0,0),Sigma=matrix(c(1,   eta,   phi1,
                                              eta,   1,   phi2,
                                              phi1, phi2, 1   ),ncol=3,byrow=TRUE))  
  
  A <- mvrnorm(n,mu=c(0,0),Sigma=matrix(c(1,0,0,1),ncol=2,byrow=TRUE))
  X <- A%*%xi+ U_XY[,1:2]
  Y <- X%*%matrix(truealpha,ncol=1,byrow=TRUE) + U_XY[,3]
  
  A_1 <- "none"
  Z <- X
  
  P_A <- A%*%solve(t(A)%*%A)%*%t(A)
  YtP_AY <- t(Y)%*%P_A%*%Y
  ZtP_AZ <- t(Z)%*%P_A%*%Z
  YtP_AZ <- t(Y)%*%P_A%*%Z
  YtY <- t(Y)%*%Y
  ZtZ <- t(Z)%*%Z
  YtZ <- t(Y)%*%Z
  
  Gn <- ((n-dA)/dA) * t(X)%*%P_A%*%X%*%solve(t(X)%*%(diag(n)-P_A)%*%X)
  eigen(Gn)
  log(272)
  
  p <- 0.5
  
  q <- qchisq(1-p,df=5, ncp=0,lower.tail = TRUE,log.p = FALSE)
  Test_Statistic(K_class(1,A,Z,Y,n),A,Z,Y,n)
 
  Test_Statistic(solve(t(Z)%*%Z)%*%t(Z)%*%Y,A,Z,Y,n) 
  
PULSE(A=A,A_1=A_1,X=X,Y=Y,p=0.05,N=1000000,n=n)

solve(t(Z)%*%Z)%*%t(Z)%*%Y
FULLER_k(1,A,A_1,X,Y,n,dA) %>% K_class(.,A,Z,Y,n)
FULLER_k(4,A,A_1,X,Y,n,dA) %>% K_class(.,A,Z,Y,n)


n <- 50
dA <- 5
dX <- 1
dZ <- 1
Rsq <- 0.1
rho <- 0.3

Conc = n*Rsq/(1-Rsq)
Fstat = Conc/dA+1
Fstat

xiEntry <- sqrt(dA*(1-Rsq)*Rsq)/(dA*(1-Rsq))

xi <- matrix(c(rep(xiEntry,dA)),ncol=1)

  U <- mvrnorm(n,mu=rep(0,2),Sigma=matrix(c(1,rho,rho,1),ncol=2,byrow=TRUE))
  A <- mvrnorm(n,mu=rep(0,dA),Sigma=diag(dA))
  X <- A%*%xi+  U[,1]
  Y <- X*truealpha+ U[,2]
  
  A_1= "none"
  Z <- X
  P_A <- A%*%solve(t(A)%*%A)%*%t(A)
  Gn <- ((n-dA)/dA) * t(X)%*%P_A%*%X%*%solve(t(X)%*%(diag(n)-P_A)%*%X)
  eigen(Gn)

  Test_Statistic(K_class(1,A,Z,Y,n),A,Z,Y,n)
  Test_Statistic(K_class(0,A,Z,Y,n),A,Z,Y,n)
  
  PULSE(A=A,A_1=A_1,X=X,Y=Y,p=0.2,N=1000000,n=n)
  solve(t(Z)%*%Z)%*%t(Z)%*%Y
  FULLER_k(1,A,A_1,X,Y,n,dA) %>% K_class(.,A,Z,Y,n)
  FULLER_k(4,A,A_1,X,Y,n,dA) %>% K_class(.,A,Z,Y,n)
  