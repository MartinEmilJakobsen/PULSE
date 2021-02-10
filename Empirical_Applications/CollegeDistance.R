library(AER)

# load CollegeDistance data set
data("CollegeDistance")

Selected_Data <- CollegeDistance

Selected_Data %>%  select(wage) %>%  unique()

n <- nrow(Selected_Data)
#Target
Y <- Selected_Data %>% select(education)  %>% as.matrix
#Included Exogenous:
A_1 <- Selected_Data %>% select(gender,ethnicity,score,fcollege,mcollege,home,urban,unemp,wage,tuition, income,region) %>%  as.matrix
#All Exogenous
A <- Selected_Data %>% select(distance,gender,ethnicity,score,fcollege,mcollege,home,urban,unemp,wage,tuition, income,regio) %>% as.matrix
#Included Endogenous
X <- Selected_Data %>% select(ed76,exp,expsq) %>%  as.matrix
#Included (all)
Z <- cbind(X,A_1)

#OLS
K_class(0,A,Z,Y,n)

#IV
K_class(1,A,Z,Y,n)

#PULSE
PULSE(A,A_1,X,Y,p=0.05,N=10000,n)