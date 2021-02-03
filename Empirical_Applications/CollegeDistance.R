library(AER)

# load CollegeDistance data set
data("CollegeDistance")
# read out relevant variables
Y <- CollegeDistance$score
X <- CollegeDistance$education
I <- CollegeDistance$distance

fit <- lm(X ~ I)
summary(fit)

##
# Ordinary least squares
##

fit <- lm(Y ~ X)
resid <- residuals(fit)
summary(fit)

# correlation between X and residuals
var(resid)
cor.test(X,resid)

##
# 2-stage least squares
##

# Step 1: regress predictor X on instrument I
step1.fit <- lm(X ~ I)
X.pred <- fitted.values(step1.fit)

# Step 2: regress predicted values of step 1 on Y
step2.fit <- lm(Y ~ X.pred)
resid <- residuals(step2.fit)
summary(step2.fit)

# correlation between X and residuals
var(resid)
cor.test(X,resid)

beta = cov(Y,I)/cov(X,I)
print(beta)
