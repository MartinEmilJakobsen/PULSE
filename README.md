# Distributional robustness of K-class estimators and the PULSE
Implementation of the PULSE estimator and replication package repository for "Distributional robustness of K-class estimators and the PULSE" @ https://doi.org/10.1093/ectj/utab031 + https://arxiv.org/abs/2005.03353 . 


## The PULSE R-Package


### Installation of R-package

```R
devtools::install_github("MartinEmilJakobsen/PULSE", subdir = "PULSE")
```

### Usage

```R
PULSE(A, X, Y, p = 0.05, N = 1000, A_inc = NULL, printsummary = FALSE)
```

##### Arguments
* `A`: exogenous variables (vector/matrix/data.frame row-wise obeservations).
* `X`: included endogenous regressors (vector/matrix/data.frame).
* `Y`: target variable (vector/matrix/data.frame).
* `p`: the p_min PULSE hyperparameter. Default 0.05 (numeric between 0 and 1).
* `N`: reciprocal binary search precision, i.e. a precision of 1/N. Default 1000 (integer/numeric).
* `A_inc`: included exogenous regressors. Default NULL, i.e., no included exogeneous regressors. (vector/matrix/data.frame).
* `printsummary`: prints comparison summary. Default FALSE (logical).

##### Description
PULSE computes the PULSE estimate; see https://arxiv.org/abs/2005.03353.

##### Details
In the over-identified setup the call outputs a message if either (a) the TSLS esimator is rejected, in which case the PULSE reverts to the Fuller(4) estimate, or (b) the OLS estimate is accepted, in which case PULSE coincides with the OLS estimate. 

The summary output generated with printsummary=TRUE is given by the following columns: method, dim(X+A_inc) columns with coefficient estimates, K-class kappa corresponding to the estimate, the test-statistic evaluated in each estimate, and the p.value for the test that the regression residuals Y-(X,A_ind)*estimate is uncorrelated with the exogenous variables A.

p is the rejection threshold for the test of uncorrelated residuals. 1/N is the binary search precision of the K-class parameter in the lambda search space.

##### Examples

```R
### Under-identified example (from section H.3.1 of the paper) ###
n <- 500
A <- rnorm(n)
H <- rnorm(n)
X1 <- 1.2*A+ 2*H + rnorm(n)
Y <- 0.1*X1+0.8*H+rnorm(n)
X2 <- 0.5*Y + rnorm(n)
X <- cbind(X1,X2)
# Compute the PULSE estimate
PULSE(A = A, X = X, Y = Y ,p = 0.05, N = 1000)
# Printing the comparison summary
PULSE(A = A, X = X, Y = Y ,p = 0.05, N = 1000, printsummary = TRUE)
```

```R
### Just-identified example ###
# Generate data
n <- 500
A <- rnorm(n)
H <- rnorm(n)
X <- 0.8*A + H + rnorm(n)
Y <- 0.5*X + H + rnorm(n)
# Compute the PULSE estimate
PULSE(A = A, X = X, Y = Y ,p = 0.05, N = 1000)
# Printing the comparison summary
PULSE(A = A, X = X, Y = Y ,p = 0.05, N = 1000, printsummary = TRUE)
```

```R
### Over-identified example ###
# Generate data
n  <- 500
A1 <- rnorm(n)
A2 <- rnorm(n)
A3 <- rnorm(n)
A  <- cbind(A1,A2,A3)
H  <- rnorm(n)
X  <- 0.8*A1 + 0.5*A2 + H + rnorm(n)
Y  <- 0.5*X + 0.2*A3 + H + rnorm(n)
# Compute the PULSE estimate
PULSE(A = A, A_inc = A3, X = X, Y = Y ,p = 0.05, N = 1000)
# Printing the comparison summary
PULSE(A = A, A_inc = A3, X = X, Y = Y ,p = 0.05, N = 1000, printsummary = TRUE)
```

## Replication package
The replication package is located in the subdirectory "/Replication Package"

### Replicate Empirical Analysis
See readme.txt in sub-folder /Empirical_Analysis

### Replicate experiments

This repository also contains all R scripts for conducting the experiments referenced in the paper. In /Data the experiment data used in the paper is present. Thus, to simply rerun analysis and illustration generation on this data start from step 2) below and ignore comments about data-pointers.

If you want to replicate all experiments then start from step 1) and do not ignore comments about data-pointer as the subsequent analysis and illustration generation has to point to the newly generated data. Note that the code is parallelized and depending on the system you need to specify the number of availiable cores in "plan(multiprocess, workers = [NO.CORES])" in each script. On a 64 core system (Intel Xeon E7-4860 v2 @ 2.6GHz and 256Gb memory) each of the VaryingConfounding scripts has a runtime of 3 hours (2.5 billion draws from a structural model), while the FixedConfounding script has a runtime of 9 hours (7.5 billion draws).

#### 1) Data generation 

##### 1.a) Initial Experiments

To generate the experiment data: Run Script.R or manually run the following scripts. Each script saves the data at "Data/[SCRIPTNAME]_ [REPLICATIONS]_ [OBSERVATIONS]_ [EXECUTIONTIME].RDS" and saves a copy of the script at "Data/Log/[SCRIPTNAME]_ [REPLICATIONS]_ [OBSERVATIONS]_ [EXECUTIONTIME].RDS".

1. **Experiment_Univariate.R**
    * *This script runs the univariate experiment*
2. **Experiment_Multivariate_VaryingConfounding_Beta00_PULSE05.R**
    * *This scripts runs the multivariate experiment with varying confounding for $\beta=(0,0)$. It computes the performance measures for OLS, Fuller(1), Fuller(4) and PULSE with $p_{\min}=0.05$*
3. **Experiment_Multivariate_VaryingConfounding_Beta00_PULSE10.R**
    * *This scripts runs the multivariate experiment with varying confounding for $\beta=(0,0)$. It computes the performance measures for OLS, Fuller(1), Fuller(4) and PULSE with $p_{\min}=0.1$*
4. **Experiment_Multivariate_VaryingConfounding_Beta11_PULSE05.R**
    * *This scripts runs the multivariate experiment with varying confounding for $\beta=(1,1)$. It computes the performance measures for OLS, Fuller(1), Fuller(4) and PULSE with $p_{\min}=0.05$*
5. **Experiment_Multivariate_VaryingConfounding_Beta-11_PULSE05.R**
    * *This scripts runs the multivariate experiment with varying confounding for $\beta=(-1,1)$. It computes the performance measures for OLS, Fuller(1), Fuller(4) and PULSE with $p_{\min}=0.05$*
6. **Experiment_Multivariate_FixedConfounding.R**
    * *This scripts runs the multivariate experiment with fixed confounding for $\beta=(0,0)$. It computes the performance measures for OLS, Fuller(1), Fuller(4) and PULSE with $p_{\min}=0.05$*
7. **Experiment_Univariate_AlternativeSetups.R**
    * *This scripts runs the univariate experiment under the alternative setups.
7. **Experiment_UnderIdentified.R**
    * *This scripts runs the under-identified experiment.
    
##### 1.b) Run Experiment for Superior Models
This will extract the coefficients for the models where PULSE(05) is MSE superior to Fuller(4) of the "Experiment_Multivariate_VaryingConfounding_Beta00_PULSE05.R" experiment and rerun the experiment for these models with 25000 repetitions to account for selection bias. 

1. Locate the datafile "Data/Experiment_Multivariate_VaryingConfounding_Beta00_PULSE05_[REPLICATIONS]_ [OBSERVATIONS]_ [EXECUTIONTIME].RDS" of 1.a). 
2. Edit the data pointer in **Analysis_SelectingSuperiorModels.R** to this datafile and run the script. 
    * This script saves "Data/VaryingConfounding_MSESuperiorModelData_[EXECUTIONTIME].RDS" with the extracted models. 
3. Edit the datapointer in "Experiment_Multivariate_VaryingConfounding_SuperiorModels.R" to the extracted models and run the script.


#### 2) Data Analysis 

All analysis on experiment data reported in the paper can be found in

1. **Analysis_Multivariate_VaryingConfounding.R**
2. **Analysis_Multivariate_FixedConfounding.R**

Note that data-pointers in these scripts have to be adjusted to if used on new experiment data.

#### 3) Generate Illustrations 

All illustrations in the paper is generated by the scripts:

1. **Plots_Univariate.R**
2. **Plots_Multivariate.R**
3. **LevelsetPlots.R**
4. **DistributionalRobustness.R**


Note that data-pointers in script 1 & 2 have to be adjusted to if used on new experiment data.


####  Version Control and Packages
VersionControl_readme.csv contains R-version ID and packages used and all their dependencies.
