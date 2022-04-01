# Distributional Robustness of K-class Estimators and the PULSE
Repository for simulations and illustrations for "Distributional Robustness of K-class Estimators and the PULSE" @ https://arxiv.org/abs/2005.03353. 


## Install PULSE R-Package
```R
devtools::install_github("MartinEmilJakobsen/PULSE", subdir = "PULSE")
```



## Implementation of K-class estimators and the PULSE
Implementation of Estimators is found in "Estimators.R". The PULSE estimate is computed by running:

PULSE(A,A_1,X,Y,p,N)

where A,A_1,X,Y are the datamatrices (rowwise obeservations) for: all exogenous, included exogenous, included endogenous and target endogeous, respectively. If there is no included exogenous variables, then set A_1= "none". p is the $p_{min}$ for PULSE rejection threshold. 1/N is the precision of the binary search (in the lambda space).

## Replicate Empirical Analysis
See readme.txt in sub-folder /Empirical_Analysis

## Replicate experiments

This repository also contains all R scripts for conducting the experiments referenced in the paper. In /Data the experiment data used in the paper is present. Thus, to simply rerun analysis and illustration generation on this data start from step 2) below and ignore comments about data-pointers.

If you want to replicate all experiments then start from step 1) and do not ignore comments about data-pointer as the subsequent analysis and illustration generation has to point to the newly generated data. Note that the code is parallelized and depending on the system you need to specify the number of availiable cores in "plan(multiprocess, workers = [NO.CORES])" in each script. On a 64 core system (Intel Xeon E7-4860 v2 @ 2.6GHz and 256Gb memory) each of the VaryingConfounding scripts has a runtime of 3 hours (2.5 billion draws from a structural model), while the FixedConfounding script has a runtime of 9 hours (7.5 billion draws).

### 1) Data generation 

#### 1.a) Initial Experiments

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
    
#### 1.b) Run Experiment for Superior Models
This will extract the coefficients for the models where PULSE(05) is MSE superior to Fuller(4) of the "Experiment_Multivariate_VaryingConfounding_Beta00_PULSE05.R" experiment and rerun the experiment for these models with 25000 repetitions to account for selection bias. 

1. Locate the datafile "Data/Experiment_Multivariate_VaryingConfounding_Beta00_PULSE05_[REPLICATIONS]_ [OBSERVATIONS]_ [EXECUTIONTIME].RDS" of 1.a). 
2. Edit the data pointer in **Analysis_SelectingSuperiorModels.R** to this datafile and run the script. 
    * This script saves "Data/VaryingConfounding_MSESuperiorModelData_[EXECUTIONTIME].RDS" with the extracted models. 
3. Edit the datapointer in "Experiment_Multivariate_VaryingConfounding_SuperiorModels.R" to the extracted models and run the script.


### 2) Data Analysis 

All analysis on experiment data reported in the paper can be found in

1. **Analysis_Multivariate_VaryingConfounding.R**
2. **Analysis_Multivariate_FixedConfounding.R**

Note that data-pointers in these scripts have to be adjusted to if used on new experiment data.

### 3) Generate Illustrations 

All illustrations in the paper is generated by the scripts:

1. **Plots_Univariate.R**
2. **Plots_Multivariate.R**
3. **LevelsetPlots.R**
4. **DistributionalRobustness.R**


Note that data-pointers in script 1 & 2 have to be adjusted to if used on new experiment data.


##  Version Control and Packages
See the sub-folder /Version_Control_And_Packages: readme.txt contains R-version ID and packages used, Packages_and_dependencies.rar contains a copy of all packages and their dependencies.
