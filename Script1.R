library(MASS) # mvrnorm
library(matrixcalc)
library(stats)
library(lqmm) #makepositive.
library(expm) #sqrtm A=SS
library(furrr)
library(tidyverse)
library(metR)
library(gridExtra)
library(rlang)
library(stringr)

message("EXECUTING SCRIPT1")

message("EXECUTING Running_Experiment_Univariate.R")
message("REMAINING: Experiment_Multivariate_VaryingConfounding_Beta00_PULSE05.R , Experiment_Multivariate_VaryingConfounding_Beta00_PULSE10.R , Experiment_Multivariate_VaryingConfounding_Beta11_PULSE05.R")
source(file="Experiment_Univariate.R")

message("EXECUTING Experiment_Multivariate_VaryingConfounding_Beta00_PULSE05.R")
source(file="Experiment_Multivariate_VaryingConfounding_Beta00_PULSE05.R")

message("EXECUTING Experiment_Multivariate_VaryingConfounding_Beta00_PULSE10.R")
source(file="Experiment_Multivariate_VaryingConfounding_Beta00_PULSE10.R")

message("EXECUTING Experiment_Multivariate_VaryingConfounding_Beta11_PULSE05.R")
source(file="Experiment_Multivariate_VaryingConfounding_Beta11_PULSE05.R")
