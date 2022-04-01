#' K-class estimator (kappa parametrization)
#'
#' @description
#' `Kclass` computes the K-class estimate for given kappa parameter k.
#'
#' @details
#' Contains several estimators. k=0: OLS, k=1: TSLS, k=Fuller_k(alpha) Fuller(alpha), k=LIML_k: LIML.
#'
#' @examples
#' ### Just-identified example ###
#' # Generate data
#' n <- 500
#' A <- rnorm(n)
#' H <- rnorm(n)
#' X <- 0.8*A + H + rnorm(n)
#' Y <- 0.5*X + H + rnorm(n)
#' # Compute the OLS estimate
#' Kclass(k = 0, A = A, Z = X, Y = Y)
#' # Compute the TSLS estimate
#' Kclass(k = 1, A = A, Z = X, Y = Y)
#' # Compute the LIML estimate
#' Kclass(k = LIML_k(A,X,Y), A = A, Z = X, Y = Y)
#' # Compute the Fuller(4) estimate
#' Kclass(k = Fuller_k(4,A,X,Y), A = A, Z = X, Y = Y)
#'
#' ### Over-identified example ###
#' # Generate data
#' n  <- 500
#' A1 <- rnorm(n)
#' A2 <- rnorm(n)
#' A3 <- rnorm(n)
#' A  <- cbind(A1,A2,A3)
#' H  <- rnorm(n)
#' X  <- 0.8*A1 + 0.5*A2 + H + rnorm(n)
#' Y  <- 0.5*X + 0.2*A3 + H + rnorm(n)
#' # Compute the OLS estimate
#' Kclass(k = 0, A = A, Z = cbind(A3,X), Y = Y)
#' # Compute the TSLS estimate
#' Kclass(k = 1, A = A, Z = cbind(A3,X), Y = Y)
#' # Compute the LIML estimate
#' Kclass(k = LIML_k(A, X, Y, A_inc = A3), A = A, Z = cbind(A3,X), Y = Y)
#' # Compute the Fuller(4) estimate
#' Kclass(k = Fuller_k(4, A, X, Y, A_inc = A3), A = A, Z = cbind(A3,X), Y = Y)
#'
#' @param k K-class kappa (numeric)
#' @param A exogenous variables (vector/matrix/data.frame)
#' @param Z included exogenous and endogenous regressors (vector/matrix/data.frame)
#' @param Y target variable (vector/matrix/data.frame)
#' @return the K-class estimator with kappa parameter k
#' @export
Kclass <- function(k,A,Z,Y){
  A = as.matrix(A)
  Z = as.matrix(Z)
  Y = as.matrix(Y)
  n <- nrow(Y)
  if(is.null(colnames(Z))){colnames(Z) <- paste0("Z",seq(1,ncol(Z),1))}

  coefs <- solve(
    t(Z) %*% ((1-k)*diag(n)+k* A%*%solve(t(A)%*%A)%*%t(A) )  %*%Z
    ) %*% t(Z) %*% ( (1-k)*diag(n) + k* A%*%solve(t(A)%*%A)%*%t(A) ) %*% Y

  return(coefs)
}


#' K-class estimator (lambda penalty parametrization)
#' @param lambda K-class penalty parameter (numeric)
#' @param A exogenous variables (vector/matrix/data.frame)
#' @param Z included exogenous and endogenous regressors (vector/matrix/data.frame)
#' @param Y target variable (vector/matrix/data.frame)
#' @return the K-class estimator with penalty parameter lambda
#' @export
Kclass_lambda <- function(lambda,A,Z,Y){
  A = as.matrix(A)
  Z = as.matrix(Z)
  Y = as.matrix(Y)
  n <- nrow(Y)
  if(is.null(colnames(Z))){colnames(Z) <- paste0("Z",seq(1,ncol(Z),1))}

  coefs <- solve(
    t(Z) %*% ( diag(n)+ as.vector(lambda)* A%*%solve(t(A)%*%A)%*%t(A) )  %*%Z
  ) %*% t(Z) %*% ( diag(n) + as.vector(lambda)* A%*%solve(t(A)%*%A)%*%t(A) ) %*% Y

  return(coefs)
}


#' Limited information maximum likelihood (LIML) K-class kappa
#'
#' @description
#' `LIML_k` compute the LIML K-class kappa parameter.
#'
#' @examples
#' ### Just-identified example ###
#' # Generate data
#' n <- 500
#' A <- rnorm(n)
#' H <- rnorm(n)
#' X <- 0.8*A + H + rnorm(n)
#' Y <- 0.5*X + H + rnorm(n)
#' # Compute the LIML kappa parameter
#' LIML_k(A, X, Y)
#' # Compute the LIML estimate
#' Kclass(k = LIML_k(A,X,Y), A = A, Z = X, Y = Y)
#'
#' ### Over-identified example ###
#' # Generate data
#' n  <- 500
#' A1 <- rnorm(n)
#' A2 <- rnorm(n)
#' A3 <- rnorm(n)
#' A  <- cbind(A1,A2,A3)
#' H  <- rnorm(n)
#' X  <- 0.8*A1 + 0.5*A2 + H + rnorm(n)
#' Y  <- 0.5*X + 0.2*A3 + H + rnorm(n)
#' # Compute the LIML kappa parameter
#' LIML_k(A, X, Y, A_inc = A3)
#' # Compute the LIML estimate
#' Kclass(k = LIML_k(A, X, Y, A_inc = A3), A = A, Z = cbind(A3,X), Y = Y)
#'
#' @param A exogenous variables (vector/matrix/data.frame)
#' @param A_inc included exogenous regressors (vector/matrix/data.frame)
#' @param X included endogenous regressors (vector/matrix/data.frame)
#' @param Y target variable (vector/matrix/data.frame)
#' @return the K-class kappa parameter corresponding to the LIML estimator
#' @export
LIML_k <- function(A,X,Y,A_inc = NULL){
  A = as.matrix(A)
  if(!is.null(A_inc)){A_inc = as.matrix(A_inc)}
  X = as.matrix(X)
  Y = as.matrix(Y)
  n <- nrow(Y)

  W = t(cbind(Y,X))%*% ( diag(n) - A%*%solve(t(A)%*%A)%*%t(A) ) %*%cbind(Y,X)
  if(is.null(A_inc)){
    W_1 = t(cbind(Y,X))%*%diag(n)%*%cbind(Y,X)
  } else{
    W_1 = t(cbind(Y,X))%*%( diag(n)- A_inc%*%solve(t(A_inc)%*%(A_inc))%*%t(A_inc) )%*%cbind(Y,X)
  }
  LIML_kappa = min(abs(eigen(W_1%*%solve(W))$values))

  return(LIML_kappa)
}

#' Fuller K-class kappa
#'
#' @description
#' `Fuller_k` compute the Fuller K-class kappa parameter for a given fuller hyperparameter, i.e., Fuller(alpha) estimate = Kclass(Fuller_k(alpha,...),...).
#'
#' @details
#' The Fuller(alpha=1) estimator is approximately unbiased and Fuller(alpha=4) approximately yields the minimal MSE estimator among all Fuller estimators.
#'
#' @param alpha The Fuller hyperparameter  (numeric)
#' @param A exogenous variables (vector/matrix/data.frame)
#' @param A_inc included exogenous regressors (vector/matrix/data.frame)
#' @param X included endogenous regressors (vector/matrix/data.frame)
#' @param Y target variable (vector/matrix/data.frame)
#' @return the Fuller estimator with hyperparameter alpha
#' @export
Fuller_k <- function(alpha,A,X,Y,A_inc=NULL){
  A = as.matrix(A)
  if(!is.null(A_inc)){A_inc = as.matrix(A_inc)}
  X = as.matrix(X)
  Y = as.matrix(Y)
  n <- nrow(Y)
  dA <- ncol(A)

  if(is.null(A_inc)){
    LIML_kappa <- LIML_k(A,X,Y,A_inc=NULL)
  } else {
    LIML_kappa <- LIML_k(A,X,Y,A_inc)
  }
  Fuller_kappa = LIML_kappa-alpha/(n-dA)

  return(Fuller_kappa)
}

#' Test-statistic evaluation
#' @param alpha coefficient/parameter to evaluate the test-statistic in (vector)
#' @param A exogenous variables (vector/matrix/data.frame)
#' @param Z included exogenous and endogenous regressors (vector/matrix/data.frame)
#' @param Y target variable (vector/matrix/data.frame)
#' @return the test-statistic evaluated in the coefficient input alpha
#' @export
Test_Statistic <- function(alpha,A,Z,Y){
  A = as.matrix(A)
  Z = as.matrix(Z)
  Y = as.matrix(Y)
  n <- nrow(Y)

  teststat_value <- n*norm(A%*%solve(t(A)%*%A)%*%t(A) %*%(Y-Z%*%alpha),type="F")^2/norm(Y-Z%*%alpha,type="F")^2
  return(teststat_value)
}

#' The P-Uncorrelated Least Squares Estimator (PULSE)
#'
#' @description
#' `PULSE` computes the PULSE estimate; see https://arxiv.org/abs/2005.03353.
#'
#' @details
#' In the over-identified setup the call outputs a message if either (a) the TSLS esimator is rejected, in which case the PULSE reverts to the Fuller(4) estimate, or (b) the OLS estimate is accepted, in which case PULSE coincides with the OLS estimate.
#'
#' The summary output generated with printsummary=TRUE is given by the following columns: method, dim(X+A_inc) columns with coefficient estimates, K-class kappa corresponding to the estimate, the test-statistic evaluated in each estimate, and the p.value for the test that the regression residuals Y-(X,A_ind)*estimate is uncorrelated with the exogenous variables A.
#'
#' p is the rejection threshold for the test of uncorrelated residuals. 1/N is the binary search precision of the K-class parameter in the lambda search space.
#'
#' @examples
#'
#' ### Under-identified example (from section H.3.1 of the paper) ###
#' n <- 500
#' A <- rnorm(n)
#' H <- rnorm(n)
#' X1 <- 1.2*A+ 2*H + rnorm(n)
#' Y <- 0.1*X1+0.8*H+rnorm(n)
#' X2 <- 0.5*Y + rnorm(n)
#' X <- cbind(X1,X2)
#' # Compute the PULSE estimate
#' PULSE(A = A, X = X, Y = Y ,p = 0.05, N = 1000)
#' # Printing the comparison summary
#' PULSE(A = A, X = X, Y = Y ,p = 0.05, N = 1000, printsummary = TRUE)
#'
#' ### Just-identified example ###
#' # Generate data
#' n <- 500
#' A <- rnorm(n)
#' H <- rnorm(n)
#' X <- 0.8*A + H + rnorm(n)
#' Y <- 0.5*X + H + rnorm(n)
#' # Compute the PULSE estimate
#' PULSE(A = A, X = X, Y = Y ,p = 0.05, N = 1000)
#' # Printing the comparison summary
#' PULSE(A = A, X = X, Y = Y ,p = 0.05, N = 1000, printsummary = TRUE)
#'
#' ### Over-identified example ###
#' # Generate data
#' n  <- 500
#' A1 <- rnorm(n)
#' A2 <- rnorm(n)
#' A3 <- rnorm(n)
#' A  <- cbind(A1,A2,A3)
#' H  <- rnorm(n)
#' X  <- 0.8*A1 + 0.5*A2 + H + rnorm(n)
#' Y  <- 0.5*X + 0.2*A3 + H + rnorm(n)
#' # Compute the PULSE estimate
#' PULSE(A = A, A_inc = A3, X = X, Y = Y ,p = 0.05, N = 1000)
#' # Printing the comparison summary
#' PULSE(A = A, A_inc = A3, X = X, Y = Y ,p = 0.05, N = 1000, printsummary = TRUE)
#'
#' @param A exogenous variables (vector/matrix/data.frame).
#' @param A_inc included exogenous regressors. Default NULL, i.e., no included exogeneous regressors. (vector/matrix/data.frame).
#' @param X included endogenous regressors (vector/matrix/data.frame).
#' @param Y target variable (vector/matrix/data.frame).
#' @param p the p-value PULSE hyperparameter. Default 0.05 (numeric between 0 and 1).
#' @param N reciprocal binary search precision, i.e. a precision of 1/N. Default 1000 (integer/numeric).
#' @param printsummary prints comparison summary. Default FALSE (logical).
#' @return the estimated PULSE estimate.
#' @importFrom stats qchisq pchisq
#' @import tibble
#' @import dplyr
#' @import tidyr
#' @import purrr
#' @importFrom magrittr %>%
#' @export

PULSE <- function(A,X,Y, p = 0.05, N = 1000,A_inc = NULL,printsummary = FALSE){
  A = as.matrix(A)
  if(is.null(colnames(A))){colnames(A) <- paste0("A",seq(1,ncol(A),1))}
  X = as.matrix(X)
  if(is.null(colnames(X))){colnames(X) <- paste0("X",seq(1,ncol(X),1))}
  Y = as.matrix(Y)
  if(!is.null(A_inc)){
    A_inc = as.matrix(A_inc)
    if(is.null(colnames(A_inc))){colnames(A_inc) <- paste0("A_inc",seq(1,ncol(A_inc),1))}
    Z = cbind(A_inc,X)
  } else {
      Z = X
  }
  if(is.null(colnames(Z))){colnames(Z) <- paste0("Z",seq(1,ncol(Z),1))}
  n <- nrow(Y)
  dX <- ncol(X)
  dA <- ncol(A)
  dZ <- ncol(Z)
  q <- stats::qchisq(1-p,df=dA, ncp=0,lower.tail = TRUE,log.p = FALSE)

  UNDERID = FALSE
  if(dX > dA){
    UNDERID = TRUE
  }
  message <- NULL
  stopBool = FALSE
  if(UNDERID == FALSE){
    if(Test_Statistic(Kclass(1,A,Z,Y),A,Z,Y)>= q){
      message <- "Note: TSLS was rejected, reverting to Fuller(4)."
      Fuller4Kappa <- Fuller_k(4,A,X,Y,A_inc)
      coefs <- Kclass(Fuller4Kappa,A,Z,Y)
      stopBool = TRUE
    }
  }
  if (UNDERID == TRUE || stopBool == FALSE){

    if(Test_Statistic(Kclass(0,A,Z,Y),A,Z,Y)<= q){
      message <- "Note: OLS was accepted."
      coefs <- Kclass(0,A,Z,Y)
    } else {
      lmax <- 2
      lmin <- 0

      while(Test_Statistic(Kclass_lambda(lmax,A,Z,Y),A,Z,Y)>q){
        lmin <- lmax
        lmax <- lmax^2
      }

      Delta <- lmax-lmin
      while(Delta > 1/N || Accept == FALSE){
        l <- (lmin+lmax)/2
        alpha <- Kclass_lambda(l,A,Z,Y)
        TestStatistic <- Test_Statistic(alpha,A,Z,Y)
        if(TestStatistic>q){
          Accept <- FALSE
          lmin <- l

        }
        else {
          Accept <- TRUE
          lmax <- l
        }
        Delta <- lmax-lmin
      }
      coefs <- Kclass_lambda(lmax,A,Z,Y)
    }
  }

  if(printsummary == TRUE){
    if(UNDERID == FALSE){
      coefficients <- colnames(Z)
      ols <- Kclass(0,A,Z,Y)
      tsls <- Kclass(1,A,Z,Y)
      pulse <- coefs
      fuller1 <- Kclass(Fuller_k(1,A,X,Y,A_inc),A,Z,Y)
      fuller4 <- Kclass(Fuller_k(4,A,X,Y,A_inc),A,Z,Y)
      liml <- Kclass(LIML_k(A,X,Y,A_inc),A,Z,Y)

      d1 <- data.frame(coefficients,ols, pulse, tsls, liml, fuller1, fuller4) %>%
        tidyr::pivot_longer(cols=-c(coefficients),names_to="method",values_to="coef") %>%
        tidyr::pivot_wider(names_from=coefficients,values_from=coef) %>%
        dplyr::mutate(Kclass_kappa = c(0,lmax/(1+lmax),1,LIML_k(A,X,Y,A_inc),Fuller_k(1,A,X,Y,A_inc),Fuller_k(4,A,X,Y,A_inc)),
                      Test_statistic = c(Test_Statistic(ols,A,Z,Y),
                                         Test_Statistic(pulse,A,Z,Y),
                                         Test_Statistic(tsls,A,Z,Y),
                                         Test_Statistic(liml,A,Z,Y),
                                         Test_Statistic(fuller1,A,Z,Y),
                                         Test_Statistic(fuller4,A,Z,Y)
                                         ),
                      p.value =  stats::pchisq(Test_statistic,df=dA, ncp=0,lower.tail = FALSE,log.p = FALSE)) %>%
        dplyr::mutate(dplyr::across(-c(method) , ~ sprintf("%.10f", .x) ))
      print(d1)
    } else {
      coefficients <- colnames(Z)
      ols <- Kclass(0,A,Z,Y)
      pulse <- coefs
      fuller1 <- Kclass(Fuller_k(1,A,X,Y,A_inc),A,Z,Y)
      fuller4 <- Kclass(Fuller_k(4,A,X,Y,A_inc),A,Z,Y)

      d1 <- data.frame(coefficients,ols, pulse, fuller1, fuller4) %>%
        tidyr::pivot_longer(cols=-c(coefficients),names_to="method",values_to="coef") %>%
        tidyr::pivot_wider(names_from=coefficients,values_from=coef) %>%
        dplyr::mutate(Kclass_kappa = c(0,lmax/(1+lmax),Fuller_k(1,A,X,Y,A_inc),Fuller_k(4,A,X,Y,A_inc)),
                      Test_statistic = c(Test_Statistic(ols,A,Z,Y),
                                         Test_Statistic(pulse,A,Z,Y),
                                         Test_Statistic(fuller1,A,Z,Y),
                                         Test_Statistic(fuller4,A,Z,Y)
                      ),
                      p.value =  stats::pchisq(Test_statistic,df=dA, ncp=0,lower.tail = FALSE,log.p = FALSE)) %>%
        dplyr::mutate(dplyr::across(-c(method) , ~ sprintf("%.10f", .x) ))
      print(d1)
    }
    #Print a summary with estimator comparisons.
  }
  message(message)
  return(coefs)
}
