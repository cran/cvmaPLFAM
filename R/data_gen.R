#' @title Simulated data
#'
#' @description Simulate sample data for illustration, including a \code{M0}-column design matrix of scalar predictors,
#'      a \code{100}-column matrix of the functional predictor, a one-column vector of \code{mu}, a one-column vector of \code{Y},
#'      and a one-column vector of \code{testY}.
#'
#' @param R A scalar of value ranging from \code{0.1} to \code{0.9}. The ratio of \code{var(mu)/var(Y)}.
#' @param K A scalar. The number of replications.
#' @param n A scalar. The sample size of training data.
#' @param ntest A scalar. The sample size of test data.
#' @param M0 A scalar. True dimension of scalar predictors.
#' @param typ A scalar of value \code{1} - \code{2}. Type of the effect for the functional predictor.
#' @param design A scalar of value \code{1} - \code{3}. Correspond to simulation studies.
#'
#' @return A \code{list} of \code{K} simulated training data sets and \code{K} simulated test data sets. Each data set is of \code{matrix} type,
#'      whose first \code{M0} columns corresponds to the design matrix of scalar predictors, followed by the
#'      recording/measurement matrix of the functional predictor, and vectors \code{mu}, \code{Y}.
#'
#'
#' @export
#'
#' @import MASS
#' @examples
#' library(MASS)
#' # Example: Design 1 in simulation study
#' set.seed(22)
#' data1 <- data_gen(R = 0.6, K = 2, n = 10, ntest = 5, M0 = 4, typ = 1, design = 1)
#' str(data1)
#' # List of 4
#' #$ : num [1:10, 1:106] -0.501 -1.266 -0.564 -0.563 -0.395 ...
#' #$ : num [1:10, 1:106] -1.207 -0.089 -0.782 0.123 0.66 ...
#' #$ : num [1:5, 1:106] 0.816 0.679 0.816 -0.563 -1.367 ...
#' #$ : num [1:5, 1:106] -0.089 -0.785 0.899 -0.785 -0.445 ...
#'
#'
#' # Example: Design 2 in simulation study
#' data_gen(R = 0.3, K = 3, n = 10, ntest = 5, M0 = 20, typ = 1, design = 2)
#'
#' # Example: Design 3 in simulation study
#' data_gen(R = 0.9, K = 5, n = 20, ntest = 10, M0 = 4, typ = 2, design = 3)
#'
#'

data_gen <- function(R, K, n, ntest, M0, typ, design)
{
  s0 = 1:M0
  Sigma1 = matrix(rep(1:(M0+1), each = M0+1), M0 + 1, M0 + 1)
  Sigma1 = 0.5^abs(t(Sigma1) - Sigma1) # covariance matrix for linear part X

  datalist = vector('list', K*2)

  if(design == 1){
    # uncorrelated, homo
    K0 = 5
    tt = seq(0, 1, by = 0.01)[-1]
    beta = c(2.1, -0.7, 0.4, 0.9)
    data = matrix(0, n, M0 + length(tt) + 2)
    testdata = matrix(0, ntest, M0 + length(tt) + 2)

    for(g in 1:K){
      X = mvrnorm(n, rep(0, M0), Sigma1[-(M0+1),-(M0+1)])   # n*M0
      Z = Zt_gen(n, K0, xi01 = NULL, tt, v = 0.2)
      gz = gg(Z$xi, typ)
      mu = c(X%*%beta) + gz                # 1*n
      eta = sqrt(var(mu)*(1-R) / R)
      e = rnorm(n) * eta             # homo
      Y = mu + e

      data[,1:M0] = X
      data[,(M0+1):(M0+length(tt))] = Z$Zt
      data[,M0+length(tt)+1] = mu
      data[,M0+length(tt)+2] = Y
      datalist[[g]] = data

      # test data
      index.boots = sample(1:n, ntest, replace = TRUE)
      X.boots = X[index.boots,]            # ntest * M0
      Zt.boots = (Z$Zt)[index.boots,]      # ntest * length(tt) matrix
      gz.boots = gz[index.boots]           # 1*ntest vector
      mu.boots = mu[index.boots]           # 1*ntest
      e1 = rnorm(ntest) * eta               # homo
      Y.test = mu.boots + e1

      testdata[,1:M0] = X.boots
      testdata[,(M0+1):(M0+length(tt))] = Zt.boots
      testdata[,M0+length(tt)+1] = mu.boots
      testdata[,M0+length(tt)+2] = Y.test
      datalist[[K+g]] = testdata

    } # end g

  }else if(design == 2){
    # correlated, heuristic
    K0 = 10
    tt = seq(0, 1, by = 0.01)[-1]
    beta = s0^(-3/2)
    data = matrix(0, n, M0 + length(tt) + 2)
    testdata = matrix(0, ntest, M0 + length(tt) + 2)

    for(g in 1:K){
      Xxi = mvrnorm(n, rep(0, M0+1), Sigma1)   # n*M0
      X = Xxi[,1:M0]
      Z = Zt_gen(n, K0, Xxi[,M0+1], tt, v = 0.2)
      gz = gg(Z$xi, typ)
      mu = c(X%*%beta) + gz                      # 1*n
      eta = sqrt(var(mu)*(1-R)/(R*var(X[,2])))   # 1*1
      e = rnorm(n)*eta*sqrt(X[,2]^2 + 0.01)
      Y = mu + e

      data[,1:M0] = X
      data[,(M0+1):(M0+length(tt))] = Z$Zt
      data[,M0+length(tt)+1] = mu
      data[,M0+length(tt)+2] = Y
      datalist[[g]] = data

      # test data
      index.boots = sample(1:n, ntest, replace = T)
      X.boots = X[index.boots,]
      Zt.boots = (Z$Zt)[index.boots,]
      gz.boots = gz[index.boots]
      mu.boots = mu[index.boots]                      # 1*n
      e1 = rnorm(ntest) * eta * sqrt(X.boots[,2]^2 + 0.01)
      Y.test = mu.boots + e1

      testdata[,1:M0] = X.boots
      testdata[,(M0+1):(M0+length(tt))] = Zt.boots
      testdata[,M0+length(tt)+1] = mu.boots
      testdata[,M0+length(tt)+2] = Y.test
      datalist[[K+g]] = testdata

    } # end g


  }else if(design == 3){
    ## linear model
    # uncorrelated, homo
    K0 = 5
    tt = seq(0, 1, by = 0.01)[-1]
    beta = c(2.1, -0.7, 0.4, 0.9)
    data = matrix(0, n, M0 + length(tt) + 2)
    testdata = matrix(0, ntest, M0 + length(tt) + 2)

    for(g in 1:K){
      X = mvrnorm(n, rep(0, M0), Sigma1[-(M0+1),-(M0+1)])   # n*M0
      Z = Zt_gen(n, K0, xi01 = NULL, tt, v = 0.2)
      gz = gg(Z$Zt, typ)
      mu = c(X%*%beta) + gz                # 1*n
      eta = sqrt(var(mu)*(1-R) / R)
      e = rnorm(n) * eta             # homo
      Y = mu + e

      data[,1:M0] = X
      data[,(M0+1):(M0+length(tt))] = Z$Zt
      data[,M0+length(tt)+1] = mu
      data[,M0+length(tt)+2] = Y
      datalist[[g]] = data

      # test data
      index.boots = sample(1:n, ntest, replace = TRUE)
      X.boots = X[index.boots,]            # ntest * M0
      Zt.boots = (Z$Zt)[index.boots,]      # ntest * length(tt) matrix
      gz.boots = gz[index.boots]           # 1*ntest vector
      mu.boots = mu[index.boots]           # 1*ntest
      e1 = rnorm(ntest) * eta               # homo
      Y.test = mu.boots + e1

      testdata[,1:M0] = X.boots
      testdata[,(M0+1):(M0+length(tt))] = Zt.boots
      testdata[,M0+length(tt)+1] = mu.boots
      testdata[,M0+length(tt)+2] = Y.test
      datalist[[K+g]] = testdata

    } # end g

  }

  return(datalist)

}

