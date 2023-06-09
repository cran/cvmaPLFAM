#' @title Simulated data
#'
#' @description Simulate sample data for illustration, including a \code{M0}-column design matrix of scalar predictors,
#'      a \code{100}-column matrix of the functional predictor, a one-column vector of \code{mu}, a one-column vector of \code{Y},
#'      and a one-column vector of \code{testY}.
#'
#' @param R A scalar of value ranging from \code{0.1} to \code{0.9}. The ratio of \code{var(mu)/var(Y)}.
#' @param K A scalar. The number of replications.
#' @param n A scalar. The sample size of simulated data.
#' @param M0 A scalar. True dimension of scalar predictors.
#' @param typ A scalar of value \code{1} or \code{2}. Type of the additive function for the functional predictor.
#' @param design A scalar of value \code{1}, \code{2}, or \code{3}. Designs 1, 2, 3 corresponding to simulation studies.
#'
#' @return A \code{list} of \code{K} simulated data sets. Each data set is of \code{matrix} type,
#'      whose first \code{M0} columns corresponds to the design matrix of scalar predictors, followed by the
#'      recording/measurement matrix of the functional predictor, and vectors \code{mu}, \code{Y}, \code{testY}.
#'
#'
#' @export
#'
#' @import MASS
#' @examples
#' library(MASS)
#' # Example: Design 1 in simulation study
#' data_gen(R = 0.6, K = 2, n = 10, typ = 1, design = 1)
#'
#' # Example: Design 2 in simulation study
#' data_gen(R = 0.3, K = 3, n = 10, typ = 2, design = 2)
#'
#' # Example: Design 3 in simulation study
#' data_gen(R = 0.9, K = 5, n = 20, typ = 1, design = 3)
#'
#'

data_gen <- function(R, K, n, M0 = 50, typ, design)
{
  s0 = 1:M0
  Sigma1 = matrix(0, M0 + 1, M0 + 1)
  for(i in 1:(M0+1)){
    for(j in 1:(M0+1)){
      Sigma1[i,j] = 0.5^abs(i-j)
    }
  } # covariance matrix for linear part X

  datalist = vector('list', K)

  if(design == 1){

    # uncorrelated, homo
    K0 = 50
    tt = seq(0, 1, by = 0.01)[-1]
    beta = s0^(-3/2)
    data = matrix(0, n, M0 + length(tt) + 3)

    for(g in 1:K){
      X = mvrnorm(n, rep(0, M0), Sigma1[-(M0+1),-(M0+1)])   # n*M0
      Z = Zt_gen(n, K0, tt, xi01 = NULL, v = 0.2)
      gz = gg(Z$xi, typ)                   # typ
      mu = c(X%*%beta) + gz                # 1*n
      eta = sqrt(var(mu)*(1-R) / R)
      e = rnorm(n) * eta             # homo
      Y = mu + e
      testY = mu + rnorm(n) * eta
      data[,1:M0] = X
      data[,(M0+1):(M0+length(tt))] = Z$Zt
      data[,M0+length(tt)+1] = mu
      data[,M0+length(tt)+2] = Y
      data[,M0+length(tt)+3] = testY
      datalist[[g]] = data
    } # end g

  }else if(design == 2){

    # correlated, heuristic
    K0 = 20
    tt = seq(0, 10, by = 0.1)[-1]
    beta = s0^(-1/2)
    data = matrix(0, n, M0 + length(tt) + 3)

    for(g in 1:K){
      Xxi = mvrnorm(n, rep(0, M0+1), Sigma1)   # n*M0+1
      X = Xxi[,1:M0]
      Z = Zt_gen(n, K0, Xxi[,M0+1], tt, v = 0.2)
      gz = gg(Z$xi, typ)                   # typ = 1 or 2
      mu = c(X%*%beta) + gz                   # 1*n
      U = 2*runif(n) - 1
      eta = sqrt(var(mu)*(1-R) / (R*var(U)))
      e = rnorm(n) * eta * sqrt(U^2 + 0.01)   # heuristic
      Y = mu + e
      U1 = 2*runif(n) - 1
      eta1 = sqrt(var(mu)*(1-R) / (R*var(U1)))
      testY = mu + rnorm(n) * eta1 * sqrt(U1^2 + 0.01)
      data[,1:M0] = X
      data[,(M0+1):(M0+length(tt))] = Z$Zt
      data[,M0+length(tt)+1] = mu
      data[,M0+length(tt)+2] = Y
      data[,M0+length(tt)+3] = testY
      datalist[[g]] = data
    } # end g

  }else{

    # correlated, heuristic
    K0 = 50
    tt = seq(0, 1, by = 0.01)[-1]
    beta = 1/s0
    data = matrix(0, n, M0 + length(tt) + 3)

    for(g in 1:K){
      Xxi = mvrnorm(n, rep(0, M0+1), Sigma1)   # n*M0
      X = Xxi[,1:M0]
      Z = Zt_gen(n, K0, Xxi[,M0+1], tt, v = 0.2)
      gz = gg(Z$xi, typ)                   # typ = 1 or 2 or 3
      mu = c(X%*%beta) + gz                      # 1*n
      eta = sqrt(var(mu)*(1-R)/(R*var(X[,2])))   # 1*1
      e = rnorm(n)*eta*sqrt(X[,2]^2 + 0.01)
      Y = mu + e
      testY = mu + rnorm(n)*eta*sqrt(X[,2]^2 + 0.01)
      data[,1:M0] = X
      data[,(M0+1):(M0+length(tt))] = Z$Zt
      data[,M0+length(tt)+1] = mu
      data[,M0+length(tt)+2] = Y
      data[,M0+length(tt)+3] = testY
      datalist[[g]] = data
    } # end g

  }

  return(datalist)

}

