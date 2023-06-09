#' @title Output the prediction risks of each method for partial linear functional additive models (PLFAMs)
#' @description Calculate the estimated weights for averaging across all candidate models and the corresponding
#'      mean squared prediction error risk. The methods include AIC, BIC, SAIC, SBIC, and CVMA for PLFAMs.
#'
#' @param M The number of candidate models.
#' @param nump The number of scalar predictors in candidate models.
#' @param numq The number of funtional principal components (FPCs) in candidate models.
#' @param a2 The number of FPCs in each candidate model. See \code{\link{modelspec}}.
#' @param a3 The index for each component in each candidate model. See \code{\link{modelspec}}.
#' @param nfolds The number of folds used in cross-validation.
#' @param XX.train The training data of predictors processed.
#' @param Y.train The training data of response variable.
#' @param XX.pred The test data of predictors processed.
#' @param Y.pred The test data of response variable.
#'
#' @return A \code{list} of
#'     \item{aic}{Mean squared error risk in training data set, produced by AIC model selection method.}
#'     \item{bic}{Mean squared error risk in training data set, produced by BIC model selection method.}
#'     \item{saic}{Mean squared error risk in training data set, produced by SAIC model averaging method.}
#'     \item{sbic}{Mean squared error risk in training data set, produced by SBIC model averaging method.}
#'     \item{cv}{Mean squared error risk in training data set, produced by CVMA method.}
#'     \item{ws}{A \code{list} of weights estimator for five methods.}
#'     \item{predaic}{Mean squared prediction error risk in test data set, produced by AIC model selection method.}
#'     \item{predbic}{Mean squared prediction error risk in test data set, produced by BIC model selection method.}
#'     \item{predsaic}{Mean squared prediction error risk in test data set, produced by SAIC model averaging method.}
#'     \item{predsbic}{Mean squared prediction error risk in test data set, produced by SBIC model averaging method.}
#'     \item{predcv}{Mean squared prediction error risk in test data set, produced by CVMA method.}
#'
#' @export
#'
#' @import fda
#' @import quadprog
#' @import mgcv
#' @import stats
#' @importFrom utils combn
#'
#'
#'

predRisk <- function(M, nump, numq, a2, a3, nfolds, XX.train, Y.train, XX.pred, Y.pred)
{
  # n.train >= 1
  n.train = length(Y.train)
  n.pred = length(Y.pred)

  muhat.train = matrix(0, n.train, M)   # train y
  ehat.train = matrix(0, n.train, M)    #
  cvehat.train = matrix(0, n.train, M)  #
  muhat.pred = matrix(0, n.pred, M)     # pred y
  prederr = matrix(0, n.pred, M)        #
  edf = rep(0, M)                       # tr(P)
  folds = cvfolds(nfolds, n.train)

  ######## main ##########
  for(j in 1:M){
    ### data
    temp1 = a3[j, 1:nump]               # variable index of linear part
    temp2 = a3[j, (nump+1):(nump+numq)] # variable index of non-linear part
    temp1 = temp1[temp1 != 0]
    temp2 = temp2[temp2 != 0]
    nx = length(temp1)
    nz = length(temp2)
    X.train = as.matrix(XX.train[,temp1])        # n*length(temp1)
    Z.train = as.matrix(XX.train[,temp2])        # n*length(temp2)
    X.pred = as.matrix(XX.pred[,temp1])
    Z.pred = as.matrix(XX.pred[,temp2])

    ### fitting
    fit = plam.fit(X.train, Z.train, Y.train, folds, X.pred, Z.pred, Y.pred)
    muhat.train[,j] = c(fit$muhat.train)
    ehat.train[,j] = Y.train - muhat.train[,j]
    muhat.pred[,j] = c(fit$muhat.pred)
    prederr[,j] = Y.pred - muhat.pred[,j]
    cvehat.train[,j] = c(fit$cvehat.train)
    edf[j] = fit$edf
  } # end j, M

  sigmahat = apply(ehat.train, 2, crossprod) / n.train
  sigmahat1 = ehat.train[,1]^2     # 1 * nt, max model residual square

  ###### methods and risks ########
  aic = n.train*log(sigmahat) + 2*edf     # 1*M
  T1 = min(aic); t1 = which.min(aic)
  aic = aic - T1
  bic = n.train*log(sigmahat) + log(n.train)*edf   # 1*M
  T2 = min(bic); t2 = which.min(bic)
  bic = bic - T2

  ##### AIC
  muhataic =  muhat.train[,t1]            # 1*nt
  RiskAIC = mean((Y.train-muhataic)^2)    # 1*1
  predAIC = mean(prederr[,t1]^2)          # 1*1
  ##### BIC
  muhatbic = muhat.train[,t2]
  RiskBIC = mean((Y.train-muhatbic)^2)    # 1*1
  predBIC = mean(prederr[,t2]^2)          # 1*1
  ##### SAIC
  wa = exp(-aic/2) / sum(exp(-aic/2))   # 1*M
  wa[wa < 1e-8] = 0
  wa = wa / sum(wa)
  muhatsaic = c(muhat.train%*%wa)           # n*M*1 -> n*1 matrix
  RiskSAIC = mean((Y.train-muhatsaic)^2)    # 1*1
  predSAIC = mean(c(prederr%*%wa)^2)        # 1*1
  ##### SBIC
  wb = exp(-bic/2) / sum(exp(-bic/2))   # 1*M
  wb[wb < 1e-8] = 0
  wb = wb / sum(wb)
  muhatsbic = c(muhat.train%*%wb)           # n*M*1 -> n*1 matrix
  RiskSBIC = mean((Y.train-muhatsbic)^2)    # 1*1
  predSBIC = mean(c(prederr%*%wb)^2)        # 1*1
  ##### Cross-Validation
  cvehat.train[is.nan(cvehat.train)] = 0    # nt*M
  Qn = t(cvehat.train)%*%cvehat.train       # M*M
  Qn = (Qn + t(Qn)) / 2
  if(qr(Qn)$rank < M){ Qn = Qn + diag(1e-8, M, M) }
  Am2 = matrix(0, M, M+1)
  Am2[,1] = rep(1, M)            # t(1)*ww = 1
  Am2[,2:(M+1)] = diag(M)        # I*ww >= 0
  bm2 = c(1, rep(0,M))
  wcv = solve.QP(Qn/n.train, rep(0, M), Am2, bm2, meq = 1)$solution  # quadratic programming
  wcv[wcv < 1e-8] = 0
  wcv = wcv / sum(wcv)
  muhatcv = c(muhat.train%*%wcv)        # nt*1 matrix->1*n vector
  Riskcv = mean((Y.train-muhatcv)^2)    # 1*1
  predcv = mean(c(prederr%*%wcv)^2)     # 1*1

  return(list(aic = RiskAIC, bic = RiskBIC,
              saic = RiskSAIC, sbic = RiskSBIC, cv = Riskcv,
              muhat = muhat.train,
              ws = list(t1, t2, wa, wb, wcv),
              predaic = predAIC, predbic = predBIC,
              predsaic = predSAIC, predsbic = predSBIC, predcv = predcv,
              mutesthat = muhat.pred))
}
