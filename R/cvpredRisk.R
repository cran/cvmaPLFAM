#' @title Output the prediction risks of the cross-validation model averaging (CVMA) method for partial linear functional additive models (PLFAMs)
#' @description Calculate the estimated weights for averaging across all candidate models and the corresponding
#'      mean squared prediction error risk.
#'
#' @param M The number of candidate models.
#' @param nump The number of scalar predictors in candidate models.
#' @param numq The number of funtional principal components (FPCs) in candidate models.
#' @param a2 The number of FPCs in each candidate model. See \code{\link{modelspec}}.
#' @param a3 The index for each component in each candidate model. See \code{\link{modelspec}}.
#' @param nfolds The number of folds used in cross-validation.
#' @param X.train The training data of scalar predictors.
#' @param ZZ.train The training data of the functional predictor.
#' @param Y.train The training data of response variable.
#' @param X.pred The test data of scalar predictors.
#' @param ZZ.pred The test data of the functional predictor.
#' @param Y.pred The test data of response variable.
#' @param nbasis The number of basis functions used for spline approximation.
#' @param tt The vector of recording/measurement points for the functional predictor.
#'
#' @return A \code{list} of
#'     \item{cv}{Mean squared error risk in training data set, produced by CVMA method.}
#'     \item{ws}{A \code{vector} of weights estimator.}
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

cvpredRisk <- function(M, nump, numq, a2, a3, nfolds, X.train, ZZ.train, Y.train,
                     X.pred, ZZ.pred, Y.pred, nbasis, tt)
{
  # n.train >= 1
  n.train = length(Y.train)
  n.pred = length(Y.pred)

  # muhat.train = matrix(0, n.train, M)   # train y
  # ehat.train = matrix(0, n.train, M)    #
  cvehat.train = matrix(0, n.train, M)  #
  # muhat.pred = matrix(0, n.pred, M)     # pred y
  # prederr = matrix(0, n.pred, M)        #
  # edf = rep(0, M)                       # tr(P)
  folds = cvfolds(nfolds, n.train)

  ##### fitting
  fit = plam.fit(M, nump, numq, a3, X.train, ZZ.train, Y.train, X.pred, ZZ.pred, Y.pred, nbasis, tt)
  muhat.train = fit$muhat.train     # nt * M
  ehat.train = fit$ehat.train       # nt * M
  muhat.pred = fit$muhat.pred       # np * M
  prederr = fit$prederr             # np * M
  edf = fit$edf                     # 1 * M

  ##### cv error
  for(ii in 1:length(folds)){
    # first split records data into training set and test set
    # then transfer to plam.fit function
    # further apply fpcscores function to extract fpc scores
    # this ensures FPCA will not learn from the test fold data
    fit1 <- plam.fit(M, nump, numq, a3, X.train = X.train[-folds[[ii]],], ZZ.train = ZZ.train[-folds[[ii]],],
                      y.train = Y.train[-folds[[ii]]], X.pred = X.train[folds[[ii]],],
                      ZZ.pred = ZZ.train[folds[[ii]],], y.pred = Y.train[folds[[ii]]], nbasis, tt)
    # cvpred <- fit1$muhat.pred     # one fold length-by-M matrix
    cvehat.train[folds[[ii]], ] <- fit1$prederr
  }

  sigmahat = apply(ehat.train, 2, crossprod) / n.train
  sigmahat1 = ehat.train[,1]^2     # 1 * nt, max model residual square

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
  wcv[wcv < 0] = 0
  wcv = wcv / sum(wcv)
  muhatcv = c(muhat.train%*%wcv)        # nt*1 matrix->1*n vector
  Riskcv = mean((Y.train-muhatcv)^2)    # 1*1
  predcv = mean(c(prederr%*%wcv)^2)     # 1*1

  return(list(cv = Riskcv, ws = wcv, predcv = predcv))

}
