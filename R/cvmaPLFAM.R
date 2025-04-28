#' @title Cross-Validation Model Averaging (CVMA) for Partial Linear Functional Additive Models (PLFAMs)
#' @description Summarize the estimate of weights for averaging across all candidate models for PLFAMs, using multi-fold cross-validation criterion,
#'      and the corresponding mean squared prediction error risk.
#'
#' @param Y The vector of the scalar response variable.
#' @param scalars The design matrix of scalar predictors.
#' @param functional The matrix including records/measurements of the functional predictor.
#' @param Y.test Test data: The vector of the scalar response variable.
#' @param scalars.test Test data: The design matrix of scalar predictors.
#' @param functional.test Test data: The matrix including records/measurements of the functional predictor.
#' @param tt The vector of recording/measurement points for the functional predictor.
#' @param nump The number of scalar predictors in candidate models.
#' @param numfpcs The number of functional principal components (FPCs) for the functional predictor in candidate models.
#' @param nbasis The number of basis functions used for spline approximation.
#' @param nfolds The number of folds used in cross-validation.
#' @param ratio.train The ratio of data for training, if test data are \code{NULL}.
#'
#' @return A \code{list} of
#'     \item{cv}{Mean squared error risk in training data set, produced by CVMA method.}
#'     \item{wcv}{The weights for each candidate model by CVMA method.}
#'     \item{predcv}{Mean squared prediction error risk in test data set, produced by CVMA method.}
#' @export
#'
#' @import fda
#' @import quadprog
#' @import mgcv
#' @import stats
#' @importFrom utils combn
#'
#' @examples
#' # Generate simulated data
#' simdata = data_gen(R = 0.7, K = 1, n = 50, ntest = 10, M0 = 4, typ = 1, design = 1)
#' train_dat = simdata[[1]]
#' scalars.train = train_dat[,1:4]
#' fd.train = train_dat[,5:104]
#' Y.train = train_dat[,106]
#'
#'
#' test_dat = simdata[[2]]
#' scalars.test = test_dat[,1:4]
#' fd.test = test_dat[,5:104]
#' Y.test = test_dat[,106]
#'
#' tps = seq(0, 1, length.out = 100)
#'
#' # Estimation
#' res = cvmaPLFAM(Y=Y.train, scalars = scalars.train, functional = fd.train,
#' Y.test = Y.test, scalars.test = scalars.test, functional.test = fd.test, tt = tps,
#'        nump = 2, numfpcs = 3, nbasis = 50, nfolds = 5)
#' # Weights estimated by CVMA method
#' res$wcv
#' # Prediction error risk on test data set
#' res$predcv
#'

cvmaPLFAM <- function(Y, scalars, functional, Y.test = NULL, scalars.test = NULL,
                  functional.test = NULL, tt, nump, numfpcs, nbasis, nfolds, ratio.train = NULL)
{
  if(is.null(Y))
    stop("The response vector \"Y\" is missing!")
  if(is.null(scalars))
    stop("The design matrix of scalar predictors \"scalars\" is missing!")
  if(is.null(functional))
    stop("The matrix of the functional predictor \"functional\" is missing!")
  if(is.null(nump))
    stop("The number of scalar predictors in candidate models \"nump\" is missing!")
  if(is.null(numfpcs))
    stop("The number of functional principal components \"numfpcs\" in candidate models is missing!")
  if(is.null(nbasis))
    stop("The number of basis functions for spline approximation \"nbasis\" is missing!")
  if(is.null(nfolds))
    stop("The number of folds for cross-validation \"nfolds\" is missing!")

  # sample size
  ntotal = length(Y)
  if(ntotal <= 1)
    stop("The sample size of training data should be greater than one!")

  if(nfolds <= 1 || nfolds > ntotal)
    stop("A wrong number of folds \"nfolds\" is given!")
  if(is.null(Y.test) && is.null(ratio.train))
    stop("A ratio of training data should be given when test data are NULL!")

  if(!is.matrix(scalars))
    scalars = as.matrix(scalars)
  if(!is.matrix(functional))
    functional = as.matrix(functional)

  if(nrow(scalars) != ntotal)
    stop("Observations of response variable and scalar predictors unmatched!")
  if(nrow(functional) != ntotal)
    stop("Observations of response variable and the functional predictor unmatched!")
  if(ncol(functional) != length(tt))
    stop("Recording/measurement points in \"functional\" and \"tt\" are unmatched!")

  vars = modelspec(nump, numfpcs)   # specify candidate models
  M = length(vars$a1)               # number of candidate models

  ####### data preparing
  if(is.null(Y.test)){
    trains = sample(ntotal, ceiling(ntotal*ratio.train))
    Y.train = Y[trains]
    Y.test = Y[-trains]
    scalars.train = scalars[trains,]
    scalars.test = scalars[-trains,]
    functional.train = functional[trains,]
    functional.test = functional[-trains,]
  }else{
    Y.train = Y
    scalars.train = scalars
    functional.train = functional
  }

  cvmalist = cvpredRisk(M, nump, numfpcs, vars$a2, vars$a3, nfolds = nfolds,
                        X.train = scalars.train, ZZ.train = functional.train,
                        Y.train = Y.train, X.pred = scalars.test, ZZ.pred = functional.test,
                        Y.pred = Y.test, nbasis, tt)

  RiskCV = cvmalist$cv
  wcv = cvmalist$ws
  predCV = cvmalist$predcv

  return(list(cv = RiskCV, wcv = wcv, predcv = predCV))

}

