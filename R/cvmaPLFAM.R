#' @title Cross-Validation Model Averaging (CVMA) for Partial Linear Functional Additive Models (PLFAMs)
#' @description Summarize the estimate of weights for averaging across all candidate models for PLFAMs, using multi-fold cross-validation criterion,
#'      and the corresponding mean squared prediction error risk. Additionally, the results of AIC, BIC, SAIC and SBIC are delivered simultaneously.
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
#'     \item{aic}{Mean squared error risk in training data set, produced by AIC model selection method.}
#'     \item{bic}{Mean squared error risk in training data set, produced by BIC model selection method.}
#'     \item{saic}{Mean squared error risk in training data set, produced by SAIC model averaging method.}
#'     \item{sbic}{Mean squared error risk in training data set, produced by SBIC model averaging method.}
#'     \item{cv}{Mean squared error risk in training data set, produced by CVMA method.}
#'     \item{waic}{The selected candidate model by AIC model selection method.}
#'     \item{wbic}{The selected candidate model by BIC model selection method.}
#'     \item{wsaic}{The weights for each candidate model by SAIC model averaging method.}
#'     \item{wsbic}{The weights for each candidate model by SBIC model averaging method.}
#'     \item{wcv}{The weights for each candidate model by CVMA method.}
#'     \item{predaic}{Mean squared prediction error risk in test data set, produced by AIC model selection method.}
#'     \item{predbic}{Mean squared prediction error risk in test data set, produced by BIC model selection method.}
#'     \item{predsaic}{Mean squared prediction error risk in test data set, produced by SAIC model averaging method.}
#'     \item{predsbic}{Mean squared prediction error risk in test data set, produced by SBIC model averaging method.}
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
#' simdata = data_gen(R = 0.7, K = 1, n = 50, M0 = 20, typ = 1, design = 3)
#' dat1 = simdata[[1]]
#' scalars = dat1[,1:20]
#' fd = dat1[,21:120]
#' Y = dat1[,122]
#' tps = seq(0, 1, length.out = 100)
#'
#' # Estimation
#' est_res = cvmaPLFAM(Y=Y, scalars = scalars, functional = fd, tt = tps,
#'        nump = 2, numfpcs = 3, nbasis = 50, nfolds = 5, ratio.train = 0.8)
#' # Weights estimated by CVMA method
#' est_res$wcv
#' # Prediction error risk on test data set
#' est_res$predcv
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

  # Risk0 = rep(NA, M)
  RiskAIC = NA
  RiskBIC = NA
  RiskSAIC = NA
  RiskSBIC = NA
  RiskCV = NA

  # mspe0 = rep(NA, M)
  predAIC = NA
  predBIC = NA
  predSAIC = NA
  predSBIC = NA
  predCV = NA

  waic = NA        # model selected by aic
  wbic = NA        # model selected bt bic
  wsaic = rep(NA, M)    # model weights in saic
  wsbic = rep(NA, M)    # model weights in sbic
  wcv = rep(NA, M)      # model weights in cv

  ########## main #############
    ######## data preparing
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

    #### train
    trfpcalist = fpcscore(functional.train, nbasis, tt)
    trorgscs = trfpcalist$score                    # n by nbasis
    trorglam = trfpcalist$eigv                     # nbasis
    #### transformed FPCs, from (-Inf,Inf) to [-1,1]
    trZ10 = pnorm(trorgscs[,1:max(10, numfpcs)], sd = rep(sqrt(trorglam[1:max(10, numfpcs)]), each = nrow(trorgscs)))
    #### test
    ttfpcalist = fpcscore(functional.test, nbasis, tt)
    ttorgscs = ttfpcalist$score
    ttorglam = ttfpcalist$eigv
    ttZ10 = pnorm(ttorgscs[,1:max(10, numfpcs)], sd = rep(sqrt(ttorglam[1:max(10, numfpcs)]), each = nrow(ttorgscs)))

    ######## design matrix Z = (X, fpc)
    trXX = cbind(scalars.train[,1:nump], trZ10[,1:numfpcs]) # n*x matrix
    testXX = cbind(scalars.test[,1:nump], ttZ10[,1:numfpcs])

    ########### methods
    allR = predRisk(M, nump, numfpcs, vars$a2, vars$a3, nfolds = nfolds, XX.train = trXX,
                    Y.train =  Y.train, XX.pred = testXX, Y.pred = Y.test)
    #### train data
    RiskAIC = allR$aic
    RiskBIC = allR$bic
    RiskSAIC = allR$saic
    RiskSBIC = allR$sbic
    RiskCV = allR$cv
    # muhat = allR$muhat   # n*M matrix

    waic = allR$ws[[1]]
    wbic = allR$ws[[2]]
    wsaic = allR$ws[[3]]
    wsbic = allR$ws[[4]]
    wcv = allR$ws[[5]]

    #### pred data
    predAIC = allR$predaic
    predBIC = allR$predbic
    predSAIC = allR$predsaic
    predSBIC = allR$predsbic
    predCV = allR$predcv
    # mutesthat = allR$mutesthat

    # ####### Risks for all candidate models
    # for(i in 1:M){
    #   Risk0[i] = mean((Y.train - muhat[,i])^2)     # 1*1
    #   mspe0[i] = mean((Y.test - mutesthat[,i])^2)
    # }

  # Risk = min(apply(Risk0, 2, mean))    # 1*1
  # mspe = min(apply(mspe0, 2, mean))

  return(list(aic = RiskAIC, bic = RiskBIC, saic = RiskSAIC, sbic = RiskSBIC,
              cv = RiskCV,
              waic = waic, wbic = wbic, wsaic = wsaic, wsbic = wsbic,
              wcv = wcv,
              predaic = predAIC, predbic = predBIC,
              predsaic = predSAIC, predsbic = predSBIC, predcv = predCV))
}

