#' @title  Fitting partial linear functional additive model
#' @description Calculate the prediction values and prediction errors across all candidate models.
#'
#' @param M The number of candidate models.
#' @param nump The number of scalar predictors in candidate models.
#' @param numq The number of funtional principal components (FPCs) in candidate models.
#' @param a3 The index for each component in each candidate model. See \code{\link{modelspec}}.
#' @param X.train The training data of scalar predictors.
#' @param ZZ.train The training data of the functional predictor.
#' @param y.train The training data of response variable.
#' @param X.pred The test data of scalar predictors.
#' @param ZZ.pred The test data of the functional predictor.
#' @param y.pred The test data of response variable.
#' @param nbasis The number of basis functions used for spline approximation.
#' @param tt The vector of recording/measurement points for the functional predictor.
#'
#' @return A \code{list} of
#'     \item{muhat.train}{A \code{matrix} of prediction values on training data set for \code{M} candidate models.}
#'     \item{ehat.train}{A \code{matrix} of prediction errors on training data set for \code{M} candidate models.}
#'     \item{muhat.pred}{A \code{matrix} of prediction values on test data set for \code{M} candidate models.}
#'     \item{prederr}{A \code{matrix} of prediction errors on test data set for \code{M} candidate models.}
#'     \item{edf}{A \code{vector} of effective degree of freedom for \code{M} candidate models.}
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
#'

plam.fit <- function(M, nump, numq, a3, X.train, ZZ.train, y.train,
                     X.pred, ZZ.pred, y.pred, nbasis, tt)
{
  n.train = length(y.train)
  n.pred = length(y.pred)
  muhat.train  = matrix(0, n.train, M)   # train
  ehat.train   = matrix(0, n.train, M)   #
  muhat.pred = matrix(0, n.pred, M)      # pred
  prederr    = matrix(0, n.pred, M)      #
  edf = rep(0, M)     # tr(P)

  X1.train <- as.matrix(X.train[,1:nump])
  X1.pred  <- as.matrix(X.pred[,1:nump])
  ZZ.train <- as.matrix(ZZ.train)
  ZZ.pred  <- as.matrix(ZZ.pred)

  ######### training
  fpcalist.train = fpcscore(ZZ.train, nbasis, tt)
  orgscs.train = fpcalist.train$score                    # n*nbasis
  orglam.train = fpcalist.train$eigv                     # 1*nbasis
  Z1.train = pnorm(orgscs.train[,1:numq], sd = rep(sqrt(orglam.train[1:numq]), each = nrow(orgscs.train)))       # transform to [-1,1]
  ######### predict
  fpcalist.pred = fpcscore(ZZ.pred, nbasis, tt)
  orgscs.pred = fpcalist.pred$score                    # nt*nbasis
  orglam.pred = fpcalist.pred$eigv                     # 1*nbasis
  Z1.pred = pnorm(orgscs.pred[,1:numq], sd = rep(sqrt(orglam.pred[1:numq]), each = nrow(orgscs.pred)))       # transform to [-1,1]

  #######
  for(j in 1:M){
    temp1 = a3[j, 1:nump]    # linear part, variable index
    temp2 = a3[j, -(1:nump)]    # non-linear part, variable index
    temp1 = temp1[temp1 != 0]
    temp2 = temp2[temp2 != 0]-nump
    X.train = as.matrix(X1.train[,temp1])        # n*length(temp1)
    Z.train = as.matrix(Z1.train[,temp2])        # n*length(temp2)
    X.pred = as.matrix(X1.pred[,temp1])
    Z.pred = as.matrix(Z1.pred[,temp2])

    nx <- ncol(X.train)
    nz <- ncol(Z.train)
    #######
    traindata <- data.frame(cbind(y.train, X.train, Z.train))
    colnames(traindata) <- xznam <- c("y", paste0("x", 1:nx), paste0("z", 1:nz))
    preddata <- data.frame(cbind(y.pred, X.pred, Z.pred))
    colnames(preddata) <- xznam

    plam <- as.formula(paste("y ~ ", paste(xznam[2:(1+nx)], collapse = " + "), " + ", paste0( "s(", paste(xznam[(2+nx):(1+nx+nz)], collapse = ", k = 5) + s(" ), ", k = 5)")))
    fit0 <- gam(plam, data = traindata, family = gaussian)
    pred <- predict(fit0, type = "response", newdata = preddata)

    muhat.train[,j] = c(fit0$fitted.values)
    ehat.train[,j] = y.train - c(fit0$fitted.values)
    muhat.pred[,j] = as.vector(pred)
    prederr[,j] = y.pred - as.vector(pred)
    edf[j] = sum(fit0$edf)

  }

  return(list(muhat.train = muhat.train, ehat.train = ehat.train,
              muhat.pred = muhat.pred, prederr = prederr, edf = edf))

}

