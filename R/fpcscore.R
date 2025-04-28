#' @title Calculate functional principal component (fpc) scores
#' @description Conduct functional principal component analysis (FPCA) on the observation matrix of the functional predictor.
#'
#' @param Z An \code{n} by \code{nT} matrix. The recording/measurement matrix of the functional predictor.
#' @param nbasis The number of basis functions used for spline approximation.
#' @param tt The vector of recording/measurement points for the functional predictor.
#'
#' @return A \code{list} of
#'      \item{score}{An \code{n} by \code{nbasis} matrix. The estimated functional principal component scores.}
#'      \item{eigv}{A vector of estimated eigen-values related to FPCA.}
#'      \item{varp}{A vector of percents of variance explained related to FPCA.}
#'
#' @export
#'
#' @import fda
#'
#' @examples
#' # Generate a recording/measurement matrix of the functional predictor
#' fddata = matrix(rnorm(1000), nrow = 10, ncol = 100)
#' tpoints = seq(0, 1, length.out = 100)
#'
#' library(fda)
#' # Using 20 basis functions for spline approximation
#' fpcscore(fddata, nbasis = 20, tt = tpoints)
#'
#'

fpcscore <- function(Z, nbasis, tt)
{
  Z = t(Z)        # nT * n matrix
  # Create basis functions
  basisT = create.bspline.basis(range(tt), nbasis)
  smoothlist = smooth.basis(tt, Z, fdParobj = basisT)

  # Fitted functional curves for design matrix Z
  zfd = smoothlist$fd
  # Functional principal component analysis
  zpca = pca.fd(zfd, nharm = nbasis)

  return(list(score = zpca$scores, eigv = zpca$values, varp = zpca$varprop))
}

