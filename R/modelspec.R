#' @title Generate candidate models
#' @description Specify non-nested or nested candidate models, according to the prescribed number of scalar predictors and the number of functional principal components (FPCs).
#'      Each candidate model comprises at least one scalar predictor and one FPC.
#'
#' @param nump The number of scalar predictors used in candidate models.
#' @param numq The number of functional principal components (FPCs) used in candidate models.
#' @param method A character string or NULL.
#' If \code{NULL}, candidate models are generated under a non-nested structure.
#' If \code{"nested"}, candidate models are generated under a nested structure.
#' Otherwise, an error will be raised.
#'
#' @return A \code{list} of
#'     \item{a1}{The number of scalar predictors in each candidate model.}
#'     \item{a2}{The number of FPCs in each candidate model.}
#'     \item{a3}{The index for each component in each candidate model.}
#' @export
#'
#' @importFrom utils combn
#'
#' @examples
#' # Example 1: non-nested models
#' # Given nump = 2 and numq = 2, resulting in 9 candidate models
#' modelspec(2, 2)
#' #$a1
#' #[1] 2 2 2 1 1 1 1 1 1
#' #$a2
#' #[1] 2 1 1 2 1 1 2 1 1
#' #$a3
#' #      [,1] [,2] [,3] [,4]
#' # [1,]    1    2    3    4
#' # [2,]    1    2    3    0
#' # [3,]    1    2    0    4
#' # [4,]    1    0    3    4
#' # [5,]    1    0    3    0
#' # [6,]    1    0    0    4
#' # [7,]    0    2    3    4
#' # [8,]    0    2    3    0
#' # [9,]    0    2    0    4
#'
#' # Example 2: nested models
#' # Given nump = 2 and numq = 3, resulting in 6 candidate models
#' modelspec(2, 3, method = "nested")
#' #$a1
#' # [1] 2 2 2 1 1 1
#' #$a2
#' # [1] 3 2 1 3 2 1
#' #$a3
#' #      [,1] [,2] [,3] [,4] [,5]
#' # [1,]    1    2    3    4    5
#' # [2,]    1    2    3    4    0
#' # [3,]    1    2    3    0    0
#' # [4,]    1    0    3    4    5
#' # [5,]    1    0    3    4    0
#' # [6,]    1    0    3    0    0
#'
#'

modelspec <- function(nump, numq, method = NULL)
{
  if(is.null(method)){
    # non-nested
    k1 = choose(nump, nump:max(1, min(0,nump-2)))
    sk1 = sum(k1)
    k2 = choose(numq, numq:max(1, min(0,numq-2)))
    sk2 = sum(k2)
    a1 = rep(nump:max(1, min(0,nump-2)), times = k1*sk2)
    a2 = rep(rep(numq:max(1, min(0,numq-2)), times = k2), sk1)
    a3 = matrix(NA, nrow = sk1*sk2, ncol = nump + numq)

    indl = 0
    indr = choose(nump, nump)
    nm = nump
    while( nm > 0 ){
      a3[(1+indl*sk2):(indr*sk2), 1:nump] = apply(t((1:nump)*combn(1:nump, nm, tabulate, nbins = nump)), 2, rep, each = sk2)
      nm = nm - 1
      indl = indr
      indr = indr + choose(nump, nm)
    }

    a33 = matrix(NA, nrow = sk2, ncol = numq)
    nl = 0
    nr = choose(numq, numq)
    m = numq
    while( m > 0 ){
      a33[(1+nl):nr,] = t( ((nump+1):(nump+numq)) * combn(1:numq, m, tabulate, nbins = numq) )
      m = m - 1
      nl = nr
      nr = nr + choose(numq, m)
    }
    a3[, (nump+1):(nump+numq)] = apply(a33, 2, rep, times = sk1)

  }else if(method == "nested"){
    ### nested
    a1 = rep(nump:1, each = numq)
    a2 = rep(numq:1, nump)
    a3 = matrix(0, nrow = nump*numq, ncol = nump + numq)
    a33 = matrix(0, numq, numq)
    for(i in 1:numq){
      a33[i,] = c(rep(1, numq+1-i), rep(0,i-1))
    }
    a33 = matrix(nump+(1:numq), nrow = numq, ncol = numq, byrow = T)*a33

    for(i in 1:nump){
      a3[(1+(i-1)*numq):(i*numq),] = cbind(matrix(c(1:(nump+1-i), rep(0, i-1)), nrow = numq, ncol = nump, byrow = T), a33)
    }
  }else {
    ### otherwise
    stop("Invalid value for 'method'. Must be NULL or 'nested'.")
  }

  return(list(a1 = a1, a2 = a2, a3 = a3))
}
