#' @title Generate cross-validation folds
#'
#' @description Randomly split the data indexes into \code{nfolds} folds.
#'
#' @param nfolds The number of folds used in cross-validation.
#' @param datasize The sample size.
#'
#' @return A \code{list}. Each element contains the index vector of sample data included in this fold.
#'
#' @export
#'
#' @examples
#' # Given sample size 20, generate 5 folds
#' set.seed(1212)
#' cvfolds(5, 20)
#' #[[1]]
#' # [1]  6 11 14 16
#' #[[2]]
#' # [1]  3  5 10 18
#' #[[3]]
#' # [1]  4  7  8 19
#' #[[4]]
#' # [1]  2  9 12 15
#' #[[5]]
#' # [1]  1 13 17 20
#'

cvfolds <- function(nfolds, datasize)
{
  cvlist = list()
  lab = rep(1:nfolds, ceiling(datasize/nfolds))[1:datasize]
  temp = sample(lab, datasize)
  x = 1:nfolds
  dataseq = 1:datasize
  cvlist = lapply(x, function(x){ dataseq[temp == x] })
  return(cvlist)
}
