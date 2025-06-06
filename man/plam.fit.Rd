% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plam.fit.R
\name{plam.fit}
\alias{plam.fit}
\title{Fitting partial linear functional additive model}
\usage{
plam.fit(
  M,
  nump,
  numq,
  a3,
  X.train,
  ZZ.train,
  y.train,
  X.pred,
  ZZ.pred,
  y.pred,
  nbasis,
  tt
)
}
\arguments{
\item{M}{The number of candidate models.}

\item{nump}{The number of scalar predictors in candidate models.}

\item{numq}{The number of funtional principal components (FPCs) in candidate models.}

\item{a3}{The index for each component in each candidate model. See \code{\link{modelspec}}.}

\item{X.train}{The training data of scalar predictors.}

\item{ZZ.train}{The training data of the functional predictor.}

\item{y.train}{The training data of response variable.}

\item{X.pred}{The test data of scalar predictors.}

\item{ZZ.pred}{The test data of the functional predictor.}

\item{y.pred}{The test data of response variable.}

\item{nbasis}{The number of basis functions used for spline approximation.}

\item{tt}{The vector of recording/measurement points for the functional predictor.}
}
\value{
A \code{list} of
    \item{muhat.train}{A \code{matrix} of prediction values on training data set for \code{M} candidate models.}
    \item{ehat.train}{A \code{matrix} of prediction errors on training data set for \code{M} candidate models.}
    \item{muhat.pred}{A \code{matrix} of prediction values on test data set for \code{M} candidate models.}
    \item{prederr}{A \code{matrix} of prediction errors on test data set for \code{M} candidate models.}
    \item{edf}{A \code{vector} of effective degree of freedom for \code{M} candidate models.}
}
\description{
Calculate the prediction values and prediction errors across all candidate models.
}
