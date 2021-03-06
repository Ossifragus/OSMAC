% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/twostep.r
\name{twostep}
\alias{twostep}
\title{The twostep algorithm}
\usage{
twostep(X, Y, r1, r2, method = c("mvc", "mmse", "uni"))
}
\arguments{
\item{X}{the input covariate matrix}

\item{Y}{the input response vector}

\item{r1}{the first step sample size}

\item{r2}{the second step sample size}

\item{method}{the method to use}
}
\description{
This function implement the OSMAC method for the input covariate matrix @param X, response vector Y, first step sample size r1, the second step sample size r2, and the method to use.
It returns a list with three elements: par, the weighted MLE; se, the standard errors; msg, the fitting message; iter, the number of itterations used; method, the method used.
}
\examples{
library(OSMAC)
dat <- adult.train
X <- as.matrix(dat[,c(1,3,5,12:13)])
X <- t(t(X) / apply(X, 2, sd))
X <- cbind(1, X)
Y <- as.numeric(dat[,15]) - 1
set.seed(0)
twostep(X, Y, 200, 800, "mmse")
}
\keyword{getMLE}
\keyword{twostep}

