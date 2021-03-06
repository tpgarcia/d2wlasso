% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/main_functions.R
\name{store.x}
\alias{store.x}
\title{Matrix to store results for covariates x.}
\usage{
store.x(XX)
}
\arguments{
\item{XX}{(n by m) matrix of main covariates where m is the number of covariates and n is the sample size.}
}
\value{
(m by 1) matrix of zeroes
}
\description{
Forms a matrix that will be used to store regression results of response y on covariates x.
}
\examples{

xx = matrix(rnorm(100*5, 0, 1),100,5)
colnames(xx) <- paste0("X",1:ncol(xx))
store.x(xx)
}
