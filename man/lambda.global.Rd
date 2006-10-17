\name{lambda.global}
\alias{lambda.global}
\title{Global Thresholding Value}
\description{
This function calculates global thresholding value for spherical wavelet estimator.
}
\usage{
lambda.global(swd, policy, nthresh, value, Q)
}
\arguments{
\item{swd}{an object of class `swd'}
\item{policy}{threshold technique. At present the possible policies are
    `"universal"', `"probability"', `"fdr"', `"Lorentz"' and `"sure"'.}
\item{nthresh}{the number of resolution levels to be thresholded in the decomposition}
\item{value}{the user supplied threshold represented by quantile level of `"probability"' policy.}
\item{Q}{parameter for the false discovery rate of `"fdr"' policy.}
}
\value{ 
\item{lam}{global thresholding value}
}
\references{
Donoho, D.~L. and Johnstone, I.~M. (1994) Ideal spatial
adaptation by wavelet shrinkage. \emph{Biometrika}, \bold{81},
425--455.

Oh, H-S. (1999)  Spherical wavelets and their statistical analysis with applications to meteorological data. Ph.D. Thesis, 
Department of Statistics, Texas A\&M University, College Station.

Oh, H-S. and Li, T-H. (2004) Estimation of global temperature fields from scattered observations by 
a spherical-wavelet-based spatially adaptive method. \emph{Journal of the Royal Statistical Society
Ser.} B, \bold{66}, 221--238.
}
\seealso{
\code{\link{mrs.comp.thresh.global}}
}
\keyword{nonparametric}
