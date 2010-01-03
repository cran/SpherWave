\name{mrs.comp.thresh.global}
\alias{mrs.comp.thresh.global}
\title{Global Thresholding of SW Coefficients}
\description{
This function calculates global thresholded SW coefficients.
}
\usage{
mrs.comp.thresh.global(coef, site, netlab, eta, K, lam, type) 
}
\arguments{
\item{coef}{coefficients of multi-scale SBF's}
\item{site}{grid points of observation sites in radian}
\item{netlab}{vector of labels representing sub-networks} 
\item{eta}{bandwidth parameters for Poisson kernel} 
\item{K}{the number of resolution levels to be thresholded in the decomposition}
\item{lam}{thresholding value}
\item{type}{the type of thresholding. This can be `"universal"', `"probability"', `"fdr"', `"Lorentz"' and `"sure"'.}
}
\details{
For selective reconstruction, this function performs thresholding SW coefficients according to several approaches.  
}
\value{ 
\item{talpha}{global thresholded SW coefficients}
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
\code{\link{mrs.comp.thresh.level}},  \code{\link{mrsfield.comp.thresh.global}}
}
\keyword{nonparametric}
