\name{mrsfield.comp.thresh.level}
\alias{mrsfield.comp.thresh.level}
\title{Generation of Detailed Fields by Level-dependent Thresholding}
\description{
This function generates detailed fields based on level-dependent thresholding of SW coefficients.
}
\usage{
mrsfield.comp.thresh.level(grid, coef, site, netlab, eta, K, 
policy, Q, type)
}
\arguments{
\item{grid}{grid points of extrapolation sites in radian}
\item{coef}{coefficients of multi-scale SBF's}
\item{site}{grid points of observation sites in radian}
\item{netlab}{vector of labels representing sub-networks} 
\item{eta}{bandwidth parameters for Poisson kernel} 
\item{K}{the number of resolution levels to be thresholded in the decomposition}
\item{policy}{threshold technique. At present the possible policies are
    `"universal"', `"fdr"' and `"Lorentz"'.}
\item{Q}{parameter for the false discovery rate of `"fdr"' policy.}
\item{type}{the type of thresholding. This can be `"hard"', `"soft"' or `"Lorentz"'.}
}
\details{
This function calculates level-dependent thresholded detailed fields.
}
\value{ 
\item{dfield}{level-dependent thresholded detailed fields}
}
\references{
Oh, H-S. (1999)  Spherical wavelets and their statistical analysis with applications to meteorological data. Ph.D. Thesis, 
Department of Statistics, Texas A\&M University, College Station.

Oh, H-S. and Li, T-H. (2004) Estimation of global temperature fields from scattered observations by 
a spherical-wavelet-based spatially adaptive method. \emph{Journal of the Royal Statistical Society
Ser.} B, \bold{66}, 221--238.
}
\seealso{
\code{\link{sbf}}, \code{\link{swd}}, \code{\link{swthresh}}, \code{\link{swr}}
}
\keyword{nonparametric}
