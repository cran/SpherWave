\name{ridge.comp}
\alias{ridge.comp}
\title{Computation of Coefficients of Multi-scale SBF's by Ridge Regression}
\description{
This function computes coefficients of multi-scale SBF's by ridge regression.
}
\usage{
ridge.comp(obs, site, ssite, snet, seta, lam)
}
\arguments{
\item{obs}{observations}
\item{site}{grid points in radian for computing coefficients}
\item{ssite}{grid points of observation sites in radian used in ridge regression}
\item{snet}{vector of labels representing sub-networks} 
\item{seta}{bandwidth parameters for Poisson kernel}
\item{lam}{smoothing parameter for ridge regression}
}

\value{
An object of class `lsfit'.
}
\references{
Oh, H-S. (1999)  Spherical wavelets and their statistical analysis with applications to meteorological data. Ph.D. Thesis, 
Department of Statistics, Texas A\&M University, College Station.

Li, T-H. (1999) Multiscale representation and analysis of spherical data by spherical wavelets. 
\emph{SIAM Journal on Scientific Computing}, \bold{21}, 924--953.

Oh, H-S. and Li, T-H. (2004) Estimation of global temperature fields from scattered observations by 
a spherical-wavelet-based spatially adaptive method. \emph{Journal of the Royal Statistical Society
Ser.} B, \bold{66}, 221--238.
}
\seealso{
\code{\link{lscoef.comp}}, \code{\link{gg.comp}}, \code{\link{ls.comp}}.
}
\keyword{nonparametric}
