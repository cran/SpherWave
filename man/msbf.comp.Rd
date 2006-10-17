\name{msbf.comp}
\alias{msbf.comp}
\title{Calculation of an Extrapolation Field with Multiscale SBF's}
\description{
This function calculates an extrapolation field evaluated on grid with multiscale SBF's. It also provides
the density of SBF on grid.  
}
\usage{
msbf.comp(grid, site, coef, netlab, eta, p0) 
}
\arguments{
\item{grid}{grid points of extrapolation sites in radian}
\item{site}{grid points of observation sites in radian}
\item{coef}{coefficients of multi-scale SBF's}
\item{netlab}{vector of labels representing sub-networks} 
\item{eta}{bandwidth parameters for Poisson kernel} 
\item{p0}{(p0+1) will be starting detailed level to be included} 
}
\value{
\item{field}{extrapolation field}
\item{density}{density of SBF}
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
\code{\link{sbf.comp}}
}
\keyword{nonparametric}
