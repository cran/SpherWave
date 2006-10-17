\name{mracoef.comp}
\alias{mracoef.comp}
\title{Computation of Global and Local Coefficients of Multiscale SBF Representation}
\description{
This function generates global and local coefficients of multiscale SBF representation.
}
\usage{
mracoef.comp(coef, site, netlab, eta)
}
\arguments{
\item{coef}{coefficients of the initial SBF representation}
\item{site}{grid points of observation sites in radian}
\item{netlab}{vector of labels representing sub-networks} 
\item{eta}{bandwidth parameters for Poisson kernel} 
}
\value{ 
\item{beta}{global coefficients of multiscale SBF representation}
\item{gamma}{local coefficients of multiscale SBF representation. That is the detailed coefficients of SW}
\item{norm}{norms of SW coefficients for sub-network}
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
\code{\link{mrafield.comp}}
}
\keyword{nonparametric}
