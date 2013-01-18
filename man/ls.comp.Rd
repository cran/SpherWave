\name{ls.comp}
\alias{ls.comp}
\title{Computation of Coefficients of Multi-scale SBF Representation}
\description{
This function computes coefficients of multi-scale SBF representation by least squares method.
}
\usage{
ls.comp(obs, site, ssite, snet, seta)  
}
\arguments{
\item{obs}{observations}
\item{site}{grid points of observation sites in radian for computing coefficients}
\item{ssite}{grid points of observation sites in radian used in least squares method}
\item{snet}{vector of labels representing sub-networks} 
\item{seta}{bandwidth parameters for Poisson kernel}
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
\code{\link{lscoef.comp}}, \code{\link{ridge.comp}}
}
\keyword{nonparametric}
