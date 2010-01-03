\name{lscoef.comp}
\alias{lscoef.comp}
\title{Computation of Interpolation Coefficients of Multi-scale SBF Representation}
\description{
This function computes interpolation coefficients of multi-scale SBF representation.
}
\usage{
lscoef.comp(site1, site2, ssite1, ssite2, snet, seta)  
}
\arguments{
\item{site1}{latitudes of observation sites in radian for computing coefficients}
\item{site2}{longitudes of observation sites in radian for computing coefficients}
\item{ssite1}{latitudes of observation sites in radian used in least squares method}
\item{ssite2}{longitudes of observation sites in radian used in least squares method}
\item{snet}{vector of labels representing sub-networks} 
\item{seta}{bandwidth parameters for Poisson kernel}
}
\value{
\item{gg}{interpolation coefficients of multi-scale SBF representation}
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
\code{\link{ls.comp}}, \code{\link{ridge.comp}}.
}
\keyword{nonparametric}
