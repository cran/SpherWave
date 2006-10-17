\name{gg.comp}
\alias{gg.comp}
\title{Computation of Design Matrix induced by Multi-scale SBF's for Ridge Regression}
\description{
This function computes design matrix induced by multi-scale SBF's for ridge regression.
}
\usage{
gg.comp(site1, site2, ssite1, ssite2, snet, seta, lam) 
}
\arguments{
\item{site1}{latitudes of observation sites in radian}
\item{site2}{longitudes of observation sites in radian}
\item{ssite1}{latitudes of observation sites in radian used in least squares method}
\item{ssite2}{longitudes of observation sites in radian used in least squares method}
\item{snet}{vector of labels representing sub-networks} 
\item{seta}{bandwidth parameters for Poisson kernel}
\item{lam}{smoothing parameter for ridge regression}
}
\value{
\item{gg}{design matrix induced by multi-scale SBF's for ridge regression.}
}
\references{
Oh, H-S. (1999)  Spherical wavelets and their statistical analysis with applications to meteorological data. Ph.D. Thesis, 
Department of Statistics, Texas A\&M University, College Station.

Li, T-H. (1999) Multiscale representation and analysis of spherical data by spherical wavelets. 
\emph{SIAM Journal on Scientific Computing}, \bold{21}, 924--953.
}
\seealso{
\code{\link{lscoef.comp}}, \code{\link{gg.comp}}, \code{\link{ridge.comp}}.
}
\keyword{nonparametric}
