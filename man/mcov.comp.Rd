\name{mcov.comp}
\alias{mcov.comp}
\title{Computation of Covariance Matrix}
\description{
This function generates full covariance matrix of data based on SBF.
}
\usage{
mcov.comp(site, netlab, eta)
}
\arguments{
\item{site}{a J x 2 matrix of grid points of observation sites in radian}
\item{netlab}{vector of labels representing sub-networks} 
\item{eta}{bandwidth parameters for Poisson kernel} 
}
\value{ 
\item{cov}{covariance matrix}
}
\references{
Oh, H-S. (1999)  Spherical wavelets and their statistical analysis with applications to meteorological data. Ph.D. Thesis, 
Department of Statistics, Texas A\&M University, College Station.

Li, T-H. (1999) Multiscale representation and analysis of spherical data by spherical wavelets. 
\emph{SIAM Journal on Scientific Computing}, \bold{21}, 924--953.
}
\seealso{
\code{\link{mcov.comp}}, \code{\link{mracoef.comp}}
}
\keyword{nonparametric}
