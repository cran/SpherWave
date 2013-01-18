\name{cov.comp}
\alias{cov.comp}
\title{Computation of Covariance Matrix}
\description{
This function generates full covariance matrix of data based on SBF.
}
\usage{
cov.comp(site1, site2, eta, netlab)
}
\arguments{
\item{site1}{latitudes of observation sites in radian}
\item{site2}{longitudes of observation sites in radian}
\item{netlab}{vector of labels representing sub-networks} 
\item{eta}{bandwidth parameters for Poisson kernel} 
}
\details{
For details, see references below.
}
\value{ 
\item{aa}{covariance matrix}
}
\references{
Oh, H-S. (1999)  Spherical wavelets and their statistical analysis with applications to meteorological data. Ph.D. Thesis, 
Department of Statistics, Texas A\&M University, College Station.

Li, T-H. (1999) Multiscale representation and analysis of spherical data by spherical wavelets. 
\emph{SIAM Journal on Scientific Computing}, \bold{21}, 924--953.
}
\seealso{
\code{\link{mcov.comp}}
}
\keyword{nonparametric}
