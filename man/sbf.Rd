\name{sbf}
\alias{sbf}
\title{Extrapolation with Multi-sale SBF's}
\description{
This function performs extrapolation with multi-sale SBF's.
}
\usage{
sbf(obs, latlon, netlab, eta, method, approx=FALSE,
    grid.size=c(50, 100), lambda=NULL, p0=0, latlim=NULL, 
    lonlim=NULL) 
}
\arguments{
\item{obs}{observations}
\item{latlon}{grid points of observation sites in degree.
    Latitude is the angular distance in degrees of a point north or south of the Equator. 
    North/South are represented by +/- sign. 
    Longitude is the angular distance in degrees of a point east or west of the Prime (Greenwich) Meridian. 
    East/West are represented by +/- sign.}
\item{netlab}{vector of labels representing sub-networks} 
\item{eta}{bandwidth parameters for Poisson kernel} 
\item{method}{extrapolation methods, `"ls"' or `"pls"'} 
\item{approx}{if TRUE, approximation is used.}
\item{grid.size}{grid size (latitude, longitude) of extrapolation site} 
\item{lambda}{smoothing parameter for penalized least squares method} 
\item{p0}{specifies starting level for extrapolation. Among resolution levels \eqn{1, \ldots, L},  
resolution levels \eqn{p0+1, \ldots, L} will be included for extrapolation.} 
\item{latlim}{range of latitudes in degree}
\item{lonlim}{range of longitudes in degree}
}
\details{
This function performs extrapolation with multi-sale SBF's.
}
\value{
An object of class `sbf'. This object is a list with the following components. 
\item{obs}{observations}
\item{latlon}{grid points of observation sites in degree}
\item{netlab}{vector of labels representing sub-networks} 
\item{eta}{bandwidth parameters for Poisson kernel} 
\item{method}{extrapolation methods, `"ls"' or `"pls"'} 
\item{approx}{if TRUE, approximation is used.}
\item{grid.size}{grid size (latitude, longitude) of extrapolation site} 
\item{lambda}{smoothing parameter for penalized least squares method} 
\item{p0}{starting level for extrapolation. Resolution levels \eqn{p0+1, \ldots, L} is used for extrapolation.} 
\item{gridlon}{longitudes of extrapolation sites in degree}
\item{gridlat}{latitudes of extrapolation sites in degree}
\item{nlevels}{the number of multi-resolution levels} 
\item{coeff}{interpolation coefficients}
\item{field}{extrapolation on grid.size}
\item{density}{density on observation's locations}
\item{latlim}{range of latitudes in degree}
\item{lonlim}{range of longitudes in degree}
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
\code{\link{swd}}, \code{\link{swthresh}}, \code{\link{swr}}.
}
\examples{
### Observations of year 1967
#data(temperature)
#names(temperature)

# Temperatures on 939 weather stations of year 1967    
#temp67 <- temperature$obs[temperature$year == 1967] 
# Locations of 939 weather stations    
#latlon <- temperature$latlon[temperature$year == 1967, ]

### Network design by BUD
#data(netlab)

### Bandwidth for Poisson kernel
#eta <- c(0.961, 0.923, 0.852, 0.723, 0.506)

### SBF representation of the observations by pls
#out.pls <- sbf(obs=temp67, latlon=latlon, netlab=netlab, eta=eta, 
#    method="pls", grid.size=c(50, 100), lambda=0.89)
}
\keyword{nonparametric}
