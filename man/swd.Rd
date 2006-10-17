\name{swd}
\alias{swd}
\title{Decomposition}
\description{
This function performs decomposition with multi-sale SBF's.
}
\usage{
swd(sbf) 
}
\arguments{
\item{sbf}{an object of class `sbf'}
}
\details{
This function performs decomposition with multi-sale SBF's.
}
\value{ 
An object of class spherical wavelet decomposition(`swd'). This object is a list with the following components. 
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
\item{density1}{density of SBF}
\item{latlim}{range of latitudes in degree}
\item{lonlim}{range of longitudes in degree}
\item{global}{List of successively smoothed data}
\item{density}{density of SW coefficients}
\item{detail}{List of details at different resolution levels}
\item{swcoeff}{SW coefficients}
\item{thresh.info}{`"None"'}
}
\references{
Oh, H-S. and Li, T-H. (2004) Estimation of global temperature fields from scattered observations by 
a spherical-wavelet-based spatially adaptive method. \emph{Journal of the Royal Statistical Society
Ser.} B, \bold{66}, 221--238.
}
\seealso{
\code{\link{sbf}}, \code{\link{swthresh}}, \code{\link{swr}}.
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

### Decomposition
#out.dpls <- swd(out.pls)
}
\keyword{nonparametric}
