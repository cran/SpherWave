\name{swthresh}
\alias{swthresh}
\title{Thresholding of Spherical Wavelet Decomposition (`swd') Object}
\description{
This function performs various ways to threshold a `swd' class object.
}
\usage{
swthresh(swd, policy, by.level, type, nthresh, value = 0.1, 
Q = 0.05) 
}
\arguments{
\item{swd}{an object of class `swd'}
\item{policy}{threshold technique. At present the possible policies are
    `"universal"', `"probability"', `"fdr"', `"Lorentz"' and `"sure"'.}
\item{by.level}{If FALSE, then perform a global threshold. If TRUE, a thresholding value
    is computed and applied separately to each resolution level.}
\item{type}{the type of thresholding. This can be `"hard"', `"soft"' or `"Lorentz"'.}
\item{nthresh}{the number of resolution levels to be thresholded in the decomposition}
\item{value}{the user supplied threshold represented by quantile level for `"probability"' policy}
\item{Q}{parameter for the false discovery rate of `"fdr"' policy}
}
\details{
This function thresholds or shrinks details stored in a `swd' object
and returns the thresholded details in a modified `swd' object.
For level-dependent thresholding, `"universal"', `"Lorentz"' and `"fdr"' are provided.
Only hard type thresholding is proper for `"probability"' thresholding.
Also note that only soft type thresholding is proper for `"sure"' thresholding.
}
\value{ 
An object of class `swd'. This object is a list with the following components. 
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
\item{swcoeff}{spherical wavelet coefficients}
\item{thresh.info}{thresholding information}
}
\references{
Oh, H-S. and Li, T-H. (2004) Estimation of global temperature fields from scattered observations by 
a spherical-wavelet-based spatially adaptive method. \emph{Journal of the Royal Statistical Society
Ser.} B, \bold{66}, 221--238.
}
\seealso{
\code{\link{sbf}}, \code{\link{swd}}, \code{\link{swr}}.
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

### Thresholding
#out.univ <- swthresh(out.dpls, policy="universal", by.level=TRUE, 
#    type="hard", nthresh=4) 
}
\keyword{nonparametric}
