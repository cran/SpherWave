\name{sw.plot}
\alias{sw.plot}
\title{Plot of Observation, Network Design, Field, SW Coefficient, Decomposition or Reconstruction Result}
\description{
This function performs plotting of observation, network design, field, SW coefficients, 
decomposition or reconstruction result.
}
\usage{
sw.plot(sw = NULL, z = NULL, latlon = NULL, latlim = NULL, 
    lonlim = NULL, type = "field", nlevel = 256, pch = NULL, 
    cex = NULL, \dots)
}
\arguments{
\item{sw}{\code{sbf} or \code{swd} object}
\item{z}{observations, network design labels or reconstruction}
\item{latlon}{grid points of observation sites in degree}
\item{latlim}{range of latitudes in degree}
\item{lonlim}{range of longitudes in degree}
\item{type}{specifies the type "obs", "network", "field", "swcoeff", "decom" or "recon"}
\item{nlevel}{number of color levels used in legend strip}
\item{pch}{either an integer specifying a symbol or a single character to be used as the default in plotting points}
\item{cex}{a numerical value giving the amount by which plotting text and symbols should be scaled relative to the default}
\item{\dots}{the usual arguments to the image function or plot function}
}
\details{
This function plots spherical wavelet results. Possible types are

    `"obs"' for observations
    
    `"network"' for network design
    
    `"field"' for field
    
    `"swcoeff"' for spherical wavelet coefficients
    
    `"decom"' for decomposition results
    
    `"recon"' for reconstruction result.
    
For `sw', \code{sbf} or \code{swd} object must be provided. For \code{sbf} object, type `"obs"', `"network"', `"field"'
are possible whereas all types are possible for \code{swd} object. Or specify `z' and `latlon'
without `sw'. 
}
\examples{
### Observations of year 1967
data(temperature)
names(temperature)

# Temperatures on 939 weather stations of year 1967    
temp67 <- temperature$obs[temperature$year == 1967] 
# Locations of 939 weather stations    
latlon <- temperature$latlon[temperature$year == 1967, ]

### Draw the temperature data
sw.plot(z=temp67, latlon=latlon, type="obs")

### Network design by BUD
data(netlab)
sw.plot(z=netlab, latlon=latlon, type="network")

### SBF representation of the observations
#eta <- c(0.961,0.923,0.852,0.723,0.506)
#out.pls <- sbf(obs=temp67, latlon=latlon, netlab=netlab, eta=eta,
#    method="pls", grid.size=c(50, 100), lambda=0.89)
# observation
#sw.plot(out.pls, type="obs")
# network design
#sw.plot(out.pls, type="network")
# field 
#sw.plot(out.pls, type="field")  

### Decomposition
#out.dpls <- swd(out.pls)
# observation
#sw.plot(out.dpls, type="obs")
# network design
#sw.plot(out.dpls, type="network")   
# SBF representation of the observations
#sw.plot(out.dpls, type="field")
# sw coefficient
#sw.plot(out.dpls, type="swcoeff")
# decomposition result
#sw.plot(out.dpls, type="decom")

# Thresholding  
#out.univ <- swthresh(out.dpls, policy="universal", by.level=TRUE, 
#    type="hard", nthresh=4)  
#par(oma=c(0,0,3.5,0))
#sw.plot(out.univ, type="decom") 
#mtext("Decomposition & Threshold", side = 3, outer = TRUE, 
#    cex = 1.2, line = 1)

# Reconstruction
#out.rec <- swr(out.univ)
#sw.plot(z=out.rec, type="recon", xlab="", ylab="")
}
\keyword{nonparametric}
