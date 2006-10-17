\name{swr}
\alias{swr}
\title{Spherical Wavelet Reconstruction of `swd' Object} 
\description{
This function performs spherical wavelet reconstruction.
}
\usage{
swr(swd) 
}
\arguments{
\item{swd}{an object of class `swd'}
}
\details{
This function performs spherical wavelet reconstruction.
}
\value{ 
\item{recon}{the spherical wavelet reconstruction}
}
\references{
Oh, H-S. and Li, T-H. (2004) Estimation of global temperature fields from scattered observations by 
a spherical-wavelet-based spatially adaptive method. \emph{Journal of the Royal Statistical Society
Ser.} B, \bold{66}, 221--238.
}
\seealso{
\code{\link{sbf}}, \code{\link{swd}}, \code{\link{swthresh}}
}
\examples{
### Observations of year 1967
data(temperature)
names(temperature)

# Temperatures on 939 weather stations of year 1967    
temp67 <- temperature$obs[temperature$year == 1967] 
# Locations of 939 weather stations    
latlon <- temperature$latlon[temperature$year == 1967, ]

### Network design by BUD
data(netlab)

### Bandwidth for Poisson kernel
eta <- c(0.961, 0.923, 0.852, 0.723, 0.506)

### SBF representation of the observations by pls
out.pls <- sbf(obs=temp67, latlon=latlon, netlab=netlab, eta=eta, 
    method="pls", grid.size=c(50, 100), lambda=0.89)

### Decomposition
out.dpls <- swd(out.pls)

### Thresholding
out.univ <- swthresh(out.dpls, policy="universal", by.level=TRUE, 
    type="hard", nthresh=4)  

### Reconstruction
out.rec <- swr(out.univ)
sw.plot(z=out.rec, type="recon", xlab="", ylab="")
}
\keyword{nonparametric}
