\name{gcv.lambda}
\alias{gcv.lambda}
\title{Calculation of Generalized Cross-validation}
\description{
This function calculates generalized cross-validation for ridge regression.
}
\usage{
gcv.lambda(obs, latlon, netlab, eta, approx=FALSE, lambda) 
}
\arguments{
\item{obs}{observations}
\item{latlon}{grid points of observation sites in degree}
\item{netlab}{vector of labels representing sub-networks} 
\item{eta}{bandwidth parameters for Poisson kernel} 
\item{approx}{if TRUE, approximation is used.}
\item{lambda}{smoothing parameter for penalized least squares method} 
}
\value{
\item{gcv}{generalized cross-validation for ridge regression.}
}
\seealso{
\code{\link{ridge.diacomp}}, \code{\link{ridge.comp}}.
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

### Select smoothing parameter lambda by generalized cross-validation
lam <- seq(0.1, 0.9, ,9)
gcv <- NULL
for(i in 1:length(lam))
    gcv <- c(gcv, gcv.lambda(obs=temp67, latlon=latlon, 
        netlab=netlab, eta=eta, lambda=lam[i])$gcv)
lam[gcv == min(gcv)]
}
\keyword{nonparametric}
