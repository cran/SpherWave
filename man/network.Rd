\name{network}
\alias{network}
\alias{reg.grid}
\alias{red.grid}
\alias{gotreg.grid}
\alias{gotred.grid}
\alias{hsreg.grid}
\alias{hsred.grid}
\title{Generation of Network Levels}
\description{
produces multi-resolution network by G\"{o}ttlemann's method, modified G\"{o}ttlemann's method
or standard method. 
}
\usage{
reg.grid(x, latlon)    # for modified G\"{o}ttlemann's regular grid
red.grid(x, latlon)    # for modified G\"{o}ttlemann's reduced grid
gotreg.grid(x, latlon) # for G\"{o}ttlemann's regular grid
gotred.grid(x, latlon) # for G\"{o}ttlemann's reduced grid
hsreg.grid(x, latlon)  # for standard regular grid
hsred.grid(x, latlon)  # for standard reduced grid
}
\arguments{
\item{x}{radius of territory in degree}
\item{latlon}{grid points of observation sites in degree}
}
\details{
This function partitions the grid points of observations into networks.
Each network corresponds resolution level and level 1 is the most detailed level.
}
\value{
\item{netlab}{vector of network labels indicating level of multi-resolution.}
}
\references{
Oh, H-S. (1999)  Spherical wavelets and their statistical analysis with applications to meteorological data. Ph.D. Thesis, 
Department of Statistics, Texas A\&M University, College Station.
}
\examples{
### Observations of year 1967
data(temperature)
names(temperature)

# Locations of 939 weather stations    
latlon <- temperature$latlon[temperature$year == 1967, ]

netlab <- reg.grid(x=3, latlon)
netlab <- red.grid(x=3, latlon)
netlab <- gotreg.grid(x=2, latlon)
netlab <- gotred.grid(x=2, latlon)
netlab <- hsreg.grid(x=5, latlon)
netlab <- hsred.grid(x=5, latlon)
}
\keyword{nonparametric}
