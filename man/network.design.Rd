\name{network.design}
\alias{network.design}
\title{Network Design}
\description{
produces multi-resolution network.
}
\usage{
network.design(latlon, method = "Oh", type = "reduce", nlevel, x) 
}
\arguments{
\item{latlon}{grid points of observation sites in degree}
\item{method}{network design method, "cover", "ModifyGottlemann", "Gottlemann" or "Oh"}
\item{type}{specifies grid type, "regular" or "reduced"}
\item{nlevel}{the number of observations in each resolution when using the method `"cover"'}
\item{x}{radius in degree}
}
\details{
This function partitions the grid points of observations into networks.
Each network corresponds resolution level and level 1 is the most detailed level.
Possible methods are

    `"cover"' for utilizing `"cover.design"' in the package `"field"'
    
    `"ModifyGottlemann"' for modified G\"{o}ttlemann's method
    
    `"Gottlemann"' for G\"{o}ttlemann's method
    
    `"Oh"' for Oh's method.
    
For `"ModifyGottlemann"', `"Gottlemann"' and `"Oh"', two types of design, `"regular"' and `"reduced"' are provided.
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

#netlab <- network.design(latlon=latlon, method="cover",
#    nlevel=c(507, 244, 117, 49, 22))
netlab <- network.design(latlon=latlon, method="ModifyGottlemann", 
    type="regular", x=3)
netlab <- network.design(latlon=latlon, method="Gottlemann", 
    type="regular", x=2)
netlab <- network.design(latlon=latlon, method="Oh", 
    type="reduce", x=5)
}
\keyword{nonparametric}
