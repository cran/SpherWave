\name{temperature}
\docType{data}
\alias{temperature}
\title{The Surface Air Temperature}
\description{
the surface air temperature in Celsius observed by a network of weather stations in 1961 through 1990.
}
\usage{
data(temperature)
}
\format{
A list of year, latlon (global grid point in degree), obs (temperature)
}
\source{This data set was organized by Jones, Raper, Cherry, Goodess, Wigley, Santer, and Kelly (1991). 
The primary sources of this
data can be obtained from \url{http://cdiac.esd.ornl.gov/ftp}.
}
\examples{
### Observations of year 1967
data(temperature)
names(temperature)

# Temperatures on 939 weather stations of year 1967    
temp67 <- temperature$obs[temperature$year == 1967] 
# Locations of 939 weather stations    
latlon <- temperature$latlon[temperature$year == 1967, ]

# Plot of the observations
sw.plot(z = temp67, latlon=latlon, type="obs", xlab="", ylab="")
}
\keyword{datasets}
