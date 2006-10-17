\name{multi-levels}
\alias{multi-levels}
\alias{modnetlevel}
\alias{gotnetlevel}
\alias{netlevel}
\title{Deciding the Number of Multi-Resolution Levels}
\description{
This function decides the number of multi-resolution levels by G\"{o}ttlemann's method, modified G\"{o}ttlemann's method
or standard method.}
}
\usage{
modnetlevel(angle) # for modified G\"{o}ttlemann's grid
gotnetlevel(angle) # for G\"{o}ttlemann's grid
netlevel(angle)    # for standard grid
}
\arguments{
\item{angle}{radius (geodesic distance) from locations of data within a territory to the center point of the territory}
}
\value{
\item{nlevel}{the number of multi-resolution levels.}
}
\references{
Oh, H-S. (1999)  Spherical wavelets and their statistical analysis with applications to meteorological data. Ph.D. Thesis, 
Department of Statistics, Texas A\&M University, College Station.
}
\examples{
modnetlevel(3 * pi/180) # for modified G\"{o}ttlemann's grid
gotnetlevel(2 * pi/180) # for G\"{o}ttlemann's grid
netlevel(5 * pi/180)    # for standard grid
}
\keyword{nonparametric}
