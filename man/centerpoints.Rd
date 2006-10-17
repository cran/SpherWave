\name{centerpoints}
\alias{centerpoints}
\alias{modregcenter}
\alias{modredcenter}
\alias{gotregcenter}
\alias{gotredcenter}
\alias{regcenter}
\alias{redcenter}
\title{Generation of Center Points}
\description{
Given a multi-resolution level, generate center points of each grid box by G\"{o}ttlemann's method, 
modified G\"{o}ttlemann's method or standard method. 
}
\usage{
modregcenter(l) # for modified G\"{o}ttlemann's regular grid
modredcenter(l) # for modified G\"{o}ttlemann's reduced grid
gotregcenter(l) # for G\"{o}ttlemann's regular grid
gotredcenter(l) # for G\"{o}ttlemann's reduced grid
regcenter(l)    # for standard regular grid
redcenter(l)    # for standard reduced grid
}
\arguments{
\item{l}{an integer which denotes the index of multi-resolution level}
}
\details{
It is for obtaining the center points of each grid box from the grid according to multi-resolution level. The values will be used for 
network design.
}
\value{
\item{center}{center points from grid}
}
\references{
G\"{o}ttlemann, J. (1996) Locally supported wavelets on the sphere. Preprint, Johannes Gutenberg University, Mainz.  

Oh, H-S. (1999)  Spherical wavelets and their statistical analysis with applications to meteorological data. Ph.D. Thesis, 
Department of Statistics, Texas A\&M University, College Station.
}
\seealso{
\code{\link{network}}
}
\examples{
modregcenter(3)$center
modredcenter(3)$center
gotregcenter(3)$center
gotredcenter(3)$center
regcenter(3)$center
redcenter(3)$center
}
\keyword{nonparametric}
