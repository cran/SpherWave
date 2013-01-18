\name{mesh}
\alias{mesh}
\title{Creation of M by N grid}
\description{
This function creates M by N longitude-latitude grid.
}
\usage{
mesh(M, N)
}
\arguments{
\item{M}{M longitudes of M by N grid}
\item{N}{N latitudes of M by N grid}
}
\details{
This function creates M by N longitude-latitude grid.
}
\value{
\item{theta}{latitudes in radian of M by N longitude-latitude grid}
\item{phi}{longitudes in radian of M by N longitude-latitude grid}
}
\examples{
mesh(100, 50)
}
\keyword{nonparametric}
