\name{grid}
\alias{grid}
\alias{modreggrid}
\alias{modredgrid}
\alias{gotreggrid}
\alias{gotredgrid}
\alias{reggrid}
\alias{redgrid}
\title{Generation of Grid}
\description{
produces grid points by G\"{o}ttlemann's method, modified G\"{o}ttlemann's method
or standard method. 
}
\usage{
modreggrid(l)    # for modified G\"{o}ttlemann's regular grid
modredgrid(l)    # for modified G\"{o}ttlemann's reduced grid
gotreggrid(l)    # for G\"{o}ttlemann's regular grid
gotredgrid(l)    # for G\"{o}ttlemann's reduced grid
reggrid(l)       # for standard regular grid
redgrid(l)       # for standard reduced grid
}
\arguments{
\item{l}{levels}
}
\details{
This function generates the grid points on globe.
}
\value{
\item{grid}{grid points}
}
\references{
Oh, H-S. (1999)  Spherical wavelets and their statistical analysis with applications to meteorological data. Ph.D. Thesis, 
Department of Statistics, Texas A\&M University, College Station.
}
\examples{
rec.reg <- reggrid(3)$grid 
rec.red <- redgrid(3)$grid 
mod.reg <- modreggrid(3)$grid
mod.red <- modredgrid(3)$grid
got.reg <- gotreggrid(2)$grid
got.red <- gotredgrid(2)$grid

par(mfrow = c(3, 2), mar = c(2.1, 4.1, 4.1, 1.1))
world() 
points(rec.reg[, 2], rec.reg[, 1], cex = 0.7)
title(main = "(a)")
world() 
points(rec.red[, 2], rec.red[, 1], cex = 0.7)
title(main="(b)")
world() 
points(mod.reg[, 2], mod.reg[, 1], cex = 0.7)
title(main = "(c)")
world() 
points(mod.red[, 2], mod.red[, 1], cex = 0.7)
title(main = "(d)")
world() 
points(got.reg[, 2], got.reg[, 1], cex = 0.7)
title(main = "(e)")
world() 
points(got.red[, 2], got.red[, 1], cex = 0.7)
title(main = "(f)")
}
\keyword{nonparametric}
