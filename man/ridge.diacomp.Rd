\name{ridge.diacomp}
\alias{ridge.diacomp}
\title{Calculation of Generalized Cross-validation}
\description{
This function calculates generalized cross-validation for ridge regression.
}
\usage{
ridge.diacomp(out.ls, obs, lam) 
}
\arguments{
\item{out.ls}{an object of class `lsfit'}
\item{obs}{observations}
\item{lam}{smoothing parameter for penalized least squares method} 
}
\details{
This function calculates generalized cross-validation for ridge regression.
}
\value{
\item{rsq}{R-squared}
\item{gcv}{generalized cross-validation for ridge regression.}
\item{df}{degree of freedom}
}
\references{
Oh, H-S. (1999)  Spherical wavelets and their statistical analysis with applications to meteorological data. Ph.D. Thesis, 
Department of Statistics, Texas A\&M University, College Station.

Li, T-H. (1999) Multiscale representation and analysis of spherical data by spherical wavelets. 
\emph{SIAM Journal on Scientific Computing}, \bold{21}, 924--953.

Oh, H-S. and Li, T-H. (2004) Estimation of global temperature fields from scattered observations by 
a spherical-wavelet-based spatially adaptive method. \emph{Journal of the Royal Statistical Society
Ser.} B, \bold{66}, 221--238.
}
\seealso{
\code{\link{gcv.lambda}}, \code{\link{ridge.comp}}.
}
\keyword{nonparametric}
