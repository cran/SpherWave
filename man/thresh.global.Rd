\name{thresh.global}
\alias{thresh.global}
\title{Global Thresholding of SW Coefficients}
\description{
This function calculates global thresholded SW coefficients.
}
\usage{
thresh.global(x, lam, type)
}
\arguments{
\item{x}{coefficients of multiscale SBF's}
\item{lam}{thresholding value}
\item{type}{the type of thresholding. This can be `"hard"', `"soft"' or `"Lorentz"'.}
}
\value{ 
\item{tgamma}{global thresholded SW coefficients}
}
\references{
Donoho, D.~L. and Johnstone, I.~M. (1994) Ideal spatial
adaptation by wavelet shrinkage. \emph{Biometrika}, \bold{81},
425--455.

Oh, H-S. (1999)  Spherical wavelets and their statistical analysis with applications to meteorological data. Ph.D. Thesis, 
Department of Statistics, Texas A\&M University, College Station.

Oh, H-S. and Li, T-H. (2004) Estimation of global temperature fields from scattered observations by 
a spherical-wavelet-based spatially adaptive method. \emph{Journal of the Royal Statistical Society
Ser.} B, \bold{66}, 221--238.
}
\seealso{
\code{\link{mrs.comp.thresh.global}}.
}
\keyword{nonparametric}
