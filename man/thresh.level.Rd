\name{thresh.level}
\alias{thresh.level}
\title{Level-dependent Thresholding of SW Coefficients}
\description{
This function calculates level-dependent thresholded SW coefficients.
}
\usage{
thresh.level(x, norm, policy, Q, type) 
}
\arguments{
\item{x}{coefficients of multiscale SBF's}
\item{norm}{norm of multiscale SBF's (SW)}
\item{policy}{threshold technique. At present the possible policies are
    `"universal"', `"fdr"' and `"Lorentz"'.}
\item{Q}{parameter for the false discovery rate of `"fdr"' policy.}
\item{type}{the type of thresholding. This can be `"hard"', `"soft"' or `"Lorentz"'.}
}
\value{ 
\item{tgamma}{level-dependent thresholded SW coefficients}
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
\code{\link{mrs.comp.thresh.level}}
}
\keyword{nonparametric}
