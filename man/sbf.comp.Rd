\name{sbf.comp}
\alias{sbf.comp}
\title{Calculation of Field and Density with Multi-sale SBF's}
\description{
This function calculates field and density with Multi-sale SBF's.
}
\usage{
sbf.comp(point1, point2, site1, site2, coef, netlab, eta, p0)  
}
\arguments{
\item{point1}{latitude of extrapolation sites in radian}
\item{point2}{longitude of extrapolation sites in radian}
\item{site1}{latitude of observation sites in radian}
\item{site2}{longitude of observation sites in radian}
\item{coef}{coefficients of multi-scale SBF's}
\item{netlab}{vector of labels representing sub-networks} 
\item{eta}{bandwidth parameters for Poisson kernel} 
\item{p0}{specifies starting level for extrapolation. Among resolution levels \eqn{1, \ldots, L},  
resolution levels \eqn{p0+1, \ldots, L} will be included for extrapolation.} 
}
\details{
For a given field, this function provides a multiscale SBF representation
\deqn{
T(x)=\sum_{l}\sum_{j}\beta_{l,j}G_{l}(x \cdot x_j),
}{
T(x) = \sum_{l} \sum_{j} \beta_{l,j} G_{l}(<x, x_j>),
}
where \eqn{G_{l}(\cdot)}{G_{l}()} denotes a SBF with bandwidth \eqn{h_l} at multiresolution level \eqn{l}.
}
\value{
\item{aa}{a multiscale SBF field on observation's locations}
\item{bb}{density on observation's locations}
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
\code{\link{msbf.comp}}.
}
\keyword{nonparametric}
