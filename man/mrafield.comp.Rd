\name{mrafield.comp}
\alias{mrafield.comp}
\title{Decomposition of a Field}
\description{
This function computes global and local components (fields) on grid from an initial field. 
}
\usage{
mrafield.comp(grid, coeff, site, netlab, eta, field, density) 
}
\arguments{
\item{grid}{grid points of extrapolation sites in radian}
\item{coeff}{coefficients of multi-scale SBF's}
\item{site}{grid points of observation sites in radian}
\item{netlab}{vector of labels representing sub-networks} 
\item{eta}{bandwidth parameters for Poisson kernel} 
\item{field}{extrapolation on grid}
\item{density}{density of locations induced from an initial field}
}
\details{
This function generates decomposition of a field,
\deqn{
T_1(x)=T_l(x)+D_{l-1}(x)+\ldots+D_1(x), \qquad l=2,\ldots,L
}{
T_1(x) = T_l(x) + D_{l-1}(x) + \ldots + D_1(x), l = 2,\ldots, L
}
where a global component \eqn{T_{l+1}(x) \in \mathcal{V}_{l+1}}{T_{l+1}(x) \in V_{l+1}}
and a local component \eqn{D_l(x) \in \mathcal{W}_l}{D_l(x) \in W_l}. The 
corresponding space are nested as \eqn{\mathcal{V}_l \supset \mathcal{V}_{l+1}}{V_l \supset V_{l+1}}, 
so that \eqn{\mathcal{V}_l=\mathcal{V}_{l+1} \oplus \mathcal{W}_l}{V_l = V_{l+1} + W_l}. 
}
\value{ 
\item{global}{matrix of successively smoothed data}
\item{local}{matrix of difference of successively smoothed data}
\item{density}{density of locations in global and local fields}
\item{swcoeff}{spherical wavelet coefficients}
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
\code{\link{sbf}}, \code{\link{swd}}, \code{\link{swthresh}}, \code{\link{swr}}
}
\keyword{nonparametric}
