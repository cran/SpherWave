\name{coefmatrix}
\alias{coefmatrix}
\title{Computation of Coefficients of SBF and SW}
\description{
This function generates several coefficients such as coefficients of SBF in spherical wavelets (SW), 
coefficients of SBF after removing subnet \eqn{l}, and coefficients of SW for subnet \eqn{l}.
}
\usage{
coefmatrix(beta1, fullcov, netlab, l)
}
\arguments{
\item{beta1}{coefficients of SBF from previous SBF representation}
\item{fullcov}{covariance matrix of all observation sites}
\item{netlab}{vector of labels representing sub-networks} 
\item{l}{resolution level} 
}
\details{
The multiresolution analysis based on SBF is derived from the problem of characterizing the loss in an SBF representation 
as the number of observations are more larger. This function provides the coefficients of basis functions of multiresolution levels. 
For details, see references below.
}
\value{ 
\item{wcoef}{coefficients of SBF in SW}
\item{beta2}{coefficients of SBF after removing sub-network \eqn{l}}
\item{gamma1}{coefficients of SW for sub-network \eqn{l}}
\item{alpha1}{detailed coefficients of SBF for sub-network \eqn{l}}
\item{norm}{norms of SW for sub-network \eqn{l}}
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
\code{\link{mracoef.comp}}
}
\keyword{nonparametric}
