\name{eta.comp}
\alias{eta.comp}
\title{Calculation of Bandwidth Parameters of Poisson Kernel/SBF's}
\description{
This function calculates bandwidth parameters of SBF's induced by Poisson Kernel.
}
\usage{
eta.comp(netlab) 
}
\arguments{
\item{netlab}{the index vector of network level} 
}
\details{
The bandwidths \eqn{h_l} of SBF's are chosen as follows
\deqn{
h_{L-1}=e^{-\rho_l},\quad l=1,2,\ldots,L-1, 
}{
h_{L-1} = e^{-\rho_l}, l = 1, 2,\ldots, L-1, 
} 
where \eqn{\rho_l=\rho^{*}/2^l}{\rho_l = \rho^{*} / 2^l}. The \eqn{\rho^*} is obtained from the bandwidth of the
coarsest network level \eqn{L}, \eqn{h_L}, that is 
\deqn{
\rho^*=-\log h_L.
}{
\rho^* = -\log h_L.
} 
Note that \eqn{h_L} is obtained from the function bandwidth.
}
\value{
\item{eta}{bandwidth parameters for Poisson kernel/SBF's} 
}
\references{
Oh, H-S. (1999)  Spherical wavelets and their statistical analysis with applications to meteorological data. Ph.D. Thesis, 
Department of Statistics, Texas A\&M University, College Station.
}
\seealso{
\code{\link{bandwidth}}
}
\examples{
eta.comp(5)
}
\keyword{nonparametric}
