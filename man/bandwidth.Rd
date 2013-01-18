\name{bandwidth}
\alias{bandwidth}
\title{Calculation of Bandwidth Parameter}
\description{
Given the number of observations at a level, this function calculates bandwidth parameter of spherical basis function (SBF). 
}
\usage{
bandwidth(n) 
}
\arguments{
\item{n}{the number of observations at a level} 
}
\details{
This function is used for obtaining the bandwidth of the coarsest network level \eqn{L}, \eqn{h_L}. 
From geometry, the surface area covered by surface mass distribution with variance \eqn{\sigma^2} over 
unit sphere \eqn{\Omega} is \eqn{2\pi(1-\sqrt{1-\sigma^2})}{2 \pi (1 - \sqrt{1 - \sigma^2})}. 
Since the total surface area of the unit sphere is \eqn{4\pi}{4 \pi} and 
the variance of SBF induced from Poisson kernel is 
\eqn{\sigma^2=\left(\frac{1-h^2}{1+h^2}\right)^2}{\sigma^2=\left((1 - h^2)/(1 + h^2)\right)^2}, 
the surface covered are is \eqn{2\pi\left(1-\sqrt{1-(\frac{1-h^2}{1+h^2})^2}\right)}{2 \pi \left(1 - \sqrt{1 - ((1 - h^2)/(1 + h^2))^2}\right)}. 
Under the assumption that the observations are 
distributed equally over the sphere, it can be easily known how many observation are needed in order to cover 
the whole sphere with fixed \eqn{h}, and how large the \eqn{h} is needed to cover the sphere when the number 
of observations are fixed as follows:
\deqn{
\# \mbox{ of observations}=n=\frac{2}{1-\sqrt{1-(\frac{1-h^2}{1+h^2})^2}} \mbox{ and } h=\sqrt{\frac{1-a_n}{1+a_n}}
}{
\# of observations = n = 2/(1 - \sqrt{1 - ((1 - h^2)/(1 + h^2))^2}) and h = \sqrt{(1 - a_n)/(1 + a_n)}
}
where \eqn{a_n=\sqrt{1-(1-\frac{2}{n})^2}}{a_n = \sqrt{1 - (1 - 2/n)^2}}.
}
\value{
\item{h}{bandwidth parameter at a level} 
}
\references{
Oh, H-S. (1999)  Spherical wavelets and their statistical analysis with applications to meteorological data. Ph.D. Thesis, 
Department of Statistics, Texas A\&M University, College Station.
}
\seealso{
\code{\link{eta.comp}}
}
\examples{
bandwidth(20)
}
\keyword{nonparametric}
