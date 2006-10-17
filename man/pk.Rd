\name{pk}
\alias{pk}
\title{Calculation of Normalized Poisson Kernel}
\description{
This function calculates normalized Poisson kernel as a function of angle.
}
\usage{
pk(theta, eta)
}
\arguments{
\item{theta}{angle}
\item{eta}{bandwidth parameter for Poisson kernel}
}
\details{
This function calculates normalized Poisson kernel as a function of angle. 
}
\value{
vector of normalized Poisson kernel.
}
\examples{
theta <- c(-100:100)/100

par(mfrow=c(1,1), pty="m", mar=c(4,4,7,1)+0.1)
plot(theta, pk(theta, 0.9), type="l", xlab="angle (x pi)", ylab="SBF", 
    xlim=c(-1, 1), ylim=c(0, 1), lab=c(9, 7, 7), cex=1)
lines(theta, pk(theta, 0.7), lty=2)
lines(theta, pk(theta, 0.5), lty=4)
legend(0.4, 0.8, legend = c("ETA = 0.9","ETA = 0.7","ETA = 0.5"), 
    lty=c(1,2,4), cex=0.7)
}   
\keyword{nonparametric}
