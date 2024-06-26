
\name{dnorMmix}
\alias{dnorMmix}
\title{Density from Multivariate Normal Mixture Distribution}
\description{
  Calculates the probability density function of the multivariate normal distribution.
}
\usage{
dnorMmix(x, nMm)
}
\arguments{
  \item{x}{a vector or matrix of multivariate observations}
  \item{nMm}{a \code{"\link{norMmix}"} object}
}
%% \details{
%% }
\value{
  Returns the density of \code{nMm} at point \code{x}. Iterates over components of the
  mixture and returns weighted sum of \code{\link[mvtnorm]{dmvnorm}}.
}
%% \references{
%% }
\author{Nicolas Trutmann}
\seealso{
  \code{\link{rnorMmix}}
}
%% \examples{
%% }
