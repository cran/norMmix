\name{nMm2par}
\alias{nMm2par}
\title{Multivariate Normal Mixture Model to parameter for MLE}
\description{
  From a \code{"\link{norMmix}"}(-like) object, return the numeric
  parameter vector in our MLE parametrization.
}
\usage{
  nMm2par(obj,
        model = c("EII", "VII", "EEI", "VEI", "EVI",
                  "VVI", "EEE", "VEE", "EVV", "VVV"),
        meanFUN = mean.default,
        checkX = FALSE)
}
\arguments{
  \item{obj}{a \code{\link{list}} containing \describe{
      \item{\code{sig}:}{covariance matrix array,}
      \item{\code{mu}:}{mean vector matrix,}
      \item{\code{w}:}{= weights,}
      \item{\code{k}:}{= number of components,}
      \item{\code{p}:}{= dimension}
    }
  }
  \item{model}{a \code{\link{character}} string specifying the (Sigma)
    model, one of those listed above.}
  \item{meanFUN}{a \code{\link{function}} to compute a mean (of variances typically).}
  \item{checkX}{a boolean. check for positive definiteness of covariance matrix.}
}
\details{
    This transformation forms a vector from the parameters of a normal
    mixture. These consist of weights, means and covariance matrices.

    Covariance matrices are given as D and L from the LDLt decomposition
}
\value{
  vector containing encoded parameters of the mixture.
  first, the centered log ratio of the weights, then the means, and then the model specific encoding
  of the covariances.
}
\seealso{
    the \emph{inverse} function of \code{\link{nMm2par}()} is \code{\link{par2nMm}()}.
}

\examples{
A <- MW24
nMm2par(A, model = A$model)
# [1] -0.3465736  0.0000000  0.0000000  0.0000000  0.0000000  0.0000000
# [7] -2.3025851

## All MW* models in {norMmix} pkg:
pkg <- "package:norMmix"
lMW <- mget(ls(pattern = "^MW", pkg), envir=as.environment(pkg))
lM.par <- lapply(lMW, nMm2par)
## but these *do* differ  ___ FIXME __ 
modMW <- vapply(lMW, `[[`, "model", FUN.VALUE = "XYZ")
cbind(modMW, lengths(lM.par),  npar = sapply(lMW, npar))[order(modMW),]
}
