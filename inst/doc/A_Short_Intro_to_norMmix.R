## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  message = FALSE,
  warning = FALSE,
  error = FALSE,
  tidy = FALSE,
  cache = FALSE
)

## ----setup--------------------------------------------------------------------
library(norMmix)
set.seed(2020)

## ----dev-and-compile, echo = FALSE, eval = FALSE------------------------------
#  devtools::load_all()
#  rmarkdown::render("~/R/Pkgs/norMmix/vignettes/A_Short_Intro_to_norMmix.Rmd")

## -----------------------------------------------------------------------------
faith <- norMmixMLE(faithful, 3, model="VVV", initFUN=claraInit)

## -----------------------------------------------------------------------------
plot(faith)

## ----norMmix------------------------------------------------------------------
w <- c(0.5, 0.3, 0.2)
mu <- matrix(1:6, 2, 3)
sig <- array(c(2,1,1,2,
               3,2,2,3,
               4,3,3,4), c(2,2,3))
nm <- norMmix(mu, Sigma=sig, weight=w)
plot(nm)

## ----norMmix_panels-----------------------------------------------------------
plot(MW32)

## ----norMmix_data-------------------------------------------------------------
x <- rnorMmix(500, nm)
plot(nm, xlim = c(-5,10), ylim = c(-5, 12),
     main = "500 observations from a mixture of 3 bivariate Gaussians")
points(x)

## -----------------------------------------------------------------------------
ret <- norMmixMLE(x, 3, model="VVV", initFUN=claraInit)
ret # -> print.norMmixMLE(ret)

## ---- plot-MLE----------------------------------------------------------------
plot(ret)

## -----------------------------------------------------------------------------
# suppose we wanted some mixture model, let
mu <- matrix(1:6, 2,3)      # 2x3 matrix -> 3 means of dimension 2
w <- c(0.5, 0.3, 0.2)       # needs to sum to 1
diags <- c(4, 3, 5)         # these will be the entries of the diagonal of the covariance matrices (see below)

nm <- norMmix(mu, Sigma=diags, weight=w)
print(nm)
str(nm)

## -----------------------------------------------------------------------------
nm$mu
nm$Sigma

