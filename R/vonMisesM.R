#' @title Random number generation from von Mises mixture distribution
#' @description \code{rvonMisesM} Generate random numbers from von Mises mixture distribution
#'
#' @importFrom CircStats rvm
#' @param n number of observations.
#' @param mu mean direction of each component distribution.
#' @param kappa non-negative numeric value for the concentration parameter of each component distribution.
#' @param theta simplex vector for mixing ratios.
#' @return Random number of length n.
#' @export
#' @examples
#' #rvonMisesM(10, c(2.1, 3.4), c(1.0, 1.0), c(0.5, 0.5))

rvonMisesM <- function(n, mu, kappa, theta) {
  y <- rep(NA, n)
  K <- length(mu)
  for(i in 1:n) {
    g <- sample(1:K, 1, prob = theta)
    y[i] <- CircStats::rvm(1, mu[g], kappa[g])
  }

  y
}

#' @title Density function for the von Mises mixture distribution
#' @description \code{dvonMisesM} Density function for the von Mises mixture distribution
#'
#' @importFrom CircStats dvm
#' @param x random variables.
#' @param mu mean direction of each component distribution.
#' @param kappa non-negative numeric value for the concentration parameter of each component distribution.
#' @param theta simplex vector for mixing ratios
#' @return Probability density of the same length as object x
#' @export
#' @examples
#' #rvonMisesM(pi, c(2.1, 3.4), c(1.0, 1.0), c(0.5, 0.5))
dvonMisesM<- function(x, mu, kappa, theta) {
  N <- length(x)
  pdfs <- rep(NA, N)
  for(i in 1:N) {
    pdf <- CircStats::dvm(x[i], mu, kappa)
    pdfs[i] <- sum(pdf * theta)
  }

  pdfs
}
