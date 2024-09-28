# log likelihood ---------------------------------------------------------------
log_likelihood <- function(y, mu, kappa, theta) {
  N <- length(y)
  tmp <- rep(NA, N)
  for (n in 1:N) {
    tmp[n] <- log(sum(dvm(y[n], mu, kappa) * theta))
  }
  sum(tmp)
}

# E step -----------------------------------------------------------------------
responsibility_calc <- function(y, mu, ka, th) {
  th * dvm(y, mu, ka)/sum(th * dvm(y, mu, ka))
}

E_step <- function(N, K, y, mu, kappa, theta) {
  r <- array(NA, dim = c(N, K))
  for(i in 1:N) {
    r[i,] <- responsibility_calc(y[i], mu, kappa, theta)
  }

  r
}

# M step -----------------------------------------------------------------------
mu_est <- function(C, S) {
  if((C > 0) && (S >= 0)) {
    atan(S/C)
  }else if((C == 0) && (S > 0)) {
    pi/2.0
  }else if(C < 0) {
    atan(S/C) + pi
  }else if((C == 0) && (S < 0)) {
    (3*pi)/2
  }else if((C > 0) && (S < 0)) {
    atan(S/C) + 2 * pi
  }else{
    print("S or C values are abnormal ...")
  }
}

M_step <- function(y, r, K) {
  N_k <- apply(r, 2, sum)

  # Initialize parameters
  mu_new <- rep(NA, K)
  ka_new <- rep(NA, K)
  th_new <- rep(NA, K)

  for(k in 1:K) {
    # update mu ---
    cos_sum   <- sum(cos(y) * r[,k]) / N_k[k]
    sin_sum   <- sum(sin(y) * r[,k]) / N_k[k]
    mu_new[k] <- mu_est(cos_sum, sin_sum)# + 2*pi

    # update kappa ---
    d <- sqrt(cos_sum^2 + sin_sum^2)

    if(d < 0.53) {
      ka_new[k] <- (2.0 * d) + (d^3) + (5*(d^5))/6
    }else if((0.53 <= d) && (d < 0.85)) {
      ka_new[k] <- -0.4 + 1.39 * d + 0.43/(1 - d)
    }else if(d >= 0.85) {
      ka_new[k] <- 1/((d^3) - 4*(d^2) + 3*d)
    }else{
      stop("d <- sqrt(C^2 + S^2) values are abnormal ...")
    }
  }
  # update theta ---
  th_new <- N_k/sum(N_k)

  list(mu_new = mu_new, ka_new = ka_new, th_new = th_new)
}

# EM algorithm -----------------------------------------------------------------
#' @title Parameter estimation of von Mises Mixture model using EM algorithm
#' @description \code{fit.EM} XXXX
#' @importFrom stats rgamma
#' @importFrom stats runif
#' @param N_step XXXX.
#' @param K XXX.
#' @param y XXX.
#' @return XXX.
#' @export
#' @examples
#' #fit.EM(10, K, y)

fit.EM <- function(N_step, K, y) {
  N <- length(y)

  # Initialize parameters
  mu    <- runif(K, 0.0, 2*pi)
  kappa <- runif(K, 0.0, 10.0)
  theta <- rgamma(K, 100, 10)
  theta <- theta/sum(theta)

  # log-likelihood
  ll <- rep(NA, N_step)

  for(t in 1:N_step) {
    # E step
    r <- E_step(N, K, y, mu, kappa, theta)

    # M step
    pars  <- M_step(y, r, K)
    mu    <- pars$mu_new
    kappa <- pars$ka_new
    theta <- pars$th_new

    # log-likelihood
    ll[t] <- log_likelihood(y, mu, kappa, theta)
  }

  return(list(pars = pars, log_likelihood = ll))
}

