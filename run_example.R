# setting for parameters -------------------------------------------------------
mu_true    <- c(2.0, 4.0)
theta_true <- c(0.5, 0.5)
kappa_true <- c(3.2, 4.3)

# data -------------------------------------------------------------------------
N <- 2000
y <- rvonMisesM(N, mu_true, kappa_true, theta_true)

# EM algorithm -----------------------------------------------------------------
K      <- 2
N_step <- 20
res    <- fit.EM(N_step, K, y)

# results ----------------------------------------------------------------------
plot(res$log_likelihood, type = "l")

mu_new <- res$pars$mu_new
ka_new <- res$pars$ka_new
th_new <- res$pars$th_new

hist(y, xlim = c(0, 2*pi), ylim = c(0, 0.8), probability = T, main = "")
curve(dvonMisesM(x, mu_true, kappa_true, theta_true), from = 0, to = 2*pi, add = T, col = "blue")
curve(dvonMisesM(x, mu_new , ka_new    , th_new    ), from = 0, to = 2*pi, add = T, col = "red")
