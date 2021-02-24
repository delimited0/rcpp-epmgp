d <- 1000
lb <- rep(0, d)
ub <- rep(1, d)
mu <- rep(0, d)
Sigma <- diag(d)

result <- epmgp::axisepmgp(mu, Sigma, lb, ub)
result$logZ

# Eric example ----
load("testing/patrick_v1.RData")
result <- epmgp::axisepmgp(mu_beta, Q_beta_inv, lb, ub)
epmgp::pmvn(lb, ub, mu_beta, Q_beta_inv)
result_TN <- 
result_poly <- epmgp::epmgp(mu_beta, Q_beta_inv, diag(length(mu_beta)), lb, ub, 200)

result_poly$mu
result_poly$Sigma[1:5, 1:5]

result$Sigma[1:5, 1:5]


load("testing/patrick_v2.RData")
result <- epmgp::axisepmgp(mu_beta, Q_beta_inv, lb, ub)
TruncatedNormal::pmvnorm(mu_beta, Q_beta_inv, lb, ub)

# trunc norm moments Eric example nans:
dat <- read.csv("testing/trunc_norm_moments_failure.csv")

result <- epmgp::trunc_norm_moments(lb, ub, dat$nu_cavity / dat$tau_cavity, 1 / dat$tau_cavity)

result <- epmgp::trunc_norm_moments(lb[1], ub[1], 
                                    dat$nu_cavity[1] / dat$tau_cavity[1], 
                                    1 / dat$tau_cavity[1])

result <- epmgp::trunc_norm_moments(lb[3], ub[3], 
                                    dat$nu_cavity[3] / dat$tau_cavity[3], 
                                    1 / dat$tau_cavity[3])

# miwa craig ----
d < 2
lb <- rep(0, d)
ub <- rep(Inf, d)
mu <- rep(0, d)
Sigma <- .5 * diag(d) + .5 * rep(1, d) %*% t(rep(1, d))
epmgp::axisepmgp(mu, Sigma, lb, ub)
