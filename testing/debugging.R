d <- 200
# mu <- rep(0, d)
Sigma <- .5 * diag(d) + .5 * rep(1, d) %*% t(rep(1, d))
# Sigma <- diag(d)
# Sigma <- 2 * diag(d) - .5 * matrix(1, d, d)
lb <- rep(0, d)
ub <- rep(Inf, d)
n <- 1000
initial <- rep(1.1, d)
# lb <- rep(-1, d)
# ub <- rep(2, d)

samples <- epmgp::rtmvn(n, mu, Sigma, lb, ub, initial = initial)

ess_samples <- lincongauss::rtmvn(n, mu, Sigma, lb, ub, x_init = initial)

vaidya_samples <- gausspolytopewalk::rtmvn(1000, mu, Sigma, lb, ub, "vaidya",
                                           r = .5, lazy = 0, initial = rep(1.1, d))

plot(samples[, 1:2], pch = 20)
plot(ess_samples[, 1:2])
plot(vaidya_samples[, 1:2])

epmgp::pmvn(lb, ub, mu, Sigma)

moments <- epmgp::moments(lb, ub, rep(0, d), diag(d))

constr <- get_polytope_constraints(lb, ub, mu, t(chol(Sigma)))
n <- 10000
samples <- epmgp::sample_epess(n, moments$mu, chol(moments$Sigma), 
             constr$A, constr$b, J = 1, N = 1, initial = rep(.1, d))
epmh_samples <- epmgp::sample_epmh(n, moments$mu, chol(moments$Sigma),
                                   A, b, rep(.1, d))

plot(t(epmh_samples)[, 1:2])
plot(t(samples)[, 1:2])

botev_samples <- TruncatedNormal::rtmvnorm(n, mu, Sigma, lb, ub)



plot(botev_samples[, 1:2])


# test pseudo likelihood ----
d <- 2
mu <- rep(0, d)
Sigma <- .5 * diag(d) + .5 * rep(1, d) %*% t(rep(1, d))
# Sigma <- diag(d)
lb <- rep(1, d)
ub <- rep(2, d)

ep_moments <- epmgp::moments(lb, ub, mu, Sigma)
A <- rbind(diag(d), -diag(d))
b <- c(-lb, ub)
ep_chol <- chol(ep_moments$Sigma)
x = c(-.21591, -.1714)
  
epmgp::pseudo_llik(c(.5, .5), ep_moments$mu, ep_chol, A, b, 1, 1)


y = c(1.5, 1.5)
llik(y, ep_moments$mu, ep_chol, A, b, 1, 1)

# test wall hitting ----
d <- 2
mu <- rep(0, d)
Sigma <- .5 * diag(d) + .5 * rep(1, d) %*% t(rep(1, d))
lb <- rep(0, d)
ub <- rep(1, d)
A <- rbind(diag(d), -diag(d))
b <- c(-lb, ub)
moments <- epmgp::moments(lb, ub, mu, Sigma)
test_wall_hit(rep(.1, d), moments$mu, chol(moments$Sigma),
              A, b, J = 2, N = 2)

d <- 2
A <- diag(d)
mu <- rep(0, d)
Sigma <- diag(d)
lb <- rep(0, d)
ub <- rep(Inf, d)

epmgp::pmvn2(mu, Sigma, lb, ub, A, log = FALSE)


d <- 2
cube <- volesti::gen_cube(d, "H")
A <- cube$A
b <- cube$b

mu <- rep(0, d)
Sigma <- diag(d)
lb <- b[1:(length(b)/2)]
ub <- b[(length(b)/2 + 1):length(b)]

epmgp::pmvn2(mu, Sigma, rep(Inf, ), b, A)
