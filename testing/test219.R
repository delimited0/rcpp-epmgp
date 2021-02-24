load("testing/patrick0219.RData")

k = 2

u_k = unname(unlist(psi_df[k,1:D]))

# H_k = pracma::hessian(slow_psi, u_k, params = params) # numerical hessian
H_k = hess_gwish(u_k, params) # closed form computation of hessian -> precision
H_k_inv = chol2inv(chol(H_k)) # inverse precision -> covariance

# lambda_k = pracma::grad(slow_psi, u_k, params = params) # numerical gradient
lambda_k = grad_gwish(u_k, params) # closed form computation of gradient

b_k = H_k %*% u_k - lambda_k
m_k = H_k_inv %*% b_k

lb = bounds[k, seq(1, 2 * D, 2)] %>% unname %>% unlist
ub = bounds[k, seq(2, 2 * D, 2)] %>% unname %>% unlist

# this next line fails
epmgp::pmvn(lb, ub, m_k, H_k_inv, log = TRUE)

# something is wrong with trunc_norm_moments again
tau_cavity <- c(1.8823,
                1.0000,
                2.0472,
                1.0000,
                1.0000,
                1.9634,
                1.0000,
                1.0000,
                1.0000,
                2.0169,
                1.0000,
                1.0000,
                1.0000,
                1.0000,
                1.9679)
nu_cavity <- c(4.5960e+01,
               1.9766e-01,
               4.5572e+01,
               9.3553e-15,
               -3.9175e-01,
               4.3553e+01,
               -4.8540e-02,
               3.1137e-02,
               9.2677e-03,
               4.4701e+01,
               -2.2014e-01,
               -6.1844e-03,
               1.0442e+00,
               3.0310e-01,
               4.4192e+01)

trunc_norm_moments(lb[4], ub[4], nu_cavity[4] / tau_cavity[4], 1 / tau_cavity[4])

a <- (lb[4] - nu_cavity[4] / tau_cavity[4]) / sqrt(2 * 1)
b <- (ub[4] - nu_cavity[4] / tau_cavity[4]) / sqrt(2 * 1)
erfcx(a) - exp(-(b^2 - a^2)) * erfcx(b)

G_k[k] = epmgp::pmvn(lb, ub, m_k, H_k_inv, log = TRUE)

log_terms[k] = D / 2 * log(2 * pi) - 0.5 * log_det(H_k) - 
    psi_df$psi_u[k] + sum(lambda_k[,k] * u_k) - 0.5 * t(u_k) %*% H_k %*% u_k + 
    0.5 * t(m_k) %*% H_k %*% m_k + G_k[k]



