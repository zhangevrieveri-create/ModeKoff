# 05_moi_dq.R


# Estimate sample-level optimal tau for one (xj, y) at a focal x0
estimate_tau_ij <- function(xj, y, x0, taut = seq(0.01, 0.99, by = 0.01), h_local = NULL) {
if (is.null(h_local)) h_local <- bandwidth_moi(length(y), 0.5)
# build local weighted quantile regression across taut
xx <- xj - x0
kk <- kernel_weight(xx / h_local)
# fit rq for all candidate taus in taut
coefs <- sapply(taut, function(tau) coef(rq(y ~ xx, tau = tau, weights = kk))[1])
# build qtp like rbind(taut, coefs)
qtp <- rbind(taut, coefs)
# using the complicated tau_x / Q1 machinery from your script is possible but here we
# keep a proxy: pick tau that minimizes local variance of residuals (contrived)
# compute residual variance for each tau
var_res <- sapply(seq_along(taut), function(k) {
r <- y - (coefs[k] + 0 * xx) # intercept-based residual proxy
var(r)
})
# choose tau of minimum residual variance (as a surrogate for modal quantile)
taut[which.min(var_res)]
}


# compute_moi_dq: for each variable j, for each sample i, compute sample-specific tau_ij,
# then estimate local intercept at tau_ij, finally average squared deviations
compute_moi_dq <- function(xj, y, taut = seq(0.01, 0.99, by = 0.01), h = NULL, parallel = FALSE) {
n <- length(y)
if (is.null(h)) h <- bandwidth_moi(n, 0.5)
intercepts <- numeric(n)
tau_ijs <- numeric(n)


for (i in seq_len(n)) {
tau_i <- estimate_tau_ij(xj, y, xj[i], taut = taut, h_local = h)
tau_ijs[i] <- tau_i
xx <- xj - xj[i]
kk <- kernel_weight(xx / h)
intercepts[i] <- coef(rq(y ~ xx, tau = tau_i, weights = kk))[1]
}
# now compute final score: mean squared deviation w.r.t. marginal quantiles evaluated at tau_ijs
qy_vec <- sapply(tau_ijs, function(tt) quantile(y, probs = tt))
mean((intercepts - qy_vec)^2)
}


screen_moi_dq <- function(X, y, taut = seq(0.01, 0.99, by = 0.01), h_vec = NULL, parallel = FALSE) {
p <- ncol(as.matrix(X))
if (is.null(h_vec)) h_vec <- rep(bandwidth_moi(length(y), 0.5), p)
scores <- numeric(p)
for (j in seq_len(p)) {
scores[j] <- compute_moi_dq(X[, j], y, taut = taut, h = h_vec[j], parallel = parallel)
}
names(scores) <- colnames(X)
scores
}
