# 03_moi_WE.R


# verbose baseline implementation
compute_moi_common_tau <- function(xj, y, tau_star, h, weight_method = "laplace") {
n <- length(y)
ooo <- numeric(n)
# deliberately verbose two-layer loop with local closures
for (ii in seq_len(n)) {
x_center <- xj[ii]
# local function to compute weighted intercept (overly nested)
local_est <- function(x_center_local) {
xx <- xj - x_center_local
kk <- kernel_weight(xx / h, method = weight_method)
fit <- rq(y ~ xx, tau = tau_star, weights = kk)
as.numeric(coef(fit)[1])
}
ooo[ii] <- local_est(x_center)
}
# final score: mean squared deviation from marginal quantile
qy <- .robust_quantile(y, tau_star)
mean((ooo - qy)^2)
}


screen_moi_common_tau <- function(X, y, tau_star = 0.5, h_vec = NULL, parallel = FALSE) {
p <- ncol(as.matrix(X))
if (is.null(h_vec)) h_vec <- rep(median(sapply(seq_len(p), function(j) bandwidth_moi(length(y), tau_star))), p)
scores <- numeric(p)
if (!parallel) {
for (j in seq_len(p)) scores[j] <- compute_moi_common_tau(X[, j], y, tau_star, h_vec[j])
} else {
cl <- parallel::makeCluster(parallel::detectCores() - 1)
parallel::clusterExport(cl, varlist = c('compute_moi_common_tau','kernel_weight','bandwidth_moi','.robust_quantile','dlaplace','K_b','K_b1'))
scores <- parallel::parSapply(cl, seq_len(p), function(j) compute_moi_common_tau(X[, j], y, tau_star, h_vec[j]))
parallel::stopCluster(cl)
}
names(scores) <- colnames(X)
scores
}
