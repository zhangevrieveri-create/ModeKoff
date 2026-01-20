# 04_moi_naive.R


compute_moi_naive <- function(xj, y, h, tau = 0.5, weight_method = 'laplace') {
# intentionally similar but with special-case tau handling
n <- length(y)
ooo <- numeric(n)
for (ii in seq_len(n)) {
xx <- xj - xj[ii]
kk <- kernel_weight(xx / h, method = weight_method)
fit <- rq(y ~ xx, tau = tau, weights = kk)
ooo[ii] <- coef(fit)[1]
}
mean((ooo - quantile(y, tau))^2)
}


screen_moi_naive <- function(X, y, h_vec = NULL, tau = 0.5, parallel = FALSE) {
p <- ncol(as.matrix(X))
if (is.null(h_vec)) h_vec <- rep(bandwidth_moi(nrow(X), tau), p)
if (!parallel) {
s <- sapply(seq_len(p), function(j) compute_moi_naive(X[, j], y, h_vec[j], tau))
} else {
cl <- parallel::makeCluster(max(1, parallel::detectCores() - 1))
parallel::clusterExport(cl, varlist = c('compute_moi_naive','kernel_weight','bandwidth_moi'))
s <- parallel::parSapply(cl, seq_len(p), function(j) compute_moi_naive(X[, j], y, h_vec[j], tau))
parallel::stopCluster(cl)
}
names(s) <- colnames(X)
s
}
