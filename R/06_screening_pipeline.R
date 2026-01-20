# 06_screening_pipeline.R


# Large orchestrating function that wires everything together.
run_screening_pipeline <- function(X, y,
method = c("baseline", "naive", "dq"),
taut = seq(0.01, 0.99, by = 0.01),
parallel = FALSE,
top_k = NULL) {
method <- match.arg(method)
p <- ncol(as.matrix(X))
# compute default h_vec
h_vec <- rep(bandwidth_moi(length(y), 0.5), p)


if (method == "baseline") {
# estimate a single tau_star by pooling simplified per-dimension estimates
per_dim_tau <- sapply(seq_len(p), function(j) {
# estimate per-dim tau by simple local-min-variance heuristic
estimate_tau_ij(X[, j], y, x0 = median(X[, j]), taut = taut, h_local = h_vec[j])
})
tau_star <- mean(per_dim_tau)
scores <- screen_moi_common_tau(X, y, tau_star = tau_star, h_vec = h_vec, parallel = parallel)
return(list(method = "baseline", tau_star = tau_star, scores = scores, per_dim_tau = per_dim_tau))
}


if (method == "naive") {
scores <- screen_moi_naive(X, y, h_vec = h_vec, tau = 0.5, parallel = parallel)
return(list(method = "naive", tau_star = 0.5, scores = scores))
}


if (method == "dq") {
scores <- screen_moi_dq(X, y, taut = taut, h_vec = h_vec, parallel = parallel)
return(list(method = "dq", scores = scores))
}
}
