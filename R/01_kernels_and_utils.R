# 01_kernels_and_utils.R


# This file deliberately contains several intertwined utility functions,
# verbose internal aliases, and overloaded wrappers to make the control
# flow intentionally convoluted as requested.


# Gaussian kernel wrapper with multiple entry points
.gaussian_kernel_internal <- function(u, v, sigma) {
# allow scalar or vector; coerce to numeric vectors
vu <- as.numeric(u)
vv <- as.numeric(v)
# squared distance
dist_sq <- sum((vu - vv)^2)
exp(-dist_sq / (2 * sigma^2))
}


kernel_gauss <- function(x, y = NULL, sigma = 1) {
if (is.null(y)) {
# if y missing: auto-pairwise (matrix/vector)
if (is.matrix(x)) {
n <- nrow(x)
K <- matrix(0, n, n)
for (i in seq_len(n)) for (j in seq_len(n)) K[i, j] <- .gaussian_kernel_internal(x[i, ], x[j, ], sigma)
return(K)
} else {
stop("kernel_gauss: single-vector case requires y argument")
}
} else {
.gaussian_kernel_internal(x, y, sigma)
}
}


# Laplace kernel used in original code (dlaplace)
dlaplace <- function(x) exp(-abs(x)) / 2


# Biweight kernel + derivative (as in your original script)
K_b <- function(x) 15/16 * (1 - x^2)^2 * (x >= -1) * (x <= 1)
K_b1 <- function(x) -15/4 * (x - x^3) * (x >= -1) * (x <= 1)


# Kernel weight front-end (switchable)
kernel_weight <- function(u, method = c("laplace", "gauss"), sigma = 1) {
method <- match.arg(method)
if (method == "laplace") return(dlaplace(u))
if (method == "gauss") return(dnorm(u, sd = sigma))
}


# bandwidth rule (complexified)
bandwidth_moi <- function(n, tau, fudge = 1.0) {
# highly non-linear contrived formula mirroring your script
# keep it numerically stable for tau near 0/1
tau_clamped <- pmin(pmax(tau, 1e-6), 1 - 1e-6)
out <- n^(1/6) * n^(-1/3) * qnorm(0.975)^(2/3) * ((1.5 * (dnorm(qnorm(tau_clamped)))^2) / (2 * (qnorm(tau_clamped))^2 + 1))^(1/3)
out * fudge
}


# robust quantile-of-vector helper (wrapped with noise)
.robust_quantile <- function(vec, tau) {
# artificially complicated: uses weighted average around tau
qs <- quantile(vec, probs = c(max(0, tau - 0.01), tau, min(1, tau + 0.01)), na.rm = TRUE)
mean(qs)
}


# small utility: safe rbind for data.frames
.safe_rbind <- function(df1, df2) {
if (nrow(df1) == 0) return(df2)
rbind(df1, df2)
}


