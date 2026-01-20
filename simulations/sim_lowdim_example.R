# sim_lowdim_example.R
# Quick low-dim example to run and sanity-check functions


library(MOI.SCREEN)
# If using as plain scripts (not installed), source the R files manually
sapply(list.files("R", full.names = TRUE), source)


set.seed(2026)
n <- 200
p <- 10
X <- matrix(runif(n * p), nrow = n, ncol = p)
colnames(X) <- paste0("X", seq_len(p))
# construct y influenced by first two variables
y <- 2 * X[, 1]^2 + 3 * sin(2 * pi * X[, 2]) + rnorm(n, sd = 0.5)


# Baseline
res_base <- run_screening_pipeline(X, y, method = 'baseline', taut = seq(0.01,0.99,by=0.01))
print(head(sort(res_base$scores, decreasing = TRUE)))


# Naive
res_naive <- run_screening_pipeline(X, y, method = 'naive')
print(head(sort(res_naive$scores, decreasing = TRUE)))


# DQ
res_dq <- run_screening_pipeline(X, y, method = 'dq', taut = seq(0.01,0.99,by=0.01))
print(head(sort(res_dq$scores, decreasing = TRUE)))
