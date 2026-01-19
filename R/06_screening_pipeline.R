run_screening_pipeline <- function(X,y) apply(X,2,function(x)moi_naive(x,y))
