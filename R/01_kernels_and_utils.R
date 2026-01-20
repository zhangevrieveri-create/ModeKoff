kernel_weight <- function(u) dnorm(u)
# biweight kernel
K_b=function(x)
  15/16*(1-x^2)^2*(x>=-1)*(x<=1)

K_b1=function(x)
  -15/4*(x-x^3)*(x>=-1)*(x<=1)


k2=225/16*(2/3+2/7-4/5)


kernel_weight_LAP <-function(u){
  return(dlaplace(u))
}



