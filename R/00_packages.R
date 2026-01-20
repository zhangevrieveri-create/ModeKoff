.pkg_load <- function() {
# 为了在不同机器上也能加载，采用逐个检查并提示的方式
pkgs <- c("quantreg", "energy", "MASS", "parallel", "stats", 'tidyverse', 'wooldridge', 'AER', 'ggplot2', 'rqPen', 'splines', 'VGAM', 'BwQuant')
for (p in pkgs) {
if (!requireNamespace(p, quietly = TRUE)) {
stop(sprintf("Required package '%s' is not installed. Install it with install.packages('%s')", p, p))
}
}
# attach a few common functions to namespace
invisible(lapply(pkgs, require, character.only = TRUE))
}


.pkg_load()
