## Demo: ward_residual() should give W ~ 0 (no interaction) and W != 0 (with A*M).
## Small K/B here for speed; raise for real use.
source("ward_residual.R")
set.seed(7); expit <- function(x) 1/(1+exp(-x))
gen <- function(n, tint, binaryY) {
  C <- rnorm(n); A <- rbinom(n, 1, expit(-0.2 + 0.8*C))
  M <- 0.6*A + 0.4*C + rnorm(n)
  lin <- 0.5*A + 0.7*M + tint*A*M + 0.3*C
  Y <- if (binaryY) rbinom(n, 1, expit(lin)) else lin + rnorm(n)
  data.frame(A, M, Y, C)
}
cat("=== Continuous Y, NO interaction (expect W~0, CI covers 0) ===\n")
print(ward_residual(gen(1500, 0.0, FALSE), K=40, B=50))
cat("\n=== Binary Y, WITH A*M interaction (expect W!=0) ===\n")
print(ward_residual(gen(2500, 0.9, TRUE), K=40, B=50))
