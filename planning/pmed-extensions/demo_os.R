## Verify one-step EIF estimator vs the bootstrap plug-in.
## Expect: no-interaction W~0 (CI covers 0); A*M interaction W!=0 (CI ~ [0.24,0.43]).
source("ward_residual_os.R")
set.seed(7); expit <- function(x) 1/(1+exp(-x))
gen <- function(n, tint, binaryY) {
  C <- rnorm(n); A <- rbinom(n, 1, expit(-0.2 + 0.8*C))
  M <- 0.6*A + 0.4*C + rnorm(n)
  lin <- 0.5*A + 0.7*M + tint*A*M + 0.3*C
  Y <- if (binaryY) rbinom(n, 1, expit(lin)) else lin + rnorm(n)
  data.frame(A, M, Y, C)
}
cat("=== Continuous Y, NO interaction (expect W~0) ===\n")
print(ward_os(gen(3000, 0.0, FALSE)))
cat("\n=== Binary Y, WITH A*M interaction (expect W!=0, CI ~ [0.24,0.43]) ===\n")
print(ward_os(gen(4000, 0.9, TRUE)))
