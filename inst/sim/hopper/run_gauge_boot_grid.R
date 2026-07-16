## Gauge bootstrap-CI coverage grid (Hopper). Analytic (Wald) vs percentile bootstrap
## for the ratios W=R/OE and P_med=IIE/OE. nsim=2000 per cell (8 chunks x 250).
## Gates: per-cell MCSE, weak-ID (Wald-vs-percentile divergence), failed-run logging,
## known-answer (tau=0 => W=0), SE-vs-estimate data. SEE SIM-REQUEST-gauge-bootstrap-grid.md
suppressMessages(library(probmed))
## --- stale-package guard (A1): abort unless the bootstrap arm is live ---
.chk <- local({ set.seed(99); n<-300; C<-rnorm(n); A<-rbinom(n,1,plogis(0.3*C))
  M<-0.5*A+0.3*C+rnorm(n); Y<-0.4*A+0.6*M+rnorm(n); d<-data.frame(A,M,Y,C)
  a<-ward_residual(d,se_method="analytic"); b<-ward_residual(d,se_method="bootstrap",B=50L)
  !isTRUE(all.equal(a@W_ci,b@W_ci)) })
if (!isTRUE(.chk)) stop("STALE probmed: se_method='bootstrap' is a no-op. Check R_LIBS ordering.")

expit <- function(x) 1/(1+exp(-x))
gen <- function(n,tint,binY){C<-rnorm(n);A<-rbinom(n,1,expit(-0.2+0.8*C));M<-0.6*A+0.4*C+rnorm(n)
  lin<-0.5*A+0.7*M+tint*A*M+0.3*C; Y<-if(binY) rbinom(n,1,expit(lin)) else lin+rnorm(n); data.frame(A,M,Y,C)}
truth <- function(tint,binY,N=2e6){C<-rnorm(N)
  th<-function(a,ap){M<-0.6*ap+0.4*C+rnorm(N);lin<-0.5*a+0.7*M+tint*a*M+0.3*C;mean(if(binY) expit(lin) else lin)}
  t11<-th(1,1);t10<-th(1,0);t01<-th(0,1);t00<-th(0,0);OE<-t11-t00
  c(W=(OE-(t10-t00)-(t01-t00))/OE, P=(t01-t00)/OE)}

cells <- expand.grid(n=c(800,3000), tint=c(0,0.8), binY=c(FALSE,TRUE))   # 8 cells
NCHUNK <- 8L; REPS_PER <- 250L; B <- 999L
aid <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID","1"))                 # 1..64
ci <- ((aid-1L) %/% NCHUNK) + 1L                                          # cell 1..8
ch <- ((aid-1L) %%  NCHUNK) + 1L                                          # chunk 1..8
n<-cells$n[ci]; tint<-cells$tint[ci]; binY<-cells$binY[ci]
tr <- truth(tint,binY)
base <- 100000L*ci + 1000L*ch
rows <- vector("list", REPS_PER); nfail <- 0L; nsing <- 0L
for (j in seq_len(REPS_PER)) {
  seed <- base + j; set.seed(seed); d <- gen(n,tint,binY)
  fa <- tryCatch(ward_residual(d, seed=seed), error=function(e) NULL)
  fb <- tryCatch(ward_residual(d, seed=seed, se_method="bootstrap", B=B), error=function(e) NULL)
  if (is.null(fa) || is.null(fb)) { nfail<-nfail+1L; next }
  if (abs(fa@OE) < 1e-3) nsing <- nsing + 1L                              # near-singular OE (A1)
  divW <- max(abs(fa@W_ci - fb@W_ci))                                     # weak-ID divergence (A2)
  rows[[j]] <- data.frame(cell=ci, chunk=ch, seed=seed, n=n, tint=tint, binY=binY,
    trW=tr["W"], trP=tr["P"], W=fa@W, P=fa@p_med, seW_an=fa@W_se, seW_bt=fb@W_se,
    covW_an = tr["W"]>=fa@W_ci[1] && tr["W"]<=fa@W_ci[2],
    covW_pct= tr["W"]>=fb@W_ci[1] && tr["W"]<=fb@W_ci[2],
    covP_an = tr["P"]>=fa@p_med_ci[1] && tr["P"]<=fa@p_med_ci[2],
    covP_pct= tr["P"]>=fb@p_med_ci[1] && tr["P"]<=fb@p_med_ci[2],
    divW_WaldVsPct = divW, row.names=NULL)
}
out <- do.call(rbind, rows)
attr(out,"nfail")<-nfail; attr(out,"nsing")<-nsing
dir.create("~/gauge_boot/parts", recursive=TRUE, showWarnings=FALSE)
f <- sprintf("~/gauge_boot/parts/part_cell%02d_chunk%02d.rds", ci, ch)
saveRDS(out, f)
cat(sprintf("cell %d chunk %d done: %d rows, nfail=%d nsing=%d -> %s\n",
            ci, ch, nrow(out), nfail, nsing, f))
