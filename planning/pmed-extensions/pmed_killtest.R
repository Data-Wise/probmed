## Prop 1 kill-test: Incremental Mediated Elasticity / proportion-mediated curve
## Q1: does P_med^inc reduce to classical P_med in linear-Gaussian (exact, zero remainder)?
## Q2: does the curve flex (and remainder go nonzero) under nonlinearity+interaction?
set.seed(1)
expit <- function(x) 1/(1+exp(-x))

## --- parameters ---
g0<- -0.2; gc<- 0.8                 # logistic exposure propensity
b0<- 0.0; ba<- 0.6; bc<- 0.4; sM<-1 # mediator model M = b0+ba*A+bc*C+N(0,sM^2)
t0<- 0.0; ta<- 0.5; tm<- 0.7; tc<- 0.3  # outcome model
tint<- 0.8                          # A*M interaction (used in nonlinear case)

## classical natural-effect P_med (linear, no interaction)
NIE<-ba*tm; NDE<-ta; Pmed_classical<-NIE/(NIE+NDE)

## --- common random numbers ---
N<-2e6
C<-rnorm(N); g<-expit(g0+gc*C)
U_dir<-runif(N); U_med<-runif(N); epsM<-rnorm(N,0,sM); U_Y<-runif(N)
qd<-function(d) d*g/(d*g+1-g)
Adraw<-function(d,U) as.numeric(U < qd(d))
Mgen<-function(A) b0+ba*A+bc*C+epsM

## linear continuous outcome (noise cancels under common RNG; omit for precision)
Ylin<-function(Adir,Amed){ M<-Mgen(Amed); t0+ta*Adir+tm*M+tc*C }
EYlin<-function(ddir,dmed) mean(Ylin(Adraw(ddir,U_dir),Adraw(dmed,U_med)))

## nonlinear: binary outcome with A*M interaction
Ylog<-function(Adir,Amed){ M<-Mgen(Amed); p<-expit(t0+ta*Adir+tm*M+tint*Adir*M+tc*C); as.numeric(U_Y<p) }
EYlog<-function(ddir,dmed) mean(Ylog(Adraw(ddir,U_dir),Adraw(dmed,U_med)))

## interventional incremental decomposition vs baseline delta=1
decomp<-function(EY,d){
  base<-EY(1,1); OE<-EY(d,d)-base; IDE<-EY(d,1)-base; IIE<-EY(1,d)-base
  c(delta=d,OE=OE,IDE=IDE,IIE=IIE,remainder=OE-IDE-IIE,Pmed_inc=IIE/OE)
}
deltas<-c(0.5,1.5,2,3)
cat("==== LINEAR-GAUSSIAN (expect remainder~0, Pmed_inc constant = classical) ====\n")
print(round(t(sapply(deltas,function(d) decomp(EYlin,d))),5))
cat("\nClassical natural-effect P_med =",round(Pmed_classical,5),"\n\n")
cat("==== NONLINEAR (logistic Y, A*M interaction): curve should flex, remainder != 0 ====\n")
print(round(t(sapply(deltas,function(d) decomp(EYlog,d))),5))
