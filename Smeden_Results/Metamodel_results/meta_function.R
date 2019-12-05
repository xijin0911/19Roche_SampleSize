#this file is for the function used to calculate 3 prediction error metrics by meta-model


#####################################################
#########function for meta-model result based on van Smeden

MPSE <- function(n,ef,P){
  ML <- -0.59 - 1.06*log(n) + 0.36*log(ef) + 0.94*log(P)
  Firth <- -0.86  - 1.03*log(n) + 0.33*log(ef) + 0.93*log(P) 
  HS <- -0.75 - 0.97*log(n) + 0.44*log(ef) + 0.74*log(P) 
  Lasso <- -0.86 - 0.93*log(n) + 0.46*log(ef) + 0.68*log(P) 
  Ridge  <- -0.93 - 1.88*log(n) + 0.50*log(ef) + 0.49*log(P) 
  MLp <- -0.57 - 1.03*log(n) + 0.4*log(ef) + 0.96*log(P) 
  MLaic <- -0.59 - 1.05*log(n) + 0.38*log(ef) + 0.95*log(P) 
  Firthp <- -0.66 - 1.01*log(n) + 0.39*log(ef) + 0.95*log(P) 
  Firthaic <- -0.74 - 1.03*log(n) + 0.37*log(ef) + 0.95*log(P) 
  out <- exp(c(ML,Firth,HS,Lasso,Ridge, MLp,MLaic,Firthp,Firthaic))
  names(out) <-c("ML","Firth","HS","Lasso","Ridge", "MLp","MLaic","Firthp","Firthaic")
  return(out)
}

MAPE <- function(n,ef,P){
  ML <- -0.48 - 0.53*log(n) + 0.31*log(ef) + 0.48*log(P)
  Firth <- -0.61  - 0.51*log(n) + 0.29*log(ef) + 0.47*log(P) 
  HS <- -0.56 - 0.49*log(n) + 0.33*log(ef) + 0.39*log(P) 
  Lasso <- -0.59 - 0.48*log(n) + 0.34*log(ef) + 0.35*log(P) 
  Ridge  <- -0.61 - 0.45*log(n) + 0.36*log(ef) + 0.26*log(P) 
  MLp <- -0.45 - 0.52*log(n) + 0.31*log(ef) + 0.49*log(P) 
  MLaic <- -0.48 - 0.53*log(n) + 0.31*log(ef) + 0.49*log(P) 
  Firthp <- -0.49 - 0.52*log(n) + 0.30*log(ef) + 0.50*log(P) 
  Firthaic <- -0.55 - 0.52*log(n) + 0.30*log(ef) + 0.49*log(P) 
  out <- exp(c(ML,Firth,HS,Lasso,Ridge, MLp,MLaic,Firthp,Firthaic))
  names(out) <-c("ML","Firth","HS","Lasso","Ridge", "MLp","MLaic","Firthp","Firthaic")
  return(out)
}

Brier <- function(n,ef,P){
  ML <- -0.91 - 0.04*log(n) + 0.62*log(ef) + 0.04*log(P) 
  Firth <- -0.95  - 0.03*log(n) + 0.62*log(ef) + 0.03*log(P) 
  HS <- -0.93 - 0.03*log(n) + 0.62*log(ef) + 0.02*log(P) 
  Lasso <- -0.96 - 0.03*log(n) + 0.62*log(ef) + 0.02*log(P) 
  Ridge  <- -0.98 - 0.02*log(n) + 0.62*log(ef) + 0.01*log(P) 
  MLp <- -0.89 - 0.04*log(n) + 0.62*log(ef) + 0.04*log(P) 
  MLaic <- -0.91 - 0.04*log(n) + 0.62*log(ef) + 0.04*log(P) 
  Firthp <- -0.90 - 0.04*log(n) + 0.62*log(ef) + 0.04*log(P) 
  Firthaic <- -0.93 - 0.04*log(n) + 0.62*log(ef) + 0.04*log(P) 
  out <- exp(c(ML,Firth,HS,Lasso,Ridge, MLp,MLaic,Firthp,Firthaic))
  names(out) <-c("ML","Firth","HS","Lasso","Ridge", "MLp","MLaic","Firthp","Firthaic")
  return(out)
}

