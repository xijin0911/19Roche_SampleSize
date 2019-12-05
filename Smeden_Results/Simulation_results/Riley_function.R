#this file is used to illustrate functions correspond to three criteria by Riley Richard.

Size1 <- function(p,csR2_adj,S_VH){
  n1 <- p/((S_VH-1)*log(1-csR2_adj/S_VH))
  return(n1)
}

Size2 <- function(p,n,csR2_adj,E,d){
  max_csR2 <- 1-exp(2*(E*log(E/n)+(n-E)*log(1-E/n))/n)
  S_VH <- csR2_adj/(csR2_adj+d*max_csR2)
  n2<- p/((S_VH-1)*log(1-csR2_adj/S_VH))
}

Size3 <- function(phi,error_mar){
  n3<- (1.96/error_mar)^2*phi*(1-phi)
}

Size <- function(n1,n2,n3){
  n <- max(n1,n2,n3)
}
