# This file is part of combination of R-files that contain code to perform a simulation study
# The code was used for a simulation study for which the results are presented in a currently unpublished manuscript

########################################
#### Functions to estimate models   ####
#### Last modified: Feb 8 2018      ####
#### Author: M van Smeden           ####
########################################

# BEGIN CODE ---------------------------------------


########################################
#### functions to generate data     ####
########################################

gen.MVNX <- function(n,mu,sigma){           #simulation setting is based on multivariate distributed predictors with no correlation
  mvrnorm(n=n,mu=mu,Sigma=sigma)
}

#generate dependent variables
gen.binY <- function(SIMx,dgm.par){
  no.X <- ncol(SIMx)
  design.mat <- cbind(1,SIMx)
  p <- exp(design.mat%*%dgm.par)/(1+exp(design.mat%*%dgm.par))
  y <- rbinom(length(p),1,p)
  SIMxy <-data.frame(x=SIMx,y=y)
  colnames(SIMxy)[1:no.X]<-paste("X",1:no.X,sep="")  
  SIMxy
}


##############################################################
#### calculate sample size satisfying three-criteria      ####
##############################################################

size1 <- function(S,P,csR2){              #relative drop
  s1 <- P/((S-1)*log(1-csR2/S))         
  round(s1,digits = 0)
}

size2 <- function(E,n,csR2,d){          #absolute difference
  max_csR2 <- 1-exp(2*(E*log(E/n)+
                         (n-E)*log(1-E/n))/n)
  S_VH <- csR2/(csR2+d*max_csR2)      #S_vh is what we calculatebased on relative difference, different the desired shrinkage factor
  s2 <- P/((S_VH-1)*log(1-csR2/S_VH))
  round(s2,digits = 0)
}

size3 <- function(error_margin,E,n){        #precise estimation
  n3 <- ((1.96/error_mar)^2)*(E/n)*(1-(E/n))
  round(n3,digits=0)
}

##############################################################
#### Model development based on heuristic shrinkage       ####
##############################################################
tryCatch.W.E <- function(expr){
  W <- NULL
  w.handler <- function(w){
    W <<- w
    invokeRestart("muffleWarning")}
  list(value = withCallingHandlers(tryCatch(expr, error = function(e) e), warning = w.handler),warning = W)
}


heuristicshrink.lrm <- function(SIMxy){
  formu <- as.formula(paste("y~", paste(colnames(SIMxy)[-which(colnames(SIMxy)=="y")],collapse="+"),sep=""))
  SIMx <- SIMxy[,-ncol(SIMxy)]
  y <- SIMxy[,ncol(SIMxy)]
  
  TEMP.heur <- list()
  unpen.fit <- glm(formu,data=SIMxy,family=binomial(link = "logit"))
  TEMP.heur$int.re.est <- tryCatch.W.E(logistf(formu,data=SIMxy,firth = FALSE,alpha = 0.05, dataout = T))  # I set the value of alphs (significance level) as 0.1
  chisq <-unpen.fit$null.deviance-unpen.fit$deviance
  s <- (chisq-TEMP.heur$int.re.est$value$df)/chisq
  
  if(chisq!=0){
    A <- TEMP.heur$int.re.est$value$coefficients["(Intercept)"]     #original intercept 
    B <- TEMP.heur$int.re.est$value$coefficients[-which(names(TEMP.heur$int.re.est$value$coefficients)=="(Intercept)")]   # original coefficient for predictors
    
    TEMP.heur$int.re.est$value$coefficients[-which(names(TEMP.heur$int.re.est$value$coefficients)=="(Intercept)")] <- B*s
    
    TEMP.heur$int.approx <- TEMP.heur$int.re.est    
    TEMP.heur$int.approx$value$coefficients["(Intercept)"] <- (1-s)*mean(SIMxy$y)+s*A
    
    off <- data.matrix(SIMx)%*%data.matrix(B*s)
    
    TEMP.heur$int.re.est$value$coefficients["(Intercept)"] <- coefficients(glm(SIMxy$y~offset(off), family = binomial(link = "logit")))["(Intercept)"]
  }else{TEMP.heur$int.approx <- TEMP.heur$int.re.est}
  
  list(OUT=TEMP.heur,hs.lambda=s)         # lambda is the heuristic shirnkage factor 
}

#function for prediction of maximum likelihood model
predict.lrm <- function(lp){  
  as.vector(exp(lp)/(1+exp(lp)))
} 

##############################################################
#### predictive error metrics calculation                 ####
##############################################################
brier_true <- function(p,SIMxy){
  sum((SIMxy[,"y"]-p)^2)/length(p)
}

mspe_true <- function(p,SIMxy,dgm.par){
  design.mat <- as.matrix(cbind(1,SIMxy[,-which(colnames(SIMxy)=="y")]))
  p_true <- exp(design.mat%*%dgm.par)/(1+exp(design.mat%*%dgm.par))
  mean((p_true-p)^2)
}

mape_true <- function(p,SIMxy,dgm.par){
  design.mat <- as.matrix(cbind(1,SIMxy[,-which(colnames(SIMxy)=="y")]))
  p_true <- exp(design.mat%*%dgm.par)/(1+exp(design.mat%*%dgm.par))
  mean(abs(p_true-p))
}



###############################################################
#### predictive error metrics calculation by metamodel for ML #      
###############################################################
mspe_ML <- function(MSPE,Ef,P){
  s1 <- exp((-0.59 - log(MSPE) + 0.36*log(Ef) + 0.94*log(P))/1.06)
  round(s1,digits = 0 )
}

mape_ML <- function(MAPE,Ef,P){
  s2 <- exp((-0.48 - log(MAPE) + 0.31*log(Ef) + 0.48*log(P))/0.53)
  round(s2,digits = 0 )
}

brier_ML <- function(Brier,Ef,P){
  s3 <- exp((-0.91 - log(Brier) + 0.62*log(Ef) + 0.04*log(P))/0.04)
  round(s3,digits = 0 )
}  

###############################################################
#### predictive error metrics calculation by metamodel for HS #      
###############################################################

mspe_HS <- function(MSPE,Ef,P){
  s1 <- exp((-0.75 - log(MSPE) + 0.44*log(Ef) + 0.74*log(P))/(0.97))
  round(s1,digits = 0 )
}

mape_HS <- function(MAPE,Ef,P){
  s2 <- exp((-0.56 - log(MAPE) + 0.33*log(Ef) + 0.39*log(P))/(0.49))
  round(s2,digits = 0 )
}

brier_HS <- function(Brier,Ef,P){
  s3 <- exp((-0.93 - log(Brier) + 0.62*log(Ef) + 0.02*log(P))/(0.03))
  round(s3,digits = 0 )
}


# END CODE ---------------------------------------------