#this file used to extract the matrix for 3 predictive error metrix from the result of meta model.


library(readr)
setwd("/R/ssizepred/van_Smeden")
SimulationSet <- read_delim("prelim/maincond.csv",";", escape_double = FALSE, trim_ws = TRUE)

#this matrix saved the characteristics when we do simulation
#P: number of predictors
#n: sample size for model development
#ef: events fraction
SimulationSet <- cbind(
                       P=SimulationSet$P,
                       n=SimulationSet$EPV*SimulationSet$P,
                       ef=SimulationSet$prevy1)

#three metrix to save the outcomes of meta-model
mspe_meta <- matrix(0,nrow=nrow(SimulationSet),ncol=9)
mape_meta <- matrix(0,nrow=nrow(SimulationSet),ncol=9)
brier_meta <- matrix(0,nrow=nrow(SimulationSet),ncol=9)

setwd("/R/ssizepred/weekly_work/week6/simulation/Meta_model")
source("meta_function.R")
for (i in 1:nrow(SimulationSet)){
  mspe_meta[i,] <- MPSE(n=SimulationSet[i,"n"],
                        ef=SimulationSet[i,"ef"],
                        P=SimulationSet[i,"P"])
}
colnames(mspe_meta) <- c("ML","Firth","HS","Lasso","Ridge", "MLp","MLaic","Firthp","Firthaic")


for (i in 1:nrow(SimulationSet)){
  mape_meta[i,] <- MPSE(n=SimulationSet[i,"n"],
                        ef=SimulationSet[i,"ef"],
                        P=SimulationSet[i,"P"])
}
colnames(mape_meta) <- c("ML","Firth","HS","Lasso","Ridge", "MLp","MLaic","Firthp","Firthaic")


for (i in 1:nrow(SimulationSet)){
  brier_meta[i,] <- MPSE(n=SimulationSet[i,"n"],
                         ef=SimulationSet[i,"ef"],
                         P=SimulationSet[i,"P"])
}
colnames(brier_meta) <- c("ML","Firth","HS","Lasso","Ridge", "MLp","MLaic","Firthp","Firthaic")

#save these 3 matixes with respect to 3 predictive error metrics

saveRDS(brier_meta, file="brier_meta.rds")
saveRDS(mape_meta, file="mape_meta.rds")
saveRDS(mspe_meta, file="mspe_meta.rds")
#save simulation setting
saveRDS(SimulationSet, file="SimulationSet.rds")

