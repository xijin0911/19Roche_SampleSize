#this file is to extract the 4032 simulation results from van Smeden, and save these 3
#predictive error metrics seperatly.



#####################################
#preparation for data extraction#####
#####################################

library(foreign)
# create and view an object with file names and full paths
f <- file.path("/R/ssizepred/van_Smeden_Rewrite/OUT", 
               paste("main",sep="",formatC(seq(00001:04032), width = 5, flag = '0')),
               paste(
                 paste("main",sep="",formatC(seq(00001:04032), width = 5, flag = '0')),"_",
                 paste(formatC(rep(00001,length=length(seq(00001:04032))), width = 5, flag = '0')),sep="",
                 rep(".rds",length=length(seq(00001:04032)))))
#read outputs of simulation results into a list
out_list <- lapply(f, readRDS)


#names for all of the modelling strategies
strategies_categories <- unlist(lapply(out_list[1],names))  #3 categories for the modelling stageties
strategies <- c(unlist(lapply(out_list[[1]][[1]],names)),
                unlist(lapply(out_list[[1]][[2]],names)),
                unlist(lapply(out_list[[1]][[3]],names)))
names(strategies) <- c(rep(strategies_categories[1],6),
                       rep(strategies_categories[2],2),
                       rep(strategies_categories[3],2))



#names in the simulation outcomes after unlist function, which helps us to extract data
brier_names_logistf <- paste(strategies_names[1],".","validate",".",sep="",
                             unlist(lapply(out_list[[1]][[1]],names)),".","brier")
brier_names_hs <- paste(strategies_names[2],".","validate",".",sep="",
                        unlist(lapply(out_list[[1]][[2]],names)),".","brier")
brier_names_glmnet <- paste(strategies_names[3],".","validate",".",sep="",
                            unlist(lapply(out_list[[1]][[3]],names)),".","brier")
brier_names <- c(brier_names_logistf,
                 brier_names_hs,
                 brier_names_glmnet)


MAPE_names_logistf <- paste(strategies_names[1],".","validate",".",sep="",
                            unlist(lapply(out_list[[1]][[1]],names)),".","MAPE")
MAPE_names_hs <- paste(strategies_names[2],".","validate",".",sep="",
                       unlist(lapply(out_list[[1]][[2]],names)),".","MAPE")
MAPE_names_glmnet <- paste(strategies_names[3],".","validate",".",sep="",
                           unlist(lapply(out_list[[1]][[3]],names)),".","MAPE")
MAPE_names <- c(MAPE_names_logistf,
                MAPE_names_hs,
                MAPE_names_glmnet)


MSPE_names_logistf <- paste(strategies_names[1],".","validate",".",sep="",
                            unlist(lapply(out_list[[1]][[1]],names)),".","MSPE")
MSPE_names_hs <- paste(strategies_names[2],".","validate",".",sep="",
                       unlist(lapply(out_list[[1]][[2]],names)),".","MSPE")
MSPE_names_glmnet <- paste(strategies_names[3],".","validate",".",sep="",
                           unlist(lapply(out_list[[1]][[3]],names)),".","MSPE")
MSPE_names <- c(MSPE_names_logistf,
                MSPE_names_hs,
                MSPE_names_glmnet)


#####################################
################data extraction######
#####################################


setwd("/R/ssizepred/weekly_work/week6/simulation/Meta_model")
SimulationSet <- readRDS("SimulationSet.rds")

#matrixes to save the outcome
brier_simulation <- matrix(0,nrow=nrow(SimulationSet),ncol=10)
colnames(brier_simulation) <- brier_names
mape_simulation <- matrix(0,nrow=nrow(SimulationSet),ncol=10)
colnames(mape_simulation) <- MAPE_names
mspe_simulation <- matrix(0,nrow=nrow(SimulationSet),ncol=10)
colnames(mspe_simulation) <- MSPE_names

for (i in 1:4302){
  allvalues <- unlist(out_list[i])
  brier_simulation[i,] <- as.matrix(allvalues[brier_names])
  mape_simulation[i,] <- as.matrix(allvalues[MAPE_names])
  mspe_simulation[i,] <- as.matrix(allvalues[MSPE_names])
}

#name the 3 matrixes
matrix_name <- c("Firth","ML","Firthp","MLp","Firthaic","MLaic","HS1","HS2","Ridge","Lasso")
colnames(brier_simulation) <-  matrix_name 
colnames(mape_simulation) <-   matrix_name
colnames(mspe_simulation) <-   matrix_name

setwd("/R/ssizepred/weekly_work/week6/simulation/Simulation_result")
saveRDS(brier_simulation, file="brier_simulation.rds")
saveRDS(mape_simulation, file="mape_simulation.rds")
saveRDS(mspe_simulation, file="mspe_simulation.rds")
























