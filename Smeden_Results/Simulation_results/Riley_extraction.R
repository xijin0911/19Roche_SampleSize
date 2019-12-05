#this file is to extract the corresponding values of pseudo R-squared from simulaion 


library(foreign)
# create and view an object with file names and full paths
f <- file.path("~/R/ssizepred/van_Smeden_Rewrite/OUT", 
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
cox_names_logistf <- paste(strategies_categories[1],".","validate",".",sep="",
                           unlist(lapply(out_list[[1]][[1]],names)),".","Rsqrs.cox")
cox_names_hs <- paste(strategies_categories[2],".","validate",".",sep="",
                      unlist(lapply(out_list[[1]][[2]],names)),".","Rsqrs.cox")
cox_names_glmnet <- paste(strategies_categories[3],".","validate",".",sep="",
                          unlist(lapply(out_list[[1]][[3]],names)),".","Rsqrs.cox")
cox_names <- c(cox_names_logistf,
               cox_names_hs,
               cox_names_glmnet)

nagelkerke_names_logistf <- paste(strategies_categories[1],".","validate",".",sep="",
                                  unlist(lapply(out_list[[1]][[1]],names)),".","Rsqrs.nagelkerke")
nagelkerke_names_hs <- paste(strategies_categories[2],".","validate",".",sep="",
                             unlist(lapply(out_list[[1]][[2]],names)),".","Rsqrs.nagelkerke")
nagelkerke_names_glmnet <- paste(strategies_categories[3],".","validate",".",sep="",
                                 unlist(lapply(out_list[[1]][[3]],names)),".","Rsqrs.nagelkerke")
nagelkerke_names <- c(nagelkerke_names_logistf,
                      nagelkerke_names_hs,
                      nagelkerke_names_glmnet)

#S_vh_names <- paste(strategies_categories[3],".","validate",".",sep="",
#                    unlist(lapply(out_list[[1]][[3]],names)),".","shrinkage")


coxR2_simulation <- matrix(0,nrow=4032,ncol=10)
nagelkerkeR2_simulation <- matrix(0,nrow=4032,ncol=10)

#S_vh <- matrix(0,nrow=4032,ncol=2)

for (i in 1:4302){
  allvalues <- unlist(out_list[i])
  coxR2_simulation[i,] <- as.matrix(allvalues[cox_names])
  nagelkerkeR2_simulation[i,] <- as.matrix(allvalues[nagelkerke_names])
  S_vh[i,] <- as.matrix(allvalues[S_vh_names])
}

#give the name of each modelling strategies 
matrix_name <- c("Firth","ML","Firthp","MLp","Firthaic","MLaic","HS1","HS2","Ridge","Lasso")
colnames(coxR2_simulation) <-  matrix_name 
colnames(nagelkerkeR2_simulation) <-   matrix_name
#colnames(S_vh) <- c("Ridge","Lasso")

saveRDS(coxR2_simulation,file="coxR2_simulation.rds")
saveRDS(nagelkerkeR2_simulation,file="nagelkerkeR2_simulation.rds")
#saveRDS(S_vh,file="S_vh.rds")
