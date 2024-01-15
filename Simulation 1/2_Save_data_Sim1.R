####----
rm(list = ls())

Base_file = "C:/Users/pereiramacedo/Desktop/2023_06_23_Article_1_V1/Papier/Code_github/Simulation 1"
Workspace_name = "Workspace_fun_Sim1.RData"
Workspace_file = paste(Base_file,Workspace_name,sep = "/")

load(Workspace_file)

Data_file = paste(Base_file,"/Data_itt",sep = "")


#packages ----

library(ggplot2)
library(readxl)
library(cowplot)
library(psych)
library(doRNG)
library(doParallel)
library(dplyr)
library(gee)
library(geepack)
library(spind)
library(doBy)
library(arm)
library(here)
library(geesmv)
library(matrixcalc)
library(lme4)
library(glmm)
library(glmmML)
library(MCMCglmm)
library(pbkrtest)
library(modelbased)
library(parameters)
library(hglm)
library(ggh4x)
library(ggtext)
library(readr)
library(ggpubr)
library(gt)
library(knitr)
library(kableExtra)


#### ----

registerDoParallel(cores = 4) # Number of cores used for the paralellism


itt = 1000                  # Number of iteration 
Scen = rep(1:Nb_scen)

Total_time_d = Sys.time()

for (n in Scen) {
  dir.create(paste(Data_file,"/Scenario_",n,sep=""))
  
  set.seed(seed[n])
  
  Scenario_use = get(paste("Scenario",sep = "_",n))
  
  ## Simu ----
  
  
  debut = Sys.time()
  
  
  res <- foreach(i = 1:itt,
                 # .combine = rbind, #unlist
                 #.errorhandling = "remove", #s'il y a un problème enlève la ligne de
                 .packages = c("stats","arm","gee","geepack","spind","doBy","doRNG","doParallel","dplyr","here")) %dorng% fun_para_save_data(itt=1,itt_para  = i,n,Scenario = Scenario_use,Pi_int,Pi_con,rho_z,OR_int,OR_con,Data_itt_File = Data_file)
  
  
  end = Sys.time()
  time_scen = end-debut
  print(paste("Time to generated data from Scenario",n,sep = " : "))
  print(time_scen)
  
}

Total_time_e = Sys.time()
Total_time = Total_time_e - Total_time_d 
print("Total time for generating Data :")
print(Total_time)
