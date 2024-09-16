rm(list = ls())

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

# File where the workspace with the functions is ----

Base_file = "C:/Users/pereiramacedo/Desktop/2023_06_23_Article_1_V1/Papier/Code_github/Simulation 2"
Workspace_name = "Workspace_fun_Sim2.RData"
Workspace_file = paste(Base_file,Workspace_name,sep = "/")

load(Workspace_file)

Resu_file_gee = paste(Base_file,"/data_output_gee",sep = "")
Data_file = paste(Base_file,"Data_itt",sep = "/")

# G_res_file_gee = paste(Resu_file_gee,"/Data_Output_gaussid",sep = "")
# B_res_file_gee = paste(Resu_file_gee,"/Data_Output_binlogit",sep = "")
P_res_file_gee = paste(Resu_file_gee,"/Data_Output_poisslog",sep = "")

# Vector of the scenarii ----
Scen = rep(1:Nb_scen)

Total_time_d = Sys.time()

Time = 20
Cor.FG = TRUE

# Creation for each scenarii of the Excels ----
for (n in Scen) {
  
  ## Column names for empty files ----
  
  name_poisslog = data.frame(Risk_difference_Estimate_poisslog = "Risk_difference_Estimate_poisslog",
                             Number_of_itteration_poisslog = "Number_of_itteration_poisslog",
                             LL_95 = "LL_95",
                             UL_95 = "UL_95",
                             SE_itt = "SE_itt",
                             Coverage_rate_poisslog = "Coverage_rate_poisslog",
                             Abs_Bias_iteration_poisslog = "Abs_Bias_iteration_poisslog",
                             Bias_iteration_poisslog = "Bias_iteration_poisslog",
                             Relative_bias_iteration_poisslog = "Relative_bias_iteration_poisslog",
                             itt_para = "itt_para",
                             type_warning = "type_warning")
  
  
  
  ## Empty files .csv ----
  
  write.table(name_poisslog,
              file = paste(P_res_file,"/Data_output_Poiss_log_Scenario_",n,".csv",sep = ""),
              row.names = FALSE,
              col.names = FALSE,
              sep = ";")
}

# Parallelism ----
registerDoParallel(cores = Sys.getenv('SLURM_NTASKS'))

itt = 1000

for (n in Scen) {
  
  debut = Sys.time()
  
  
  res <- foreach(i = 1:itt,
                 # .combine = rbind, #unlist
                 #.errorhandling = "remove", #s'il y a un problÃ¨me enlÃ¨ve la ligne de
                 .packages = c("stats","arm","gee","geepack","spind","doBy","doRNG","doParallel","dplyr","here","geesmv","matrixcalc")) %dorng% fun_cor_poiss(i=i,Time = Time,n=n,Base_file = Base_file,Workspace_name = Workspace_name,Data_file = Data_file,Resu_file = Resu_file,Cor.FG=Cor.FG)
  
  
  end = Sys.time()
  time_scen = end-debut
  print(paste("Time to generated data from Scenario",n,sep = " : "))
  print(time_scen)
  
}



# Iterations not done after parallelism due to no convergence ----


for (n in Scen) {
  assign(paste("Data_Poiss_S",n,sep = "_"),read.csv2(paste(P_res_file,"/Data_output_Poiss_log_Scenario_",n,sep = "",".csv")))
  
  
  
  Data_Poiss_use = get(paste("Data_Poiss_S",n,sep = "_"))
  
  index_itt_sim = na.omit(unique(Data_Poiss_use$itt_para))
  length(index_itt_sim)
  itt_sim = c(1:itt)
  itt_sim_not_done = itt_sim[-index_itt_sim]
  
  if(length(itt_sim_not_done) >= 1){
    for (i in (1:length(itt_sim_not_done))) {
      res_poisslog = data_frame(Risk_difference_Estimate_poisslog = NA,
                                Number_of_itteration_poisslog = NA,
                                LL_95 = NA,
                                UL_95 = NA,
                                SE_itt = NA,
                                Coverage_rate_poisslog = NA,
                                Abs_Bias_iteration_poisslog = NA,
                                Bias_iteration_poisslog = NA,
                                Relative_bias_iteration_poisslog = NA,
                                itt_para = itt_sim_not_done[i],
                                type_warning = "No convergence")
      
      
      
      write.table(res_poisslog,file = paste(P_res_file,paste("Data_output_Poiss_log_Scenario_",sep = "",n,".csv"),sep = "/"),append = TRUE,row.names = FALSE,col.names = FALSE,sep = ";",dec = ",")
    }
  }
}

Total_time_e = Sys.time()
Total_time = Total_time_e - Total_time_d 
print("Total time for generating Data :")
print(Total_time)
