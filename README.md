# Estimating an adjusted risk difference in a cluster randomized trial with indivdual-level analyses
This repository contains all the functions and programs to reconduct the simulation study of our paper:

Jules Antoine Pereira Macedo, Bruno Giraudeau, for the ESCIENT collaborators. Estimating an adjusted risk difference in a cluster randomized trial with individual-level analyses. (Journal name). (DATE). (DOI)

## Abstract
In cluster randomized trials (CRTs) with a binary outcome, intervention effects are usually reported as odds ratios, but the CONSORT statement advocates reporting both a relative and an absolute intervention effect. With a simulation study, we assessed several methods to estimate a risk difference (RD) in the framework of a CRT with adjustment on both individual- and cluster-level covariates. We considered both a conditional approach (with the generalized linear mixed model [GLMM]) and a marginal approach (with the generalized estimating equation [GEE]). For both approaches, we considered the Gaussian, binomial and Poisson distributions. When considering the binomial or Poisson distribution, we used the g-computation method to estimate the RD. Convergence problems were observed with the GEE approach, especially with low intra-cluster coefficient correlation values, small number of clusters, small mean cluster size, high number of covariates and prevalences close to 0. All methods reported no bias. The Gaussian distribution with both approaches and binomial and Poisson distributions with the GEE approach had satisfactory results in estimating the standard error. Results for type I error and coverage rates were better with the GEE than GLMM approach. We recommend using the Gaussian distribution because of its ease of use (the RD is estimated in one step only). The GEE approach should be preferred and replaced with the GLMM approach in cases of convergence problems.

## Getting Started
### Dependencies
Data generation and analyses were run on a computing center of the Orléans-Tours university with R version 4.0, to parallelize faster but it can be run on a regular computer.
### Installing
To run all the files, start to run the package.R files to install all the package and to be able to lunch the simulation.
### Executing programs
This simulation study has 810 scenarios.

And has been divided in 2 parts with 972 scenarios each corresponding to this 2 files: 

"Simulation 1": Non-null treatment effect

"Simulation 2": No treatment effect

Notation: %method: correspond to a method “gauss” or “G” for gaussian, “bin” or “B” for binomial and “poiss” or “P” for poisson.

For each R files take care to change your file path where you want the simulation to be saved. It is named “Base_file =” in the R files (path in line 6 in R file 1_Workspace_Sim1.R, for example).

In "Simulation 1":
We have 8 files, numeroted from 1 to 4. It needs to be run in the right order (1 → 4). And for the files 3a, 3b, 3c and 4a, 4b, 4c they can be run simultaneously

"1_Workspace_Sim1":
This R file create a R workspace. You will have to load it for each of the following R files in “Simulation 1”. All the functions used are in this file, with comments, and also this file creates “list” variables "Scenarios_%n" from 1 to 972,  which correspond to the scenarios parameters used for the files in “Simulation 1”.

!Warning! 
About this simulation study, we were faced to convergence problem with the 'gee' function from the 'geepack' package. The 'gee' function using binomial and Poisson distribution for some iterations for some scenarios could not converge to find the exchangeable correlation matrix estimation. The function was running for an infinite lapse of time and we were forced to stop the R console and restart it. So, we found a way to handle this problem by creating a function ('fun_cor_%method' cf '1_Workspace_Sim1') that opens a R session for a T lapse of time, and stops this iteration at the end of the T time to manage the possible infinite time.

"2_Data_Save_Sim1":
This R file creates and saves all the dataset (1000 dataset for each of the 972 scenarios) for the Simulation 1.
The dataset of the scenarios were generated using the seed "37250".
 
"3_%a_Ana_%method_GEE_Sim1":
This R file does the analysis for GEE approach, using parallelism, for each method on G = gaussian, B = binomial, P = poisson. Using the regular analysis function for gaussian distribution and the corrected version for the binomial and poisson distribution.

"4_%a_Ana_%method_GLMM_Sim1" :
This R file does the analysis for GEE approach, using parallelism, for each method on G = gaussian, B = binomial, P = poisson. 

In "Simulation 2":

For this section is actually the same as the "Simulation 1” section. We just generated data with no treatment effect scenarios. We have to run as the same pattern than for “Simualtion 1”. Lunching the R files from 1 to 4 and the possibility to run the subdivision files a, b and c at the same time.

## Help
Users should contact pereiramacedo@univ-tours.fr (Pereira Macedo JA) if they have issues for running the programs.
## Authors
Jules Antoine Pereira Macedo: pereiramacedo@univ-tours.fr
