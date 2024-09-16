rm(list = ls())


### 1) Crea_data functions ----


## Input ##

# pI               : Prevalence in Intervention arm
# pC               : Prevalence in Control arm
# k                : Number of cluster per arm
# icc              : Intra-cluster correlation coefficient 
# nb_cov           : Number of "individual" covariates
# p_cov            : Vector of prevalence of each "individual-level" covariates
# OR_cov           : Vector of Odd Ratio between each "individual-level" covariates and the outcome
# nb_cov_clus      : Number of "cluster" covariates
# p_cov_clus       : Vector of prevalence of each "cluster-level" covariates
# OR_cov_clus      : Vector of Odd Ratio between each "cluster-level" covariates and the outcome

## Output ##

# Data        : A 'data.frame' which correspond to an entire data set

Crea_Data_V12 <- function(pI,pC,k,m,icc,nb_cov,p_cov,OR_cov,nb_cov_clus,p_cov_clus,OR_cov_clus){
  if (length(p_cov) < nb_cov || length(OR_cov) < nb_cov || length(p_cov_clus) < nb_cov_clus || length(OR_cov_clus) < nb_cov_clus ){stop('Vector of covariates probability or OR vector have lower size than the number of covariates')}
  
  ## Number of subject per cluster
  ## Suj_per_num[i,j] ~ Binom Negative(m,0.5)
  Suj_per_clus = matrix(0,nrow = k,ncol = 2)
  for (l in 1:(2*k)) {
    Suj_per_clus[l] = rnbinom(n = 1,mu = m, size = m**2/((0.6*m)**2-m))           # 0.6 is the coefficient of variation set a priori denoted 'cv' in the article
    while (Suj_per_clus[l]==0) {
      Suj_per_clus[l] = rnbinom(n = 1,mu = m, size = m**2/((0.6*m)**2-m))         # 0.6 is the coefficient of variation set a priori denoted 'cv' in the article
    }
  }
  
  
  ## Creation of the outcome and co variables 
  Data = c()
  Outcome = c()
  nb_indiv = sum(Suj_per_clus)
  if(nb_cov != 0){Covariate = matrix(0,ncol = nb_cov,nrow = nb_indiv)}
  if(nb_cov_clus != 0){Covariate_clus = matrix(0,nrow = nb_indiv,ncol = nb_cov_clus)}
  
  ## Cluster Effect
  Z_ij_int = rbinom(n = k,1,pI)
  Z_ij_con = rbinom(n = k,1,pC)
  
  Z_ij = matrix(c(Z_ij_int,Z_ij_con),nrow = k,ncol = 2)
  
  for (l in 1:(2*k)) {
    if(l <= k){p=pI}else{p=pC}
    n_ij = Suj_per_clus[l]
    
    ## Outcome : Y = (1-U)V + UZ
    ## U ~ Binom(1,sqrt(icc))
    ## V ~ Binom(1,p)  p = pI or pC (depend on which arm the individual is)
    ## Z ~ Binom(1,p)  same for the individuals within a cluster
    ## Y ~ Binom(1,p)
    
    U_ij = rbinom(n_ij,1,sqrt(icc))
    V_ij = rbinom(n_ij,1,p)
    
    Y_ij = (1-U_ij)*V_ij + U_ij*Z_ij[l]
    Outcome = c(Outcome,Y_ij)
    
    
    ## Covariates i :  X_i = A_i*W_i
    ## A_ijl ~ binom(1,alpha_i)
    ## alpha_i = p_cov[i]/p  ///  p = pI or pC (depend on which arm the individual is)
    ## W_ijl = (1-R_ijl)*T_ijl + R_ijl*Y
    ## R_ijl ~ binom(1,r_i)
    ## r_i = rho_cov[i] * ( sqrt((1-p)*(1-p_cov[i])) / (sqrt(alpha_i*alpha_i+1)*(1-p)) )
    ## T_ijl ~ binom(p)
    ## Y : Outcome
    ## W_ijl ~ binom(1,p)
    
    if(nb_cov !=0){
      for (c in 1:nb_cov) {
        
        q_m = p_cov[c]
        
        # Find the good coefficient of correlation 'rho' to create the individual covariate with the correspondant OR between the covariate and the outcome
        if(OR_cov[c]!=1){
          rho <- find_rho(p,q_m,OR_cov[c]) # find 'find_rho()' in the section '2 - Correlation verification function'
          if(length(rho) == 2){
            rdm = rbinom(1,1,1/2) + 1
            rho = rho[rdm]
          }
          if(length(rho)==0){stop('No solution to the equation to find a coefficient of correlation for the individuals covariates')}
          
          c_min_max <- corr_min_max(p,pC,q_m)
          if(rho <= c_min_max[1] | rho >= c_min_max[2]){stop("Doesn't valid the hypothesis of Prentice for the correlation between two binary variables" )}
        }
        if(OR_cov[c]==1){rho = 0}
        
        
        alpha = q_m/p
        r_m = rho * ( sqrt(1-q_m) / (sqrt(alpha)*sqrt(1-p) ))
        
        R_ijl = rbinom(n_ij,1,r_m)
        T_ijl = rbinom(n_ij,1,p)
        if(l==1){
          W_ijl = (1-R_ijl)*T_ijl + R_ijl*Outcome[1:Suj_per_clus[l]]
        }else{
          W_ijl = (1-R_ijl)*T_ijl + R_ijl*Outcome[(sum(Suj_per_clus[1:(l-1)])+1):sum(Suj_per_clus[1:l])]
        }
        
        A_ijl = rbinom(n_ij,1,alpha)
        X_ijl = A_ijl*W_ijl
        if(l==1){
          Covariate[1:Suj_per_clus[l],c] = X_ijl
        }else{
          Covariate[(sum(Suj_per_clus[1:(l-1)])+1):sum(Suj_per_clus[1:l]),c] = X_ijl
        }
      }
    }
  }   
  
  if(nb_cov_clus!=0){
    for (h in 1:nb_cov_clus){
      q = p_cov_clus[h]
      
      # Aij follow a binomial distribution with parameter alpha = q/p
      # Wij follow a binomial distribution with parameter p such as W : (1-R)T + RZ
      # Rij ~ B(r) with r find with the calcul based on the OR between the outcome and the cluster level covariate
      # Tij ~ B(p)
      # Zij ~ B(p) Same as the Z variable generated for the Outcome Y
      
      if(OR_cov_clus[h]!=1){
        # Finding the parameter r
        r_int = find_r_clus(pI,q,OR_cov_clus[h],icc) # Function 'find_r_clus()' described in section '2 - Correlation verification function'
        r_con = find_r_clus(pC,q,OR_cov_clus[h],icc) # Function 'find_r_clus()' described in section '2 - Correlation verification function'
      }
      if(OR_cov_clus[h]==1){
        r_int = 0 # If the OR between the outcome and the cluster level covariate is 1, the r parameter should be 0
        r_con = 0
      }
      
      Rij_int = rbinom(n = k,size = 1,prob = r_int)
      Rij_con = rbinom(n = k,size = 1,prob = r_con)
      
      Rij = matrix(c(Rij_int,Rij_con),nrow = k,ncol = 2)
      
      Tij_int = rbinom(n = k,size = 1,prob = pI)
      Tij_con = rbinom(n = k,size = 1,prob = pC)
      
      T_ij = matrix(c(Tij_int,Tij_con),nrow = k,ncol = 2)
      
      Aij_int = rbinom(n = k,size = 1,prob = q/pI)
      Aij_con = rbinom(n = k,size = 1,prob = q/pC)
      
      Aij = matrix(c(Aij_int,Aij_con),nrow = k,ncol = 2)
      Wij = (1-Rij)*T_ij + Rij*Z_ij
      
      Xij = Aij*Wij
      
      verif_Xij = apply(Xij,MARGIN = 2,FUN = sum)
      while (verif_Xij[1]==0 | verif_Xij[1] == k | verif_Xij[2]==0 | verif_Xij[1] == k) {
        Rij_int = rbinom(n = k,size = 1,prob = r_int)
        Rij_con = rbinom(n = k,size = 1,prob = r_con)
        
        Rij = matrix(c(Rij_int,Rij_con),nrow = k,ncol = 2)
        
        Tij_int = rbinom(n = k,size = 1,prob = pI)
        Tij_con = rbinom(n = k,size = 1,prob = pC)
        
        T_ij = matrix(c(Tij_int,Tij_con),nrow = k,ncol = 2)
        
        Aij_int = rbinom(n = k,size = 1,prob = q/pI)
        Aij_con = rbinom(n = k,size = 1,prob = q/pC)
        
        Aij = matrix(c(Aij_int,Aij_con),nrow = k,ncol = 2)
        Wij = (1-Rij)*T_ij + Rij*Z_ij
        
        Xij = Aij*Wij
        
        verif_Xij = apply(Xij,MARGIN = 2,FUN = sum)
      }
      
      Covariate_clus[,h] = rep(Xij,Suj_per_clus)
    }
  }
  
  ## Assignation of the arm effect for each individuals
  Arm = c(rep(1,apply(Suj_per_clus, 2, sum)[1]),rep(0,apply(Suj_per_clus, 2, sum)[2]))
  
  ## Assignation of the cluster number for each individuals
  cluster = c()
  for (l in 1:(2*k)) {
    cluster = c(cluster,rep(l,Suj_per_clus[l]))
  }
  
  Data = data.frame(Outcome,Arm,cluster)
  if(nb_cov!=0){Data = data.frame(Data,Covariate)}
  if(nb_cov_clus!=0){Data = data.frame(Data,Covariate_clus)}
  
  return(Data)
}


## Crea_Data_V13_confusion() ## 
# This is the function to create Our data set when we include confusion
# This function correspond to the function Crea_Data_V12() but we increase the number of participant by the mean proportion 
# of participant included due to the confusion that we chose for the simulation study.


## Input ##

# pI               : Prevalence in Intervention arm
# pC               : Prevalence in Control arm
# k                : Number of cluster per arm
# icc              : Intra-cluster correlation coefficient 
# nb_cov           : Number of "individual" covariates
# p_cov            : Vector of prevalence of each "individual-level" covariates
# OR_cov           : Vector of Odd Ratio between each "individual-level" covariates and the outcome
# nb_cov_clus      : Number of "cluster" covariates
# p_cov_clus       : Vector of prevalence of each "cluster-level" covariates
# OR_cov_clus      : Vector of Odd Ratio between each "cluster-level" covariates and the outcome
# Pi_int           : Proportion of the intervention group that we will be included to add confusion
# Pi_con           : Proportion of the control group that we will be included to add confusion 

## Output ##

# Data        : A 'data.frame' which correspond to data set before selection of participant to include confusion

Crea_Data_V13_confusion <- function(pI,pC,k,m,icc,nb_cov,p_cov,OR_cov,nb_cov_clus,p_cov_clus,OR_cov_clus,Pi_int,Pi_con){
  if (length(p_cov) < nb_cov || length(OR_cov) < nb_cov || length(p_cov_clus) < nb_cov_clus || length(OR_cov_clus) < nb_cov_clus ){stop('Vector of covariates probability or OR vector have lower size than the number of covariates')}
  
  m_int = m/Pi_int # increase the number of the mean of individual in a cluster in the intervention group because when we will apply the selection to create the confusion the number of individuals will decrease 
  m_con = m/Pi_con # increase the number of the mean of individual in a cluster in the control group because when we will apply the selection to create the confusion the number of individuals will decrease 
  
  ## Number of subject per cluster
  ## Suj_per_num[i,j] ~ Binom Negative(m,0.5)
  Suj_per_clus = matrix(0,nrow = k,ncol = 2)
  cv_int = 0.47
  cv_con = 0.47
  for (l in 1:(2*k)) {
    if(l<=k){m=m_int;cv=cv_int}else{m=m_con;cv=cv_con}
    Suj_per_clus[l] = rnbinom(n = 1,mu = m, size = m**2/((cv*m)**2-m))
    while (Suj_per_clus[l]==0) {
      Suj_per_clus[l] = rnbinom(n = 1,mu = m, size = m**2/((cv*m)**2-m))
    }
  }
  
  
  ## Creation of the outcome and co variables 
  Data = c()
  Outcome = c()
  nb_indiv = sum(Suj_per_clus)
  if(nb_cov != 0){Covariate = matrix(0,ncol = nb_cov,nrow = nb_indiv)}
  if(nb_cov_clus != 0){Covariate_clus = matrix(0,nrow = nb_indiv,ncol = nb_cov_clus)}
  
  ## Cluster Effect
  Z_ij_int = rbinom(n = k,1,pI)
  Z_ij_con = rbinom(n = k,1,pC)
  
  Z_ij = matrix(c(Z_ij_int,Z_ij_con),nrow = k,ncol = 2)
  
  for (l in 1:(2*k)) {
    if(l <= k){p=pI}else{p=pC}
    n_ij = Suj_per_clus[l]
    
    ## Outcome : Y = (1-U)V + UZ
    ## U ~ Binom(1,sqrt(icc))
    ## V ~ Binom(1,p)  p = pI or pC (depend on which arm the individual is)
    ## Z ~ Binom(1,p)  same for individuals in one cluster
    ## Y ~ Binom(1,p)
    
    U_ij = rbinom(n_ij,1,sqrt(icc))
    V_ij = rbinom(n_ij,1,p)
    
    Y_ij = (1-U_ij)*V_ij + U_ij*Z_ij[l]
    Outcome = c(Outcome,Y_ij)
    
    
    ## Covariates i :  X_i = A_i*W_i
    ## A_ijl ~ binom(1,alpha_i)
    ## alpha_i = p_cov[i]/p  ///  p = pI or pC (depend on which arm the individual is)
    ## W_ijl = (1-R_ijl)*T_ijl + R_ijl*Y
    ## R_ijl ~ binom(1,r_i)
    ## r_i = rho_cov[i] * ( sqrt((1-p)*(1-p_cov[i])) / (sqrt(alpha_i*alpha_i+1)*(1-p)) )
    ## T_ijl ~ binom(p)
    ## Y : Outcome
    ## W_ijl ~ binom(1,p)
    
    if(nb_cov !=0){
      for (c in 1:nb_cov) {
        
        q_m = p_cov[c]
        
        # Find the good coefficient of correlation 'rho' to create the individual covariate with the correspondant OR between the covariate and the outcome
        if(OR_cov[c]!=1){
          rho <- find_rho(p,q_m,OR_cov[c]) # find 'find_rho()' in the section '2 - Correlation verification function'
          if(length(rho) == 2){
            rdm = rbinom(1,1,1/2) + 1
            rho = rho[rdm]
          }
          if(length(rho)==0){stop('No solution to the equation to find a coefficient of correlation for the individuals covariates')}
          
          c_min_max <- corr_min_max(p,pC,q_m)
          if(rho <= c_min_max[1] | rho >= c_min_max[2]){stop("Doesn't valid the hypothesis of Prentice for the correlation between two binary variables" )}
        }
        if(OR_cov[c]==1){rho = 0}
        
        
        alpha = q_m/p
        r_m = rho * ( sqrt(1-q_m) / (sqrt(alpha)*sqrt(1-p) ))
        
        R_ijl = rbinom(n_ij,1,r_m)
        T_ijl = rbinom(n_ij,1,p)
        if(l==1){
          W_ijl = (1-R_ijl)*T_ijl + R_ijl*Outcome[1:Suj_per_clus[l]]
        }else{
          W_ijl = (1-R_ijl)*T_ijl + R_ijl*Outcome[(sum(Suj_per_clus[1:(l-1)])+1):sum(Suj_per_clus[1:l])]
        }
        
        A_ijl = rbinom(n_ij,1,alpha)
        X_ijl = A_ijl*W_ijl
        if(l==1){
          Covariate[1:Suj_per_clus[l],c] = X_ijl
        }else{
          Covariate[(sum(Suj_per_clus[1:(l-1)])+1):sum(Suj_per_clus[1:l]),c] = X_ijl
        }
      }
    }
  }   
  
  if(nb_cov_clus!=0){
    for (h in 1:nb_cov_clus){
      q = p_cov_clus[h]
      
      # Aij follow a binomial distribution with parameter alpha = q/p
      # Wij follow a binomial distribution with parameter p such as W : (1-R)T + RZ
      # Rij ~ B(r) with r find with the calcul based on the OR between the outcome and the cluster level covariate
      # Tij ~ B(p)
      # Zij ~ B(p) Same as the Z variable generated for the Outcome Y
      
      if(OR_cov_clus[h]!=1){
        # Finding parameter r to simulate Rij
        r_int = find_r_clus(pI,q,OR_cov_clus[h],icc) # Function 'find_r_clus()' described in section '2 - Correlation verification function'
        r_con = find_r_clus(pC,q,OR_cov_clus[h],icc) # Function 'find_r_clus()' described in section '2 - Correlation verification function'
      }
      if(OR_cov_clus[h]==1){
        r_int = 0 # If the OR between the outcome and the cluster level covariate is 1, the r parameter should be 0
        r_con = 0
      }
      
      Rij_int = rbinom(n = k,size = 1,prob = r_int)
      Rij_con = rbinom(n = k,size = 1,prob = r_con)
      
      Rij = matrix(c(Rij_int,Rij_con),nrow = k,ncol = 2)
      
      Tij_int = rbinom(n = k,size = 1,prob = pI)
      Tij_con = rbinom(n = k,size = 1,prob = pC)
      
      T_ij = matrix(c(Tij_int,Tij_con),nrow = k,ncol = 2)
      
      Aij_int = rbinom(n = k,size = 1,prob = q/pI)
      Aij_con = rbinom(n = k,size = 1,prob = q/pC)
      
      Aij = matrix(c(Aij_int,Aij_con),nrow = k,ncol = 2)
      Wij = (1-Rij)*T_ij + Rij*Z_ij
      
      Xij = Aij*Wij
      
      verif_Xij = apply(Xij,MARGIN = 2,FUN = sum)
      while (verif_Xij[1]==0 | verif_Xij[1] == k | verif_Xij[2]==0 | verif_Xij[1] == k) {
        Rij_int = rbinom(n = k,size = 1,prob = r_int)
        Rij_con = rbinom(n = k,size = 1,prob = r_con)
        
        Rij = matrix(c(Rij_int,Rij_con),nrow = k,ncol = 2)
        
        Tij_int = rbinom(n = k,size = 1,prob = pI)
        Tij_con = rbinom(n = k,size = 1,prob = pC)
        
        T_ij = matrix(c(Tij_int,Tij_con),nrow = k,ncol = 2)
        
        Aij_int = rbinom(n = k,size = 1,prob = q/pI)
        Aij_con = rbinom(n = k,size = 1,prob = q/pC)
        
        Aij = matrix(c(Aij_int,Aij_con),nrow = k,ncol = 2)
        Wij = (1-Rij)*T_ij + Rij*Z_ij
        
        Xij = Aij*Wij
        
        verif_Xij = apply(Xij,MARGIN = 2,FUN = sum)
      }
      
      Covariate_clus[,h] = rep(Xij,Suj_per_clus)
    }
  }
  
  ## Assignation of the arm effect for each individuals
  Arm = c(rep(1,apply(Suj_per_clus, 2, sum)[1]),rep(0,apply(Suj_per_clus, 2, sum)[2]))
  
  ## Assignation of the cluster number for each individuals
  cluster = c()
  for (l in 1:(2*k)) {
    cluster = c(cluster,rep(l,Suj_per_clus[l]))
  }
  
  Data = data.frame(Outcome,Arm,cluster)
  if(nb_cov!=0){Data = data.frame(Data,Covariate)}
  if(nb_cov_clus!=0){Data = data.frame(Data,Covariate_clus)}
  
  return(Data)
}

### 2) Correlation verification functions ----

# corr_min_max function
#Input :
# p1_i  : prevalence of the outcome in the intervention arm
# p1_c  : prevalence of the outcome in the intervention arm
# p2    : prevalence of the covariate

#Output :
# vector : vector of length 2, with the minimum value and maximum value of the correlation acceptable between the outcome and the covariate taking account the 
#          difference of prevalence in the two groups. Following Prentice article.

corr_min_max <- function(p1_i,p1_c,p2){
  q1_i = 1-p1_i
  q1_c = 1-p1_c
  q2 = 1-p2
  
  # Calculating for the prevalence of the outcome equals to the prevalence of the intervention arm 
  max_cor_i = min(sqrt((p1_i*q2)/(q1_i*p2)),sqrt((p2*q1_i)/(q2*p1_i)))
  min_cor_i = max(-sqrt((q1_i*q2)/(p1_i*p2)),-sqrt((p1_i*p2)/(q1_i*q2)))
  
  
  # Calculating for the prevalence of the outcome equals to the prevalence of the control arm 
  max_cor_c = min(sqrt((p1_c*q2)/(q1_c*p2)),sqrt((p2*q1_c)/(q2*p1_c)))
  min_cor_c = max(-sqrt((q1_c*q2)/(p1_c*p2)),-sqrt((p1_c*p2)/(q1_c*q2)))
  
  # Looking the max and min value of the interval for both arm conditions 
  min_cor_tot = max(min_cor_c,min_cor_i)
  max_cor_tot = min(max_cor_c,max_cor_i)
  
  return(c(min_cor_tot,max_cor_tot))
  
}

# fin_rho function
#Inpout :
# p  : prevalence of the outcome
# q  : prevalence of the individual-level covariate
# OR : odd ratio between the outcome and the covariate

#Output :
# rho : coefficient of correlation needed to create the covariate 

find_rho <- function(p,q,OR){
  ## Fin rho by using this formula of OR:
  
  # A = q*r_m+q*p-q*p*r_m
  # B = q - r_m*q - p*q + p*q*r_m
  # C = p-p*q-q*r_m+p*q*r_m
  # D = 1-q+q*r_m-p+p*q-p*q*r_m
  # 
  # OR = (A*D)/(B*C)
  
  ## OR = num / den      ==>    OR * den - num = 0
  ## with:
  # num = r_m² * num1 + r_m * num2 + num3
  # den = r_m² * den1 + r_m * den2 + den3
  # r_m = rho * coeff_rm
  
  num1 = q**2 - 2*p*(q**2) + (p**2)*(q**2)
  
  num2 = q - q**2 - 2*p*q + 3*p*(q**2) - 2*(p**2)*(q**2) + (p**2)*q
  
  num3 = p*q - p*(q**2) - (p**2)*q + (p**2)*(q**2)
  
  den1 = q**2 - 2*p*(q**2) + (p**2)*(q**2)
  
  den2 = - q**2 - p*q + 3*p*(q**2) - 2*(p**2)*(q**2) + (p**2)*q
  
  den3 = p*q - p*(q**2) - (p**2)*q + (p**2)*(q**2)
  
  coeff_rm = sqrt(1-q)/sqrt((q/p)*(1-p))
  
  # coefficient of the quadratic equation : a, b et c
  
  a = (coeff_rm**2) * (den1*OR - num1)
  b = coeff_rm * (den2*OR - num2)
  c = den3*OR - num3
  
  # delta of the quadratic equation 
  
  delta = b**2-4*a*c
  
  # 2 solutions of the quadratic equation
  
  x1 = (-b + sqrt(delta))/(2*a)
  x2 = (-b - sqrt(delta))/(2*a)
  
  rho = c(x1,x2)
  # correlation coefficient in the interval [-1;1]
  rho = rho[which(rho <= 1 & rho >= -1)]
  return(rho)
}

# find_r_clus
#Inpout :
# p     : prevalence of the outcome (corresponding to the prevalence of the outcome in each arm so different between two individual of two different arm intervention)
# q     : prevalence of the cluster-level covariate 
# OR    : Odd ratio between the cluster-level covariate and the outcome

#Output :
# r : Value of the parameter of the Rij variable to simulate the cluster-level covariate X_ij

find_r_clus <- function(p,q,OR_clus,icc){
  ## Find r by resolving this equation:
  
  # A = p**2*alpha +p*alpha*r*sqrt(icc) + p**2*alpha*r*sqrt(icc)
  # B = p*alpha - p*alpha*r*sqrt(icc) - p**2*alpha + p**2*alpha*r*sqrt(icc)
  # C = p - p*alpha*r*sqrt(icc) - p**2*alpha + p**2*alpha*r*sqrt(icc)
  # D = 1 - p - p*alpha + p*alpha*r*sqrt(icc) + p**2*alpha - p**2*alpha*r*sqrt(icc)
  # 
  # OR = (A*D)/(B*C)
  
  ## OR = num / den      ==>    OR * den - num = 0
  ## with:
  # num = r_m² * num1 + r_m * num2 + num3
  # den = r_m² * den1 + r_m * den2 + den3
  
  
  num1 = icc*(q**2) - 2*p*(q**2)*icc + (p**2)*(q**2)*icc
  
  num2 = q*sqrt(icc) - sqrt(icc)*(q**2) - 2*p*q*sqrt(icc) + 3*p*(q**2)*sqrt(icc) + (p**2)*q*sqrt(icc) - 2*(p**2)*(q**2)*sqrt(icc) 
  
  num3 = p*q - p*(q**2) - (p**2)*q + (p**2)*(q**2)
  
  den1 = icc * (q**2) - 2 * p * (q**2) * icc + (p**2)*(q**2)*icc
  
  den2 = - (q**2)*sqrt(icc) - p*q*sqrt(icc) + 3*p*(q**2)*sqrt(icc) + (p**2)*q*sqrt(icc) - 2*(p**2)*(q**2)*sqrt(icc)
  
  den3 = p*q - p*(q**2) - (p**2)*q + (p**2)*(q**2)
  
  
  # coefficient of the quadratic equation: a, b et c
  
  a = (den1*OR_clus - num1)
  b = (den2*OR_clus - num2)
  c = den3*OR_clus - num3
  
  # delta of the quadratic equation
  
  delta = b**2-4*a*c
  
  # 2 solutions of the quadratic equation
  if(delta>0){
    x1 = (-b + sqrt(delta))/(2*a)
    x2 = (-b - sqrt(delta))/(2*a)
    x = c(x1,x2)
  }
  if(delta==0){x = (-b)/(2*a)}
  if(delta<0){error(stop('No solution to the equation '))}
  
  # correlation coefficient in the interval [-1;1]
  
  x = x[which(x <= 1 & x >= 0)]
  if(length(x)==0){stop('No Solution in [0,1]')}
  return(x)
}

### 3) Prediction functions ---- 

#Prediction function because the 'predict()' function doesn't work with fit of the function 'geese' of the 'geepack' package.


#fit : model of the 'geese' function (package 'geepack')
#data : our data set giving by the function "Crea_data_V12()" or "fun_Data_confusion_V7()" (see below section "7-Confusion function")
#nb_cov_tot : nb_cov + nb_cov_clus = number of covariates in 'data' (individual and cluster level)

Pred_Trt_fun <- function(fit,data,nb_cov_tot){
  coef = fit$coefficients
  Pred_Trt = rep(coef[1]+ coef[2],nrow(data)) 
  if(nb_cov_tot>0){
    for (i in 1:nb_cov_tot) {
      Pred_Trt = Pred_Trt + data[,3+i]*coef[2+i]
    }
  }
  
  return(Pred_Trt)
}

Pred_NoTrt_fun <- function(fit,data,nb_cov_tot){
  coef = fit$coefficients
  Pred_NoTrt = rep(coef[1],nrow(data))
  if(nb_cov_tot>0){
    for (i in 1:nb_cov_tot) {
      Pred_NoTrt = Pred_NoTrt + data[,3+i]*coef[2+i]
    }
  }
  
  return(Pred_NoTrt)
}

### 4) get_data function ----

# get_data() : function to extract the interest variable of the model to calculate more easily the standard errors and Jacobians matrix in the deriv functions

#Input
# data : our data set giving by the function "Crea_data_V12()" or "fun_Data_confusion_V7()" (see below section "7-Confusion function")
# nb_cov: number of individual level covariates in 'data'
# nb_cov_clus = number of cluster level covariates in 'data'
# Output
# x: Matrix to facilitate calculation of SE and Jacobian matrix in following function section '5- Deriv function'


get_data <- function(data,nb_cov,nb_cov_clus){
  nb_cov_tot = nb_cov + nb_cov_clus
  intercept = rep(1,nrow(data))
  x = data.frame(intercept,data$Arm)
  if(nb_cov_tot>0){
    for (i in 1:nb_cov_tot) {
      x = data.frame(x,data[,3+i])
    }
  }
  return(as.matrix(x))
}
### 5) Deriv functions ----

# fun_deriv_logit_tot : To calculate the Jacobian matrix for the Delta method for the SE when we use the Binomial method
# fun_deriv_log_tot : To calculate the Jacobian matrix for the Delta method for the SE when we use the Poisson method
# dervi_logit : Derivative function of the logit

# Input

# fit : model of the 'geese' function 
# nb_cov : number of individual level covariates of the data set
# nb_cov_tot : number of cluster level covariates of the data set
# data = data set generated by the function "Crea_data_V12" or "fun_Data_confusion_V7()" (see below section "7-Confusion function")

# Output 

# d: The resulting matrix calculating for the delta method corresponding with the link function logit or log.


fun_deriv_logit_tot <- function(fit,nb_cov,nb_cov_clus,data){
  d = c()
  nb_cov_tot = nb_cov + nb_cov_clus
  x = get_data(data,nb_cov,nb_cov_clus)
  x1 <- x0 <- x
  x1[,2]<-1
  x0[,2]<-0
  
  betas <- fit$coefficients
  nb_var = length(betas)
  Trt = betas[1] + betas[2]
  NoTrt = betas[1]
  if(nb_cov_tot>0){
    for (i in 1:nb_cov_tot) {
      Trt = Trt + betas[2+i]*x[,2+i]            # it is the estimation of each individual considering in the intervention group
      NoTrt = NoTrt + betas[2+i]*x[,2+i]        # it is the estimation of each individual considering in the control group
    }
    for (j in 1:nb_var) {
      d_j = mean(deriv_logit(Trt) * x1[,j] -  deriv_logit(NoTrt) * x0[,j]) # Difference between intervention and control estimation for each individual as the function used in the delta method
      d=c(d,d_j)
    }
  }
  else{
    for (j in 1:nb_var) {
      d_j = mean(deriv_logit(Trt) * x1[,j] -  deriv_logit(NoTrt) * x0[,j])
      d=c(d,d_j)
    }
  }
  return(d)
}

fun_deriv_log_tot <- function(fit,nb_cov,nb_cov_clus,data){
  d = c()
  nb_cov_tot = nb_cov + nb_cov_clus
  x = get_data(data,nb_cov,nb_cov_clus)
  x1 <- x0 <- x
  x1[,2]<-1
  x0[,2]<-0
  
  betas <- fit$coefficients
  nb_var = length(betas)
  Trt = betas[1] + betas[2]
  NoTrt = betas[1]
  if(nb_cov_tot>0){
    for (i in 1:nb_cov_tot) {
      Trt = Trt + betas[2+i]*x[,2+i]          # it is the estimation of each individual considering in the intervention group
      NoTrt = NoTrt + betas[2+i]*x[,2+i]      # it is the estimation of each individual considering in the control group
    }
    for (j in 1:nb_var) {
      d_j = mean(exp(Trt) * x1[,j] -  exp(NoTrt) * x0[,j])
      d=c(d,d_j)
    }
  }
  else{
    for (j in 1:nb_var) {
      d_j = mean(exp(Trt) * x1[,j] -  exp(NoTrt) * x0[,j])
      d=c(d,d_j)
    }
  }
  return(d)
}

# deriv_logit: is the function corresponding to the derivate function of the inversev logit function 

deriv_logit <- function(x){
  exp(-x)/(1+exp(-x))^2
}

### 6) Auto formula for gee functions ----

# Automating the formula function by the number of covariate in our GEE model.
# Input
# nb_cov : number of individual level covariate of the data set.
# nb_cov_clus : number of cluster level covariate of the data set.
#
# Output
# Formula : An object at the formula type, getting the formula to use for gee function

fun_formula_for_gee <- function(nb_cov,nb_cov_clus){
  
  Covariable = ''
  
  #### Case nb_cov = 0 ----
  if(nb_cov==0){
    if(nb_cov_clus==0){
      Covariable = ''
    }
    if(nb_cov_clus==1){
      Covariable = paste(Covariable,'Covariate_clus',sep = '+')
    }
    if(nb_cov_clus > 1){
      for (i in 1:nb_cov_clus) {
        a = paste('X',i,sep = "")
        Covariable = paste(Covariable,a,sep = '+')
      }
    }
  }
  
  #### Case nb_cov = 1 ----
  if(nb_cov==1){
    Covariable = '+ Covariate'
    if(nb_cov_clus==1){
      Covariable = paste(Covariable,'Covariate_clus',sep = '+')
    }
    if(nb_cov_clus > 1){
      for (i in 1:(nb_cov_clus)) {
        a = paste('X',i,sep = "")
        Covariable = paste(Covariable,a,sep = '+')
      }
    }
  }
  
  #### Case nb_cov > 1 ----
  if(nb_cov > 1){
    for (i in 1:nb_cov) {
      a = paste('X',i,sep = "")
      Covariable = paste(Covariable,a,sep = '+')
    }
    if(nb_cov_clus == 0){
      Covariable = Covariable
    }
    if(nb_cov_clus == 1){
      Covariable = paste(Covariable,'Covariate_clus',sep = '+')
    }
    if(nb_cov_clus > 1){
      if(nb_cov_clus <= nb_cov){
        for (i in 1:(nb_cov_clus)) {
          a = paste('X',i,'.1',sep = "")
          Covariable = paste(Covariable,a,sep = '+')
        }
      }
      if(nb_cov_clus > nb_cov){
        for (i in 1:nb_cov) {
          a = paste('X',i,'.1',sep = "")
          Covariable = paste(Covariable,a,sep = '+')
        }
        for (i in (nb_cov+1):nb_cov_clus) {
          a = paste('X',i,sep = "")
          Covariable = paste(Covariable,a,sep = '+')
        }
      }
    }
  }
  
  
  Covariable = as.factor(Covariable)
  Formula = as.formula(paste("Outcome ~ Arm ",Covariable))
  
  return(Formula)
}

### 7) Confusion function ----

#Aim -> Include confusion in the data by selection in function of the covariates (Method developed by Leyrat) 

#Input

# data        : Data set (generated by the 'Crea_data_V13_confusion()' function)
# k           : Number of cluster per arm
# Pi_int      : Proportion of the intervention group that we will include to add confusion
# Pi_con      : Proportion of the control group that we will include to add confusion 
# rho_z       : Intraclass correlation coefficient for inclusion
# nb_cov      : Number of individual level covariate
# nb_cov_clus : Number of cluster level covariate
# OR_int      : Odd ratio fixed 'a priori' to create our variable of selection in the intervention group
# OR_con      : Odd ratio fixed 'a priori' to create our variable of selection in the control group

#Output

#data_fin     : Data set with confusion created by the selection of the individuals.

fun_Data_confusion_V7 <- function(data,pI,pC,p_cov,k,Pi_int,Pi_con,rho_z,nb_cov,OR_int,OR_con){
  
  OR_int = OR_int[1:nb_cov]
  OR_con = OR_con[1:nb_cov]
  p_cov = p_cov[1:nb_cov]
  
  # create our parameter (a,b) of the variable Tau_ij following a beta distribution Beta(a,b) in each group
  a_int = Pi_int*(1-rho_z)/rho_z
  b_int = (1-Pi_int)*(1-rho_z)/rho_z
  
  a_con = Pi_con*(1-rho_z)/rho_z
  b_con = (1-Pi_con)*(1-rho_z)/rho_z
  
  
  # Calculate Mean and Variance of the variable Tau_ij in each group
  
  E_tau_ij_int = a_int/(a_int+b_int)
  v_tau_ij_int = a_int*b_int/((a_int+b_int)**2*(a_int+b_int+1))
  
  E_tau_ij_con = a_con/(a_con+b_con)
  v_tau_ij_con = a_con*b_con/((a_con+b_con)**2*(a_con+b_con+1))
  
  
  ## Coefficient beta from de variable L_ijl
  
  Beta_pi_int = log(OR_int)
  Beta_pi_con = log(OR_con)
  
  
  # Generate the inclusion rate in each cluster of each arm
  Tau_ij_int = rbeta(k,a_int,b_int)
  Tau_ij_con = rbeta(k,a_con,b_con)
  
  # Cluster effect corrected
  
  eff_clus_ij_int = log(Tau_ij_int/(1-Tau_ij_int)) 
  eff_clus_ij_con = log(Tau_ij_con/(1-Tau_ij_con)) 
  
  # mean(eff_clus_ij_int)
  # mean(eff_clus_ij_con)
  
  
  # Cluster effect corrected
  
  
  eff_clus_ij = c(eff_clus_ij_int,eff_clus_ij_con)
  
  # Create a vector of length the number of individuals in the data set with the value of the cluster effect for each individual corresponding to their clusters 'appartenance'
  eff_clus_ij_vec = c()
  n = length(data$Outcome) # Total number of individual in the data set
  for (i in 1:(2*k)) {
    n_ij = length(which(data$cluster==i))
    eff_clus_ij_vec = c(eff_clus_ij_vec,rep(eff_clus_ij[i],n_ij))
  }
  # length(eff_clus_ij_vec)
  
  ## L function : corresponding to the regression function to generate a probability to the individual to be include in function of their covariates and cluster effect.  
  
  L_ijl = matrix(0,nrow = n,ncol = 1)
  
  if(nb_cov>1){
    for (i in 1:n) {
      if(data[i,]$Arm==1){
        for (j in 1:nb_cov) {
          L_ijl[i] = L_ijl[i] + Beta_pi_int[j]*(data[i,3+j]-p_cov[j]) 
        }
      }else{
        for (j in 1:nb_cov) {
          L_ijl[i] = L_ijl[i] + Beta_pi_con[j]*(data[i,3+j]-p_cov[j])
        }
      }
    }
  }
  
  L_ijl = L_ijl + eff_clus_ij_vec
  
  
  # Generate the probability for each individual to be include in the trial
  
  e_z_ijl = 1/(1+exp(-L_ijl))
  
  e_z_ijl_int = e_z_ijl[which(data$Arm==1)] 
  # mean(e_z_ijl_int)
  e_z_ijl_con = e_z_ijl[which(data$Arm==0)] 
  # mean(e_z_ijl_con)
  
  e_z_ijl = c(e_z_ijl_int,e_z_ijl_con)
  
  # e_z_ijl is a probability so is in [0,1]
  # So correction will be :
  max(e_z_ijl)
  min(e_z_ijl)
  
  index_need_correction = which(e_z_ijl > 1 )
  e_z_ijl[index_need_correction] = 1
  
  
  
  # Generate the inclusion variable Z for each individuals (Z_ijl = 1 -> include / Z_ijl = 0 -> No include)
  
  Z_ijl=c()
  for (i in 1:length(e_z_ijl)) {
    Z_ijl[i] = rbinom(1,1,e_z_ijl[i])
  }
  # length(Z_ijl)
  Inclus = which(Z_ijl==1)
  
  data_fin = data[Inclus,]
  
  while(length(unique(data_fin$cluster)) != (2*k)){ # while conditions in case if one cluster did not include at least 1 individual 
    # Generate the inclusion variable Z for each individuals (Z_ijl = 1 -> include / Z_ijl = 0 -> No include)
    num_clus = c(1:(2*k))
    
    num_clus_manquant = num_clus[-unique(data_fin$cluster)]
    
    for (i in num_clus_manquant) {
      index = which(data$cluster==i)
      for (j in index) {
        Z_ijl[j] = rbinom(1,1,e_z_ijl[j])
      }
    }
    
    
    # Z_ijl=c()
    # for (i in 1:length(e_z_ijl)) {
    #   Z_ijl[i] = rbinom(1,1,e_z_ijl[i])
    # }
    # length(Z_ijl)
    Inclus = which(Z_ijl==1)
    
    data_fin = data[Inclus,]
  }
  
  
  return(data_fin)
}


### 8) Saving data sets ----

# fun_para_save_dat
# function to generate and save each iteration and each data set into a Excel (.csv) file.
# This function is used with 'parallel' package, we used parallelism to generate our datas but we can use it without parallelism

#Input
# itt           : Number of iteration for this function (!care if you use parallelisation itt will be =1)
# itt_para      : It is the iteration that is performing during the parallelisation
# n             : The scenario number
# Scenario      : The scenario_use, (list type, defined in below in this file)
# Pi_int        : Proportion of inclusion in the intervention arm
# Pi_con        : Proportion of inclusion in the control arm
# rho_z         : Intraclass correlation coefficient for inclusion
# OR_int        : Odd ratio fixed 'a priori' to create our variable of selection in the intervention group
# OR_con        : Odd ratio fixed 'a priori' to create our variable of selection in the control group    
# Data_itt_File : Name of the path where you want your data save


fun_para_save_data <- function(itt,itt_para,n,Scenario,Pi_int,Pi_con,rho_z,OR_int,OR_con,Data_itt_File){
  
  
  for (i in 1:itt) {
    
    # Taking the differrent parameters set to the scenario number n
    
    pI = as.numeric(Scenario[1])                      # Prevalence in intervention arm
    pC = as.numeric(Scenario[2])                      # Prevalence in control arm
    k = as.numeric(Scenario[3])                       # Number of cluster per arm
    m = as.numeric(Scenario[4])                       # Mean cluster size
    icc = as.numeric(Scenario[5])                     # Intracluster Correlation Coefficient
    nb_cov = as.numeric(Scenario[6])                  # Number of individual level covariate
    nb_cov_clus = as.numeric(Scenario[7])             # Number of cluster level covariate
    p_cov = as.vector(unlist(Scenario[8]))            # Prevalence of the individual level covariate
    p_cov_clus = as.vector(unlist(Scenario[9]))       # Prevalence of the cluster level covariate
    OR_cov = as.vector(unlist(Scenario[10]))          # OR to measure the association between individual level covariate and the outcome
    OR_cov_clus = as.vector(unlist(Scenario[11]))     # OR to measure the association between cluster level covariate and the outcome
    
    # Depending if we have an individual level covariate or not. Because if not we do not have to include confounding, so we just have to use the regular version or the data generating function
    if(nb_cov==0){
      data = Crea_Data_V12(pI,pC,k,m,icc,nb_cov,p_cov,OR_cov,nb_cov_clus,p_cov_clus,OR_cov_clus)
    }else{
      data_itt = Crea_Data_V13_confusion(pI,pC,k,m,icc,nb_cov,p_cov,OR_cov,nb_cov_clus,p_cov_clus,OR_cov_clus,Pi_int,Pi_con)
      data = fun_Data_confusion_V7(data_itt,pI,pC,p_cov,k,Pi_int,Pi_con,rho_z,nb_cov,OR_int,OR_con)
    }
    
    # Saving the data set to an Excel file (.csv)
    write.csv2(data,here::here(paste(Data_itt_File,"/Scenario_",n,sep=""),paste("/Data_itt_",sep = "",itt_para,"_scen_",n,".csv")),row.names = FALSE)
  }
  
}
### 9) Corrected Variance Fay and Graubard Functions ----

### Here we using the function GEE.var.MD from the 'geesmv' package to correct the variance for small number of cluster
### We only change the Input 'formula' -> 'fit' to not perform a second time the 'gee' function to decrease the calculation time

fun_var_corrected_FG <- function(fit,formula,id,family=gaussian,data,corstr="independence",b=0.75){
  
  #########################################################################
  # !!!! Function GEE.var.md from the 'geesmv' package only changing 'formula' by 'fit' !!!!
  # !!!! to reduce calculation time by not doing a second gee                           !!!!
  
  # Arguments:
  # fit      specify the gee fit already performed
  # formula  specify the model of interest
  # family   "gaussian", "binomial" or "poisson"
  # data     data frame
  # corstr   Working correlation structure: "independence", "AR-M", "exchangeable", "unstructured".
  # value:   GEE returns the following elements
  #          cov.var      estimate of the variance-covariance matrix for robust variance.
  
  # !!!! Function GEE.var.md from the 'geesmv' package only changing 'formula' by 'fit' !!!!
  # !!!! to reduce calculation time by not doing a second gee                           !!!!
  #########################################################################
  # Delete the records with missing data in predictors or outcomes;
  if (is.null(data$id)){
    index <- which(names(data)==id)
    data$id <- data[,index]}
  
  ### na.action: only na.omit is used for gee;
  init <- model.frame(formula, data)
  init$num <- 1:length(init[,1])
  if(any(is.na(init))){
    index <- na.omit(init)$num
    data <- data[index,]
    ### Get the design matrix;
    m <- model.frame(formula, data)
    mt <- attr(m, "terms") 
    data$response <- model.response(m, "numeric")
    mat <- as.data.frame(model.matrix(formula, m))
  }else{
    ### Get the design matrix;
    m <- model.frame(formula, data)
    mt <- attr(m, "terms") 
    data$response <- model.response(m, "numeric")
    mat <- as.data.frame(model.matrix(formula, m))
  }
  
  ### Fit the GEE model to get the estimate of parameters \hat{\beta};
  gee.fit <- fit
  beta_est <- gee.fit$coefficient
  alpha <- gee.fit$working.correlation[1,2]
  len <- length(beta_est)
  len_vec <- len^2
  
  ### Estimate the robust variance for \hat{\beta}
  data$id <- gee.fit$id
  cluster<-cluster.size(data$id)
  ncluster<-max(cluster$n)
  size<-cluster$m
  mat$subj <- rep(unique(data$id), cluster$n)
  if(is.character(corstr)){
    var <- switch(corstr,
                  "independence"=cormax.ind(ncluster),
                  "exchangeable"=cormax.exch(ncluster, alpha),
                  "AR-M"=cormax.ar1(ncluster, alpha),
                  "unstructured"=summary(gee.fit)$working.correlation,)
  }else{
    print(corstr)
    stop("'working correlation structure' not recognized")
  }   
  if(is.character(family)){
    family <- switch(family,
                     "gaussian"="gaussian",
                     "binomial"="binomial",
                     "poisson"="poisson")
  }else{ 
    if(is.function(family)){
      family <- family()[[1]]
    }else{
      print(family)
      stop("'family' not recognized")
    }    
  }
  
  cov.beta<-unstr<-matrix(0,nrow=len,ncol=len)
  step11<-matrix(0, nrow=len, ncol=len)
  for (i in 1:size){
    y<-as.matrix(data$response[data$id==unique(data$id)[i]])
    covariate<-as.matrix(subset(mat[,-length(mat[1,])], mat$subj==unique(data$id)[i]))
    ncluster=cluster$n[i]
    var1=var[1:ncluster,1:ncluster] 
    if (family=="gaussian"){ 
      Vi=gee.fit$scale*var1
      xx<-t(covariate)%*%solve(Vi)%*%covariate
      step11<-step11+xx  
    }else if (family=="poisson") {
      D<-mat.prod(covariate, exp(covariate%*%beta_est))
      Vi <- diag(sqrt(c(exp(covariate%*%beta_est))),ncluster)%*%var1%*%diag(sqrt(c(exp(covariate%*%beta_est))),ncluster)
      xx<-t(D)%*%solve(Vi)%*%D
      step11<-step11+xx
    }else if (family=="binomial"){
      D<-mat.prod(covariate, exp(covariate%*%beta_est)/((1+exp(covariate%*%beta_est))^2))
      Vi <- diag(sqrt(c(exp(covariate%*%beta_est)/(1+exp(covariate%*%beta_est))^2)),ncluster)%*%var1%*%diag(sqrt(c(exp(covariate%*%beta_est)/(1+exp(covariate%*%beta_est))^2)),ncluster)
      xx<-t(D)%*%solve(Vi)%*%D
      step11<-step11+xx 
    }
  }
  step12<-matrix(0,nrow=len,ncol=len)
  p<-matrix(0,nrow=len_vec,ncol=size)
  
  for (i in 1:size){
    y<-as.matrix(data$response[data$id==unique(data$id)[i]])
    covariate<-as.matrix(subset(mat[,-length(mat[1,])], mat$subj==unique(data$id)[i]))
    ncluster=cluster$n[i]
    var1=var[1:ncluster,1:ncluster]
    if (family=="gaussian"){ 
      ## set up the scale parameter;
      Vi=gee.fit$scale*var1
      xx<-t(covariate)%*%solve(Vi)%*%covariate
      Qi <- xx%*%solve(step11)
      Ai<-diag((1-pmin(b,diag(Qi)))^(-0.5))
      xy<-Ai%*%t(covariate)%*%solve(Vi)%*%(y-covariate%*%beta_est)
      step12<-step12+xy%*%t(xy)
      p[,i]<-vec(xy%*%t(xy))
    }else if (family=="poisson") {
      ## set up the scale parameter;
      D<-mat.prod(covariate, exp(covariate%*%beta_est))
      Vi <- diag(sqrt(c(exp(covariate%*%beta_est))),ncluster)%*%var1%*%diag(sqrt(c(exp(covariate%*%beta_est))),ncluster)
      xx<-t(D)%*%solve(Vi)%*%D
      Qi <- xx%*%solve(step11)
      Ai<-diag((1-pmin(b,diag(Qi)))^(-0.5))
      xy<-Ai%*%t(D)%*%solve(Vi)%*%(y-exp(covariate%*%beta_est))
      step12<-step12+xy%*%t(xy)
      p[,i]<-vec(xy%*%t(xy))
    }else if (family=="binomial"){
      ## set up the scale parameter;
      D<-mat.prod(covariate, exp(covariate%*%beta_est)/((1+exp(covariate%*%beta_est))^2))
      Vi <- diag(sqrt(c(exp(covariate%*%beta_est)/(1+exp(covariate%*%beta_est))^2)),ncluster)%*%var1%*%diag(sqrt(c(exp(covariate%*%beta_est)/(1+exp(covariate%*%beta_est))^2)),ncluster)
      xx<-t(D)%*%solve(Vi)%*%D
      Qi <- xx%*%solve(step11)
      Ai<-diag((1-pmin(b,diag(Qi)))^(-0.5))
      xy<-Ai%*%t(D)%*%solve(Vi)%*%(y-exp(covariate%*%beta_est)/(1+exp(covariate%*%beta_est)))
      step12<-step12+xy%*%t(xy)
      p[,i]<-vec(xy%*%t(xy)) 
    }    
  }
  cov.beta<-solve(step11)%*%(step12)%*%solve(step11)
  return(cov.beta)
}
### 10) Function simulation parallelism ----
# One function for each distribution (cf. method)

# Input : 
# data_itt    : data set generated as a Data.frame
# itt_para    : which itteration of the paralellism is analyzed (it's too know which data we use on the 1000 simulated)
# matcor_type : correlation structure used in the gee
# Scenario    : Scenario corresponding to the the data set analyzed
# Cor.FG      : Logical object (TRUE or FALSE), TRUE if we apply the Fay and Graubard correction as defined in section '9 -Corrected Variance Fay and Graubard Functions'

# Output :
# Results of an iteration (itt_para) for one method by a table with the results : (Risk difference estimated, CI at 95%, Coverage_rage,Bias ....)

fun_sim_gauss <- function(data_itt,itt_para,matcor_type = "exchangeable",Scenario = list(pI,pC,k,m,icc,nb_cov,nb_cov_clus,p_cov,p_cov_clus,OR_cov,OR_cov_clus),Cor.FG=TRUE){
  
  # 
  
  ## Scenario parameters ----
  pI = as.numeric(Scenario[1])
  pC = as.numeric(Scenario[2])
  k = as.numeric(Scenario[3])
  m = as.numeric(Scenario[4])
  icc = as.numeric(Scenario[5])
  nb_cov = as.numeric(Scenario[6])
  nb_cov_clus = as.numeric(Scenario[7])
  p_cov = as.vector(unlist(Scenario[8]))
  p_cov_clus = as.vector(unlist(Scenario[9]))
  OR_cov = as.vector(unlist(Scenario[10]))
  OR_cov_clus = as.vector(unlist(Scenario[11]))
  
  RD_th = pI-pC
  nb_cov_tot = nb_cov+nb_cov_clus
  
  ## Formula for gee ----
  
  Formula_gee = fun_formula_for_gee(nb_cov,nb_cov_clus)
  
  ## GEE  ----
  
  ## Gaussian Identity  ----
  
  # Test if we catch an error or not
  
  cap_op  <- capture.output(suppressMessages(res <- fit_gaussid <- try(gee(formula = Formula_gee,
                                                                           data = data_itt,
                                                                           id=cluster,
                                                                           corstr = matcor_type,
                                                                           family = gaussian),silent = TRUE)))
  if(inherits(res, "try-error")==FALSE){
    
    
    summary_fit_gaussid = summary(fit_gaussid)
    
    RD_gaussid = as.numeric(fit_gaussid$coefficients[2]) # Estimand of our Risk difference
    
    #Variance Betas Corrected
    
    if(Cor.FG==TRUE){
      Mat_var_cov_corrected_FG = fun_var_corrected_FG(fit = fit_gaussid,
                                                      formula = Formula_gee,
                                                      id = 'cluster',
                                                      family = gaussian,
                                                      data = data_itt,
                                                      corstr = matcor_type,
                                                      b=0.75)
      
      se_tot_gaussid = sqrt(diag(Mat_var_cov_corrected_FG)[2])
    }
    
    #Variance Uncorrected
    if(Cor.FG==FALSE){
      se_tot_gaussid = summary_fit_gaussid$coefficients[2,4]
    }
    
    
    
    # CI RD
    
    DDF = 2*k - (2 + nb_cov + nb_cov_clus) # degree of freedom as the number of cluster minus number of parameters estimated in the model (intercept + arm + number of individual level covariates + number of cluster level covariates)
    
    LL_CI_RD_gaussid = RD_gaussid - qt(0.975,DDF)*se_tot_gaussid            # Confidence interval lower limit
    # LL_CI_RD_itt_gaussid = c( LL_CI_RD_itt_gaussid, LL_CI_RD_gaussid)
    
    UL_CI_RD_gaussid = RD_gaussid + qt(0.975,DDF)*se_tot_gaussid            # Confidence interval upper limit
    # UL_CI_RD_itt_gaussid =c(UL_CI_RD_itt_gaussid,UL_CI_RD_gaussid)
    
    ## Results Gauss id ----
    
    itt_without_error_gaussid = length(RD_gaussid)
    Coverage_rate_gaussid = length(which(RD_th > LL_CI_RD_gaussid & RD_th < UL_CI_RD_gaussid))/length(RD_gaussid)*100
    Absolute_Bias_gaussid = abs(RD_th - RD_gaussid)
    Bias_gaussid = RD_gaussid - RD_th
    Relative_Bias_gaussid = RD_gaussid/RD_th
    
    res_gaussid = data_frame(Risk_difference_Estimate_gaussid = RD_gaussid,
                             Number_of_itteration_gaussid = itt_without_error_gaussid,
                             LL_95 = LL_CI_RD_gaussid,
                             UL_95 = UL_CI_RD_gaussid,
                             SE_itt = se_tot_gaussid,
                             Coverage_rate_gaussid = Coverage_rate_gaussid,
                             Abs_Bias_iteration_gaussid = Absolute_Bias_gaussid,
                             Bias_iteration_gaussid = Bias_gaussid,
                             Relative_bias_iteration_gaussid = Relative_Bias_gaussid,
                             itt_para = itt_para,
                             type_error = fit_gaussid$error)
    
    res_gaussid[unique(which(is.na(res_gaussid) | abs(res_gaussid$Risk_difference_Estimate_gaussid)>=1 | res_gaussid$type_error !=0,arr.ind = TRUE)[,1]),1:9] = NA
  }
  else{
    res_gaussid = data_frame(Risk_difference_Estimate_gaussid = NA,
                             Number_of_itteration_gaussid = NA,
                             LL_95 = NA,
                             UL_95 = NA,
                             SE_itt = NA,
                             Coverage_rate_gaussid = NA,
                             Abs_Bias_iteration_gaussid = NA,
                             Bias_iteration_gaussid = NA,
                             Relative_bias_iteration_gaussid = NA,
                             itt_para = itt_para,
                             type_error = "error")
  }
  
  return(res_gaussid)
  
}

fun_sim_bin <- function(data_itt,itt_para,matcor_type = "exchangeable",Scenario = list(pI,pC,k,m,icc,nb_cov,nb_cov_clus,p_cov,p_cov_clus,OR_cov,OR_cov_clus),Cor.FG=TRUE){
  
  ## Scenario parameters ----
  pI = as.numeric(Scenario[1])
  pC = as.numeric(Scenario[2])
  k = as.numeric(Scenario[3])
  m = as.numeric(Scenario[4])
  icc = as.numeric(Scenario[5])
  nb_cov = as.numeric(Scenario[6])
  nb_cov_clus = as.numeric(Scenario[7])
  p_cov = as.vector(unlist(Scenario[8]))
  p_cov_clus = as.vector(unlist(Scenario[9]))
  OR_cov = as.vector(unlist(Scenario[10]))
  OR_cov_clus = as.vector(unlist(Scenario[11]))
  
  RD_th = pI-pC
  nb_cov_tot = nb_cov+nb_cov_clus
  
  # initial_value_gee = c(1,0,rep(0.5,(nb_cov_tot)))
  
  ## Formula for gee ----
  
  Formula_gee = fun_formula_for_gee(nb_cov,nb_cov_clus)
  
  ## GEE  ----
  ## Binomial Logit  ----
  
  cap_op2  <- capture.output(suppressMessages(res <- fit_binlogit <- try(gee(formula = Formula_gee,
                                                                             data = data_itt,
                                                                             id=cluster,
                                                                             corstr = matcor_type,
                                                                             family = binomial),silent = TRUE)))
  if(inherits(res, "try-error")==FALSE){
    
    
    ### ALL include ----
    
    # Margin mean
    P_Trt_est_binlogit = mean(invlogit(Pred_Trt_fun(fit_binlogit,data_itt,nb_cov_tot)))
    P_NoTrt_est_binlogit = mean(invlogit(Pred_NoTrt_fun(fit_binlogit,data_itt,nb_cov_tot)))
    
    RD_binlogit = P_Trt_est_binlogit - P_NoTrt_est_binlogit
    
    #Variance Betas Corrected
    
    if(Cor.FG==TRUE){
      Mat_var_cov_corrected_FG = fun_var_corrected_FG(fit = fit_binlogit,
                                                      formula = Formula_gee,
                                                      id = 'cluster',
                                                      family = binomial,
                                                      data = data_itt,
                                                      corstr = matcor_type,
                                                      b=0.75)
      
      #Delta method to calculate the SE of the risk difference
      
      deriv_vec_binlogit = fun_deriv_logit_tot(fit_binlogit,nb_cov,nb_cov_clus,data_itt)
      variance_binlogit <- deriv_vec_binlogit %*% Mat_var_cov_corrected_FG %*% deriv_vec_binlogit
      se_tot_binlogit <- sqrt(diag(variance_binlogit))
    }
    
    if(Cor.FG==FALSE){
      #SEs
      
      deriv_vec_binlogit = fun_deriv_logit_tot(fit_binlogit,nb_cov,nb_cov_clus,data_itt)
      variance_binlogit <- deriv_vec_binlogit %*% fit_binlogit$robust.variance %*% deriv_vec_binlogit
      se_tot_binlogit <- sqrt(diag(variance_binlogit))
    }
    
    
    # CI RD 
    
    DDF = 2*k - (2 + nb_cov + nb_cov_clus)
    
    LL_CI_RD_binlogit = RD_binlogit - qt(0.975,DDF)*se_tot_binlogit
    
    UL_CI_RD_binlogit = RD_binlogit + qt(0.975,DDF)*se_tot_binlogit
    
    ## Results Bin logit ----
    
    itt_without_error_binlogit = length(RD_binlogit)
    Coverage_rate_binlogit = length(which(RD_th > LL_CI_RD_binlogit & RD_th < UL_CI_RD_binlogit))/length(RD_binlogit)*100
    Absolute_Bias_binlogit = abs(RD_th - RD_binlogit)
    Bias_binlogit = RD_binlogit - RD_th
    Relative_Bias_binlogit = RD_binlogit/RD_th
    res_binlogit = data_frame(Risk_difference_Estimate_binlogit = RD_binlogit,
                              Number_of_itteration_binlogit = itt_without_error_binlogit,
                              LL_95 = LL_CI_RD_binlogit,
                              UL_95 = UL_CI_RD_binlogit,
                              SE_itt = se_tot_binlogit,
                              Coverage_rate_binlogit = Coverage_rate_binlogit,
                              Abs_Bias_iteration_binlogit = Absolute_Bias_binlogit,
                              Bias_iteration_binlogit = Bias_binlogit,
                              Relative_bias_iteration_binlogit = Relative_Bias_binlogit,
                              itt_para = itt_para,
                              type_error = fit_binlogit$error)
    
    res_binlogit[unique(which(is.na(res_binlogit) | abs(res_binlogit$Risk_difference_Estimate_binlogit) >= 1 | res_binlogit$type_error !=0,arr.ind = TRUE)[,1]),1:9] = NA
  }
  else{
    res_binlogit = data_frame(Risk_difference_Estimate_binlogit = NA,
                              Number_of_itteration_binlogit = NA,
                              LL_95 = NA,
                              UL_95 = NA,
                              SE_itt = NA,
                              Coverage_rate_binlogit = NA,
                              Abs_Bias_iteration_binlogit = NA,
                              Bias_iteration_binlogit = NA,
                              Relative_bias_iteration_binlogit = NA,
                              itt_para = itt_para,
                              type_error = "error")
  }
  
  
  
  
  
  return(res_binlogit)
  
}

fun_sim_poiss <- function(data_itt,itt_para,matcor_type = "exchangeable",Scenario = list(pI,pC,k,m,icc,nb_cov,nb_cov_clus,p_cov,p_cov_clus,OR_cov,OR_cov_clus),Cor.FG=TRUE){
  
  ## Scenario parameters ----
  pI = as.numeric(Scenario[1])
  pC = as.numeric(Scenario[2])
  k = as.numeric(Scenario[3])
  m = as.numeric(Scenario[4])
  icc = as.numeric(Scenario[5])
  nb_cov = as.numeric(Scenario[6])
  nb_cov_clus = as.numeric(Scenario[7])
  p_cov = as.vector(unlist(Scenario[8]))
  p_cov_clus = as.vector(unlist(Scenario[9]))
  OR_cov = as.vector(unlist(Scenario[10]))
  OR_cov_clus = as.vector(unlist(Scenario[11]))
  
  RD_th = pI-pC
  nb_cov_tot = nb_cov+nb_cov_clus
  
  # initial_value_gee = c(1,0,rep(0.5,(nb_cov_tot)))
  
  ## Formula for gee ----
  
  Formula_gee = fun_formula_for_gee(nb_cov,nb_cov_clus)
  
  ## GEE  ----
  
  ## Poisson Log  ----
  
  
  cap_op2  <- capture.output(suppressMessages(res <- fit_poisslog <- try(gee(formula = Formula_gee,
                                                                             data = data_itt,
                                                                             id=cluster,
                                                                             corstr = matcor_type,
                                                                             family = poisson),silent = TRUE)))
  if(inherits(res, "try-error")==FALSE){
    
    if(fit_poisslog$error != 0){res_poisslog = data_frame(Risk_difference_Estimate_poisslog = NA,
                                                          Number_of_itteration_poisslog = NA,
                                                          LL_95 = NA,
                                                          UL_95 = NA,
                                                          SE_itt = NA,
                                                          Coverage_rate_poisslog = NA,
                                                          Abs_Bias_iteration_poisslog = NA,
                                                          Bias_iteration_poisslog = NA,
                                                          Relative_bias_iteration_poisslog = NA,
                                                          itt_para = itt_para,
                                                          type_error = fit_poisslog$error)
    
    }
    else{
      # Margin mean
      P_Trt_est_poisslog = mean(exp(Pred_Trt_fun(fit_poisslog,data_itt,nb_cov+nb_cov_clus)))
      P_NoTrt_est_poisslog = mean(exp(Pred_NoTrt_fun(fit_poisslog,data_itt,nb_cov+nb_cov_clus)))
      
      RD_poisslog = P_Trt_est_poisslog -  P_NoTrt_est_poisslog
      
      #Variance Betas Corrected
      
      if(Cor.FG==TRUE){
        Mat_var_cov_corrected_FG = fun_var_corrected_FG(fit = fit_poisslog,
                                                        formula = Formula_gee,
                                                        id = 'cluster',
                                                        family = poisson,
                                                        data = data_itt,
                                                        corstr = matcor_type,
                                                        b=0.75)
        
        #Delta method to calculate the SE of the risk difference
        
        deriv_vec_poisslog = fun_deriv_log_tot(fit_poisslog,nb_cov,nb_cov_clus,data_itt)
        variance_poisslog <- deriv_vec_poisslog %*% Mat_var_cov_corrected_FG %*% deriv_vec_poisslog
        se_tot_poisslog <- sqrt(diag(variance_poisslog))
      }
      
      if(Cor.FG==FALSE){
        #SEs
        deriv_vec_poisslog = fun_deriv_log_tot(fit_poisslog,nb_cov,nb_cov_clus,data_itt)
        variance_poisslog <- deriv_vec_poisslog %*% fit_poisslog$robust.variance %*% deriv_vec_poisslog
        se_tot_poisslog <- sqrt(diag(variance_poisslog))
      }
      
      
      # CI RD 
      
      DDF = 2*k - (2 + nb_cov + nb_cov_clus)
      
      LL_CI_RD_poisslog = RD_poisslog - qt(0.975,DDF)*se_tot_poisslog
      
      UL_CI_RD_poisslog = RD_poisslog + qt(0.975,DDF)*se_tot_poisslog
      
      ## Results Poiss log ----
      
      itt_without_error_poisslog = length(RD_poisslog)
      Coverage_rate_poisslog = length(which(RD_th > LL_CI_RD_poisslog & RD_th < UL_CI_RD_poisslog))/length(RD_poisslog)*100
      Absolute_Bias_poisslog = abs(RD_th - RD_poisslog)
      Bias_poisslog = RD_poisslog - RD_th
      Relative_Bias_poisslog = RD_poisslog/RD_th
      res_poisslog = data_frame(Risk_difference_Estimate_poisslog = RD_poisslog,
                                Number_of_itteration_poisslog = itt_without_error_poisslog,
                                LL_95 = LL_CI_RD_poisslog,
                                UL_95 = UL_CI_RD_poisslog,
                                SE_itt = se_tot_poisslog,
                                Coverage_rate_poisslog = Coverage_rate_poisslog,
                                Abs_Bias_iteration_poisslog = Absolute_Bias_poisslog,
                                Bias_iteration_poisslog = Bias_poisslog,
                                Relative_bias_iteration_poisslog = Relative_Bias_poisslog,
                                itt_para = itt_para,
                                type_error = fit_poisslog$error)
      
      res_poisslog[unique(which(is.na(res_poisslog) | abs(res_poisslog$Risk_difference_Estimate_poisslog)>=1 ,arr.ind = TRUE)[,1]),1:9] = NA
    }
  }
  
  else{
    res_poisslog = data_frame(Risk_difference_Estimate_poisslog = NA,
                              Number_of_itteration_poisslog = NA,
                              LL_95 = NA,
                              UL_95 = NA,
                              SE_itt = NA,
                              Coverage_rate_poisslog = NA,
                              Abs_Bias_iteration_poisslog = NA,
                              Bias_iteration_poisslog = NA,
                              Relative_bias_iteration_poisslog = NA,
                              itt_para = itt_para,
                              type_error = "error")
  }
  
  
  return(res_poisslog)
  
}


### 11) Function analyzing by methods ----

# fun_para_analyse_%method : Are the function used during the parallelism
# They download the data set corresponding to the scenario and iteration corresponding during parallelism 
# and stock the result of one scenario in a Excel file created in advance but we will see that in R files number 3 and 4 of the whole simulation study.

# Input 
# itt_para       : Is the iteration of used in the parallelism mechanism, to take into account of a specific iteration (In our simulation study 1 to 1000)
# n              : Scenario number used
# matcor_type    : correlation matrix structure used in GEE model (in our case exchangeable)
# Scenario_use   : List of the parameters of the scenario 'n' used
# Cor.FG         : Logical object (TRUE or FALSE), TRUE if we apply the Fay and Graubard correction as defined in section '9 -Corrected Variance Fay and Graubard Functions'
# Data_file      : Path where the data are saved
# G_res_file     : Path where the results will be save for the Gaussian distribution
# B_res_file     : Path where the results will be save for the binomial distribution
# P_res_file     : Path where the results will be save for the Poisson distribution

# Output
# Write and stock the result for one iterations ('itt_para') in a specific Excel file for the scenario 'n'

fun_para_analyse_gauss <- function(itt_para,n,matcor_type,Scenario_use,Cor.FG=TRUE,Data_file,G_res_file){
  
  # We download the data set for the 'itt_para' iteration where all the data file are stock
  data = read.csv(paste(paste(Data_file,"/Scenario_",n,sep=""),sep ="",paste("/Data_itt_",sep = "",itt_para,"_scen_",n,".csv")), sep=";")
  
  #  Then we call the function who give us the results for each iteration
  a = fun_sim_gauss(data,itt_para,matcor_type = matcor_type,Scenario = Scenario_use,Cor.FG)
  
  # Finally we stock the result of 1 iteration ('itt_para' iteration result) that we stock in an Excel file create for each scenario independently
  write.table(a,file = paste(G_res_file,paste("/Data_output_Gauss_id_Scenario_",sep = "",n,".csv"),sep = ""),append = TRUE,row.names = FALSE,col.names = FALSE,sep = ";",dec = ",")
}

fun_para_analyse_bin <- function(itt_para,n,matcor_type,Scenario_use,Cor.FG=TRUE,Data_file,B_res_file){
  
  # We download the data set for the 'itt_para' iteration where all the data file are stock
  data = read.csv(paste(paste(Data_file,"/Scenario_",n,sep=""),sep ="",paste("/Data_itt_",sep = "",itt_para,"_scen_",n,".csv")), sep=";")
  
  #  Then we call the function who give us the results for each iteration
  a = fun_sim_bin(data,itt_para,matcor_type = matcor_type,Scenario = Scenario_use,Cor.FG)
  
  # Finally we stock the result of 1 iteration ('itt_para' iteration result) that we stock in an Excel file create for each scenario independently
  write.table(a,file = paste(B_res_file,paste("/Data_output_Binom_logit_Scenario_",sep = "",n,".csv"),sep = ""),append = TRUE,row.names = FALSE,col.names = FALSE,sep = ";",dec = ",")
}

fun_para_analyse_poiss <- function(itt_para,n,matcor_type,Scenario_use,Cor.FG=TRUE,Data_file,P_res_file){
  
  # We download the data set for the 'itt_para' iteration where all the data file are stock
  data = read.csv(paste(paste(Data_file,"/Scenario_",n,sep=""),sep ="",paste("/Data_itt_",sep = "",itt_para,"_scen_",n,".csv")), sep=";")
  
  #  Then we call the function who give us the results for each iteration
  a = fun_sim_poiss(data,itt_para,matcor_type = matcor_type,Scenario = Scenario_use,Cor.FG)
  
  # Finally we stock the result of 1 iteration ('itt_para' iteration result) that we stock in an Excel file create for each scenario independently
  write.table(a,file = paste(P_res_file,paste("/Data_output_Poiss_log_Scenario_",sep = "",n,".csv"),sep = ""),append = TRUE,row.names = FALSE,col.names = FALSE,sep = ";",dec = ",")
}



# After several test we have detected some convergence problem for the 'gee' function from the 'geepack' package
# Some iterations who did not converge could be running for an infinite laps of time
# So we created a function that for each iteration work separately from the R console
# We create a batch file that we lunch on a specific R console and it this specific console is still running after a laps of time (='Time')
# we forced the stop of the function and call to a no convergence for the iteration and then go to the following iteration
# This method adds more time of calculation but it is the only method that we found to handle this infinite time running off the 'gee' function
# Also, this method create separate R file to lunch the separate iteration to each specific R consoles

# Input
# i               : Iteration number 'i'
# Time            : Specific time chosen to decide to stop the gee function if it is running for an infinite time 
# n               : Scenario number used
# Base_file       : Path of the file where the R workspace is
# Workspace_name  : Name of the R workspace
# Data_file       : Path of the file where the data are stock
# Resu_file       : Path of the file where the results are stock

# Output
# Will stock directly the result in the excel file as the function "fun_para_analyse_%method" does



fun_cor_gauss <- function(i,Time,n,Base_file,Workspace_name,Data_file,Resu_file,Cor.FG = TRUE){
  
  # 0. File name
  Chemin = sprintf(paste(Base_file,"/correction/correction_Gauss/correction_gauss_S%d",sep = ""),n)
  dir.create(Chemin)
  setwd(Chemin)
  
  # 0.1 Name of the principal file
  FileName = paste("correction_gauss_S",n,"_itt_%d.R",sep = "")
  
  
  # 1. Creating R file with the instructions for the second session
  
  
  Stock_file = sprintf(paste(Base_file,"/correction/correction_Gauss/correction_gauss_S%d",sep = ""),n)
  
  writeLines('rm(list = ls())',
             con = sprintf(FileName,i))
  
  
  write(sprintf("Chemin  = '%s'",Stock_file),
        file  = sprintf(FileName, i),
        append = TRUE)  
  
  write("setwd(Chemin)",
        file  = sprintf(FileName, i),
        append = TRUE)
  
  # 1.1 Export ID job of the second file to kill it if it is necessary (infinite lap of time)
  Exp_id = paste("writeLines(as.character(Sys.getpid()), con = sprintf('%s/pid_%d.txt', getwd(),",i,"))",sep = "")
  write(Exp_id,
        file = sprintf(FileName, i),
        append = TRUE)
  
  # 1.2 Write the second file 
  ## Care to change the library repository, name file to change for the library command.
  ## If it is local, no need to input the "lib.loc" in library function
  
  write(paste("load('",Base_file,"/",Workspace_name,"')",sep = ""), file = sprintf(FileName, i), append = TRUE)
  
  write("library(cli, lib.loc='/home/jpereira/z.libs_Pierre/common')",file = sprintf(FileName, i), append = TRUE)
  write("library(Matrix, lib.loc='/home/jpereira/z.libs_Pierre/common')",file = sprintf(FileName, i), append = TRUE)
  write("library(MASS, lib.loc='/home/jpereira/z.libs_Pierre/common')",file = sprintf(FileName, i), append = TRUE)
  write("library(lme4, lib.loc='/home/jpereira/z.libs_Pierre/common')",file = sprintf(FileName, i), append = TRUE)
  write("library(backports, lib.loc='/home/jpereira/z.libs_Pierre/common')",file = sprintf(FileName, i), append = TRUE)
  write("library(dplyr, lib.loc='/home/jpereira/z.libs_Pierre/common')",file = sprintf(FileName, i), append = TRUE)
  write("library(parallel, lib.loc='/home/jpereira/z.libs_Pierre/common')",file = sprintf(FileName, i), append = TRUE)
  write("library(iterators, lib.loc='/home/jpereira/z.libs_Pierre/common')",file = sprintf(FileName, i), append = TRUE)
  write("library(rngtools, lib.loc='/home/jpereira/z.libs_Pierre/common')",file = sprintf(FileName, i), append = TRUE)
  write("library(foreach, lib.loc='/home/jpereira/z.libs_Pierre/common')",file = sprintf(FileName, i), append = TRUE)
  write("library(doRNG, lib.loc='/home/jpereira/z.libs_Pierre/common')",file = sprintf(FileName, i), append = TRUE)
  write("library(doParallel, lib.loc='/home/jpereira/z.libs_Pierre/common')",file = sprintf(FileName, i), append = TRUE)
  write("library(dplyr, lib.loc='/home/jpereira/z.libs_Pierre/common')",file = sprintf(FileName, i), append = TRUE)
  write("library(gee, lib.loc='/home/jpereira/z.libs_Pierre/common')",file = sprintf(FileName, i), append = TRUE)
  write("library(geepack, lib.loc='/home/jpereira/z.libs_Pierre/common')",file = sprintf(FileName, i), append = TRUE)
  write("library(spind, lib.loc='/home/jpereira/z.libs_Pierre/common')",file = sprintf(FileName, i), append = TRUE)
  write("library(doBy, lib.loc='/home/jpereira/z.libs_Pierre/common')",file = sprintf(FileName, i), append = TRUE)
  write("library(arm, lib.loc='/home/jpereira/z.libs_Pierre/common')",file = sprintf(FileName, i), append = TRUE)
  write("library(here, lib.loc='/home/jpereira/z.libs_Pierre/common')",file = sprintf(FileName, i), append = TRUE)
  write("library(geesmv, lib.loc='/home/jpereira/z.libs_Pierre/common')",file = sprintf(FileName, i), append = TRUE)
  write("library(matrixcalc, lib.loc='/home/jpereira/z.libs_Pierre/common')",file = sprintf(FileName, i), append = TRUE)  
  
  L1 = paste("n =",n)
  write(L1,
        file = sprintf(FileName,i),
        append = TRUE)
  
  write(paste("Data_file = '",Data_file,"'",sep = ""), file = sprintf(FileName, i), append = TRUE)
  write(paste("Resu_file = '",Resu_file,"'",sep = ""), file = sprintf(FileName, i), append = TRUE)
  write("G_res_file = paste(Resu_file,'/Data_Output_gaussid',sep = '')", file = sprintf(FileName, i), append = TRUE)
  
  write("matcor_type = 'exchangeable'",
        file = sprintf(FileName, i),
        append = TRUE)
  write(sprintf("Scenario_use = Scenario_%d",n),
        file = sprintf(FileName, i),
        append = TRUE)
  
  
  if(Cor.FG == TRUE){write(sprintf("fun_para_analyse_gauss(itt_para  = %d,n,matcor_type = matcor_type,Scenario = Scenario_use,Cor.FG = TRUE,Data_file,G_res_file)",i),
                           file = sprintf(FileName, i),
                           append = TRUE)}
  if(Cor.FG == FALSE){write(sprintf("fun_para_analyse_gauss(itt_para  = %d,n,matcor_type = matcor_type,Scenario = Scenario_use,Cor.FG = FALSE,Data_file,G_res_file)",i),
                            file = sprintf(FileName, i),
                            append = TRUE)}
  
  write("DONE = TRUE",file = sprintf(FileName,i),append = TRUE)
  
  # 1.3 Give the information if the iteration did converge or not
  tace_2 = paste("save(DONE, file = sprintf('res_%d.Rdata',",i,"))",sep = "")
  write(tace_2, file = sprintf(FileName, i), append = TRUE)
  
  # 2. Execute the batch file to lunch the Rscript in the T lap time in background (correction_S33_itt_x.R)
  cmd_batch = paste("R CMD BATCH --no-restore correction_gauss_S",n,"_itt_%d.R",sep = "")
  
  Start.T = Sys.time()
  system(sprintf(cmd_batch, i), wait = TRUE,timeout = Time)
  End.T = Sys.time()
  
  Diff.time = difftime(End.T,Start.T,units = "sec")
  if(Diff.time >= Time ){
    system(sprintf("tskill %d", scan(sprintf("pid_%d.txt", i), integer(),quiet = T)))
  }
}

fun_cor_bin <- function(i,Time,n,Base_file,Workspace_name,Data_file,Resu_file,Cor.FG=TRUE){
  
  # 0. File name
  Chemin = sprintf(paste(Base_file,"/correction/correction_Bin/correction_bin_S%d",sep = ""),n)
  dir.create(Chemin)
  setwd(Chemin)
  
  # 0.1 Name of the principal file
  FileName = paste("correction_bin_S",n,"_itt_%d.R",sep = "")
  
  
  # 1. Creating R file with the instructions for the second session
  
  
  Stock_file = sprintf(paste(Base_file,"/correction/correction_Bin/correction_bin_S%d",sep = ""),n)
  
  writeLines('rm(list = ls())',
             con = sprintf(FileName,i))
  
  
  write(sprintf("Chemin  = '%s'",Stock_file),
        file  = sprintf(FileName, i),
        append = TRUE)  
  
  write("setwd(Chemin)",
        file  = sprintf(FileName, i),
        append = TRUE)
  
  # 1.1 Export ID job of the second file to kill it if it is necessary (infinite lap of time)
  Exp_id = paste("writeLines(as.character(Sys.getpid()), con = sprintf('%s/pid_%d.txt', getwd(),",i,"))",sep = "")
  write(Exp_id,
        file = sprintf(FileName, i),
        append = TRUE)
  
  # 1.2 Write the second file 
  ## Care to change the library repository, name file to change for the library command.
  ## If it is local, no need to input the "lib.loc" in library function
  
  write(paste("load('",Base_file,"/",Workspace_name,"')",sep = ""), file = sprintf(FileName, i), append = TRUE)
  
  write("library(cli, lib.loc='/home/jpereira/z.libs_Pierre/common')",file = sprintf(FileName, i), append = TRUE)
  write("library(Matrix, lib.loc='/home/jpereira/z.libs_Pierre/common')",file = sprintf(FileName, i), append = TRUE)
  write("library(MASS, lib.loc='/home/jpereira/z.libs_Pierre/common')",file = sprintf(FileName, i), append = TRUE)
  write("library(lme4, lib.loc='/home/jpereira/z.libs_Pierre/common')",file = sprintf(FileName, i), append = TRUE)
  write("library(backports, lib.loc='/home/jpereira/z.libs_Pierre/common')",file = sprintf(FileName, i), append = TRUE)
  write("library(dplyr, lib.loc='/home/jpereira/z.libs_Pierre/common')",file = sprintf(FileName, i), append = TRUE)
  write("library(parallel, lib.loc='/home/jpereira/z.libs_Pierre/common')",file = sprintf(FileName, i), append = TRUE)
  write("library(iterators, lib.loc='/home/jpereira/z.libs_Pierre/common')",file = sprintf(FileName, i), append = TRUE)
  write("library(rngtools, lib.loc='/home/jpereira/z.libs_Pierre/common')",file = sprintf(FileName, i), append = TRUE)
  write("library(foreach, lib.loc='/home/jpereira/z.libs_Pierre/common')",file = sprintf(FileName, i), append = TRUE)
  write("library(doRNG, lib.loc='/home/jpereira/z.libs_Pierre/common')",file = sprintf(FileName, i), append = TRUE)
  write("library(doParallel, lib.loc='/home/jpereira/z.libs_Pierre/common')",file = sprintf(FileName, i), append = TRUE)
  write("library(dplyr, lib.loc='/home/jpereira/z.libs_Pierre/common')",file = sprintf(FileName, i), append = TRUE)
  write("library(gee, lib.loc='/home/jpereira/z.libs_Pierre/common')",file = sprintf(FileName, i), append = TRUE)
  write("library(geepack, lib.loc='/home/jpereira/z.libs_Pierre/common')",file = sprintf(FileName, i), append = TRUE)
  write("library(spind, lib.loc='/home/jpereira/z.libs_Pierre/common')",file = sprintf(FileName, i), append = TRUE)
  write("library(doBy, lib.loc='/home/jpereira/z.libs_Pierre/common')",file = sprintf(FileName, i), append = TRUE)
  write("library(arm, lib.loc='/home/jpereira/z.libs_Pierre/common')",file = sprintf(FileName, i), append = TRUE)
  write("library(here, lib.loc='/home/jpereira/z.libs_Pierre/common')",file = sprintf(FileName, i), append = TRUE)
  write("library(geesmv, lib.loc='/home/jpereira/z.libs_Pierre/common')",file = sprintf(FileName, i), append = TRUE)
  write("library(matrixcalc, lib.loc='/home/jpereira/z.libs_Pierre/common')",file = sprintf(FileName, i), append = TRUE) 
  
  L1 = paste("n =",n)
  write(L1,
        file = sprintf(FileName,i),
        append = TRUE)
  
  write(paste("Data_file = '",Data_file,"'",sep = ""), file = sprintf(FileName, i), append = TRUE)
  write(paste("Resu_file = '",Resu_file,"'",sep = ""), file = sprintf(FileName, i), append = TRUE)
  write("B_res_file = paste(Resu_file,'/Data_Output_binlogit',sep = '')", file = sprintf(FileName, i), append = TRUE)
  
  write("matcor_type = 'exchangeable'",
        file = sprintf(FileName, i),
        append = TRUE)
  write(sprintf("Scenario_use = Scenario_%d",n),
        file = sprintf(FileName, i),
        append = TRUE)
  
  
  if(Cor.FG == TRUE){write(sprintf("fun_para_analyse_bin(itt_para  = %d,n,matcor_type = matcor_type,Scenario = Scenario_use,Cor.FG = TRUE,Data_file,B_res_file)",i),
                           file = sprintf(FileName, i),
                           append = TRUE)}
  if(Cor.FG == FALSE){write(sprintf("fun_para_analyse_bin(itt_para  = %d,n,matcor_type = matcor_type,Scenario = Scenario_use,Cor.FG = FALSE,Data_file,B_res_file)",i),
                            file = sprintf(FileName, i),
                            append = TRUE)}
  
  write("DONE = TRUE",file = sprintf(FileName,i),append = TRUE)
  
  # 1.3 Give the information if the iteration did converge or not
  tace_2 = paste("save(DONE, file = sprintf('res_%d.Rdata',",i,"))",sep = "")
  write(tace_2, file = sprintf(FileName, i), append = TRUE)
  
  # 2. Execute the batch file to lunch the Rscript in the T lap time in background (correction_S33_itt_x.R)
  cmd_batch = paste("R CMD BATCH --no-restore correction_bin_S",n,"_itt_%d.R",sep = "")
  
  Start.T = Sys.time()
  system(sprintf(cmd_batch, i), wait = TRUE,timeout = Time)
  End.T = Sys.time()
  
  Diff.time = difftime(End.T,Start.T,units = "sec")
  if(Diff.time >= Time ){
    system(sprintf("tskill %d", scan(sprintf("pid_%d.txt", i), integer(),quiet = T)))
  }
}

fun_cor_poiss <- function(i,Time,n,Base_file,Workspace_name,Data_file,Resu_file,Cor.FG = TRUE){
  
  # 0. File name
  Chemin = sprintf(paste(Base_file,"/correction/correction_Poiss/correction_poiss_S%d",sep = ""),n)
  dir.create(Chemin)
  setwd(Chemin)
  
  # 0.1 Name of the principal file
  FileName = paste("correction_poiss_S",n,"_itt_%d.R",sep = "")
  
  
  # 1. Creating R file with the instructions for the second session
  
  Stock_file = sprintf(paste(Base_file,"/correction/correction_Poiss/correction_poiss_S%d",sep = ""),n)
  
  writeLines('rm(list = ls())',
             con = sprintf(FileName,i))
  
  
  write(sprintf("Chemin  = '%s'",Stock_file),
        file  = sprintf(FileName, i),
        append = TRUE)  
  
  write("setwd(Chemin)",
        file  = sprintf(FileName, i),
        append = TRUE)
  
  # 1.1 Export ID job of the second file to kill it if it is necessary (infinite lap of time)
  Exp_id = paste("writeLines(as.character(Sys.getpid()), con = sprintf('%s/pid_%d.txt', getwd(),",i,"))",sep = "")
  write(Exp_id,
        file = sprintf(FileName, i),
        append = TRUE)
  
  # 1.2 Write the second file 
  ## Care to change the library repository, name file to change for the library command.
  ## If it is local, no need to input the "lib.loc" in library function
  
  write(paste("load('",Base_file,"/",Workspace_name,"')",sep = ""), file = sprintf(FileName, i), append = TRUE)
  
  write("library(cli, lib.loc='/home/jpereira/z.libs_Pierre/common')",file = sprintf(FileName, i), append = TRUE)
  write("library(Matrix, lib.loc='/home/jpereira/z.libs_Pierre/common')",file = sprintf(FileName, i), append = TRUE)
  write("library(MASS, lib.loc='/home/jpereira/z.libs_Pierre/common')",file = sprintf(FileName, i), append = TRUE)
  write("library(lme4, lib.loc='/home/jpereira/z.libs_Pierre/common')",file = sprintf(FileName, i), append = TRUE)
  write("library(backports, lib.loc='/home/jpereira/z.libs_Pierre/common')",file = sprintf(FileName, i), append = TRUE)
  write("library(dplyr, lib.loc='/home/jpereira/z.libs_Pierre/common')",file = sprintf(FileName, i), append = TRUE)
  write("library(parallel, lib.loc='/home/jpereira/z.libs_Pierre/common')",file = sprintf(FileName, i), append = TRUE)
  write("library(iterators, lib.loc='/home/jpereira/z.libs_Pierre/common')",file = sprintf(FileName, i), append = TRUE)
  write("library(rngtools, lib.loc='/home/jpereira/z.libs_Pierre/common')",file = sprintf(FileName, i), append = TRUE)
  write("library(foreach, lib.loc='/home/jpereira/z.libs_Pierre/common')",file = sprintf(FileName, i), append = TRUE)
  write("library(doRNG, lib.loc='/home/jpereira/z.libs_Pierre/common')",file = sprintf(FileName, i), append = TRUE)
  write("library(doParallel, lib.loc='/home/jpereira/z.libs_Pierre/common')",file = sprintf(FileName, i), append = TRUE)
  write("library(dplyr, lib.loc='/home/jpereira/z.libs_Pierre/common')",file = sprintf(FileName, i), append = TRUE)
  write("library(gee, lib.loc='/home/jpereira/z.libs_Pierre/common')",file = sprintf(FileName, i), append = TRUE)
  write("library(geepack, lib.loc='/home/jpereira/z.libs_Pierre/common')",file = sprintf(FileName, i), append = TRUE)
  write("library(spind, lib.loc='/home/jpereira/z.libs_Pierre/common')",file = sprintf(FileName, i), append = TRUE)
  write("library(doBy, lib.loc='/home/jpereira/z.libs_Pierre/common')",file = sprintf(FileName, i), append = TRUE)
  write("library(arm, lib.loc='/home/jpereira/z.libs_Pierre/common')",file = sprintf(FileName, i), append = TRUE)
  write("library(here, lib.loc='/home/jpereira/z.libs_Pierre/common')",file = sprintf(FileName, i), append = TRUE)
  write("library(geesmv, lib.loc='/home/jpereira/z.libs_Pierre/common')",file = sprintf(FileName, i), append = TRUE)
  write("library(matrixcalc, lib.loc='/home/jpereira/z.libs_Pierre/common')",file = sprintf(FileName, i), append = TRUE) 
  
  L1 = paste("n =",n)
  write(L1,
        file = sprintf(FileName,i),
        append = TRUE)
  
  write(paste("Data_file = '",Data_file,"'",sep = ""), file = sprintf(FileName, i), append = TRUE)
  write(paste("Resu_file = '",Resu_file,"'",sep = ""), file = sprintf(FileName, i), append = TRUE)
  write("P_res_file = paste(Resu_file,'/Data_Output_poisslog',sep = '')", file = sprintf(FileName, i), append = TRUE)
  
  
  write("matcor_type = 'exchangeable'",
        file = sprintf(FileName, i),
        append = TRUE)
  write(sprintf("Scenario_use = Scenario_%d",n),
        file = sprintf(FileName, i),
        append = TRUE)
  
  if(Cor.FG == TRUE){write(sprintf("fun_para_analyse_poiss(itt_para  = %d,n,matcor_type = matcor_type,Scenario = Scenario_use,Cor.FG = TRUE,Data_file,P_res_file)",i),
                           file = sprintf(FileName, i),
                           append = TRUE)}
  if(Cor.FG == FALSE){write(sprintf("fun_para_analyse_poiss(itt_para  = %d,n,matcor_type = matcor_type,Scenario = Scenario_use,Cor.FG = FALSE,Data_file,P_res_file)",i),
                            file = sprintf(FileName, i),
                            append = TRUE)}
  
  
  write("DONE = TRUE",file = sprintf(FileName,i),append = TRUE)
  
  # 1.3 Give the information if the iteration did converge or not
  tace_2 = paste("save(DONE, file = sprintf('res_%d.Rdata',",i,"))",sep = "")
  write(tace_2, file = sprintf(FileName, i), append = TRUE)
  
  # 2. Execute the batch file to lunch the Rscript in the T lap time in background (correction_S33_itt_x.R)
  cmd_batch = paste("R CMD BATCH --no-restore correction_poiss_S",n,"_itt_%d.R",sep = "")
  
  
  Start.T = Sys.time()
  system(sprintf(cmd_batch, i), wait = TRUE,timeout = Time)
  End.T = Sys.time()
  
  Diff.time = difftime(End.T,Start.T,units = "sec")
  if(Diff.time >= Time ){
    system(sprintf("tskill %d", scan(sprintf("pid_%d.txt", i), integer(),quiet = T)))
  }
}




###########################################################################################################################################


### GLMM Functions ----

# 12) Delta method functions for glmer ----

#fun_se_delta_meth_logit_glmer

# Input
# fit           : fit obtained with the 'glmer' function of the 'lme4' package
# nb_cov        : Number of individual level covariates
# nb_cov_clus   : Number of cluster level covariates
# data          : data set 

# Output
# se : Standard error calculated for the risk difference with the delta method after performing g-computation with binomial distribution and logit link function

fun_se_delta_meth_logit_glmer <- function(fit,nb_cov,nb_cov_clus,data){
  
  d = c()
  
  sum_fit = summary(fit)
  
  x1 <- x0 <- data
  x1$Arm <- 1
  x0$Arm <- 0
  
  x1_mod = cbind(Intercept = rep(1,nrow(x1)),x1[,-c(1,3)])
  x0_mod = cbind(Intercept = rep(1,nrow(x0)),x0[,-c(1,3)])
  
  Trt   <- predict(fit,newdata = x1)
  NoTrt <- predict(fit,newdata = x0)
  
  fixed_beta = as.vector(sum_fit$coefficients[,1])
  nb_var = length(fixed_beta)
  
  for (j in 1:nb_var) {
    d_j = mean(deriv_logit(Trt) * x1_mod[,j] -  deriv_logit(NoTrt) * x0_mod[,j])
    d=c(d,d_j)
  }
  
  se = sqrt(d %*% vcov(fit) %*% d)
  
  return(se)
}

#fun_se_delta_meth_log_glmer

# Input
# fit           : fit obtained with the 'glmer' function of the 'lme4' package
# nb_cov        : Number of individual level covariates
# nb_cov_clus   : Number of cluster level covariates
# data          : data set 

# Output
# se : Standard error calculated for the risk difference with the delta method after performing g-computation with Poisson distribution and log link function


fun_se_delta_meth_log_glmer <- function(fit,nb_cov,nb_cov_clus,data){
  
  d = c()
  sum_fit = summary(fit)
  
  
  x1 <- x0 <- data
  x1$Arm <- 1
  x0$Arm <- 0
  
  x1_mod = cbind(Intercept = rep(1,nrow(x1)),x1[,-c(1,3)])
  x0_mod = cbind(Intercept = rep(1,nrow(x0)),x0[,-c(1,3)])
  
  Trt   <- predict(fit,newdata = x1)
  NoTrt <- predict(fit,newdata = x0)
  
  fixed_beta = as.vector(sum_fit$coefficients[,1])
  nb_var = length(fixed_beta)
  
  for (j in 1:nb_var) {
    d_j = mean(exp(Trt) * x1_mod[,j] -  exp(NoTrt) * x0_mod[,j])
    d=c(d,d_j)
  }
  
  se = sqrt(d %*% vcov(fit) %*% d)
  
  return(se)
}

# 13) Auto formula for glmer functions ----

# Automating the formula function by the number of covariate in our GLMM model.
# Input
# nb_cov : number of individual level covariate of the data set.
# nb_cov_clus : number of cluster level covariate of the data set.
#
# Output
# Formula : An object at the formula type, getting the formula to use for glmer function

fun_formula_for_glmer <- function(nb_cov,nb_cov_clus){
  
  # FIXED EFFECT ----
  
  Covariable = ''
  
  #### Case nb_cov = 0 ----
  if(nb_cov==0){
    if(nb_cov_clus==0){
      Covariable = ''
    }
    if(nb_cov_clus==1){
      Covariable = paste(Covariable,'Covariate_clus',sep = '+')
    }
    if(nb_cov_clus > 1){
      for (i in 1:nb_cov_clus) {
        a = paste('X',i,sep = "")
        Covariable = paste(Covariable,a,sep = '+')
      }
    }
  }
  
  #### Case nb_cov = 1 ----
  if(nb_cov==1){
    Covariable = '+ Covariate'
    if(nb_cov_clus==1){
      Covariable = paste(Covariable,'Covariate_clus',sep = '+')
    }
    if(nb_cov_clus > 1){
      for (i in 1:(nb_cov_clus)) {
        a = paste('X',i,sep = "")
        Covariable = paste(Covariable,a,sep = '+')
      }
    }
  }
  
  #### Case nb_cov > 1 ----
  if(nb_cov > 1){
    for (i in 1:nb_cov) {
      a = paste('X',i,sep = "")
      Covariable = paste(Covariable,a,sep = '+')
    }
    if(nb_cov_clus == 0){
      Covariable = Covariable
    }
    if(nb_cov_clus == 1){
      Covariable = paste(Covariable,'Covariate_clus',sep = '+')
    }
    if(nb_cov_clus > 1){
      if(nb_cov_clus <= nb_cov){
        for (i in 1:(nb_cov_clus)) {
          a = paste('X',i,'.1',sep = "")
          Covariable = paste(Covariable,a,sep = '+')
        }
      }
      if(nb_cov_clus > nb_cov){
        for (i in 1:nb_cov) {
          a = paste('X',i,'.1',sep = "")
          Covariable = paste(Covariable,a,sep = '+')
        }
        for (i in (nb_cov+1):nb_cov_clus) {
          a = paste('X',i,sep = "")
          Covariable = paste(Covariable,a,sep = '+')
        }
      }
    }
  }
  
  
  Covariable = as.factor(Covariable)
  Random_effect = "+ (1|cluster)"
  Formula = as.formula(paste("Outcome ~ Arm ",Covariable,Random_effect))
  
  return(Formula)
}

# 14) Function simulation for GLMM ----

# Input : 
# data_itt    : data set generated as a Data.frame
# itt_para    : which itteration of the paralellism is analyzed (it's too know which data we use on the 1000 simulated)
# Scenario    : Scenario corresponding to the the data set analyzed

# Output :
# Results of an iteration ('itt_para') for one method by a table with the results : (Risk difference estimated, CI at 95%, Coverage_rage,Bias ....)

fun_sim_gauss_glmer <- function(data_itt,itt_para,Scenario = list(pI,pC,k,m,icc,nb_cov,nb_cov_clus,p_cov,p_cov_clus,OR_cov,OR_cov_clus)){
  
  # 
  
  ## Scenario parameters ----
  pI = as.numeric(Scenario[1])
  pC = as.numeric(Scenario[2])
  k = as.numeric(Scenario[3])
  m = as.numeric(Scenario[4])
  icc = as.numeric(Scenario[5])
  nb_cov = as.numeric(Scenario[6])
  nb_cov_clus = as.numeric(Scenario[7])
  p_cov = as.vector(unlist(Scenario[8]))
  p_cov_clus = as.vector(unlist(Scenario[9]))
  OR_cov = as.vector(unlist(Scenario[10]))
  OR_cov_clus = as.vector(unlist(Scenario[11]))
  
  RD_th = pI-pC
  nb_cov_tot = nb_cov+nb_cov_clus
  
  ## Vectors empty needed ----
  RD_itt_gaussid = c()
  SE_itt_tot_gaussid = c()
  LL_CI_RD_itt_gaussid = c()
  UL_CI_RD_itt_gaussid = c()
  
  ## Formula for glmer ----
  
  Formula_glmer = fun_formula_for_glmer(nb_cov,nb_cov_clus)
  
  ## GLMMM  ----
  
  ## Gaussian Identity  ----
  
  cap_op  <- capture.output(suppressMessages(res <- fit_gaussid <- try(lmer(formula = Formula_glmer,
                                                                            data = data_itt,
                                                                            REML = TRUE),silent = TRUE)))
  sum_fit = summary(fit_gaussid)
  
  if(inherits(res, "try-error")==FALSE){
    
    
    summary_fit_gaussid = summary(fit_gaussid)
    
    RD_gaussid = as.numeric(summary_fit_gaussid$coefficients[2,1])
    RD_itt_gaussid = c(RD_itt_gaussid,RD_gaussid)
    
    #Degree of freedom Corrected
    
    df_BW = 2*k - (nb_cov_tot + 1 + 1) # Number of cluster - number of fixed effect (nb_cov_tot + Arm + intercept)
    se_tot_gaussid = summary_fit_gaussid$coefficient[2,2]
    SE_itt_tot_gaussid = c(SE_itt_tot_gaussid,se_tot_gaussid)
    
    
    # CI RD
    
    LL_CI_RD_gaussid = RD_gaussid - qt(0.975,df_BW)*se_tot_gaussid
    LL_CI_RD_itt_gaussid = c( LL_CI_RD_itt_gaussid, LL_CI_RD_gaussid)
    
    UL_CI_RD_gaussid = RD_gaussid + qt(0.975,df_BW)*se_tot_gaussid
    UL_CI_RD_itt_gaussid =c(UL_CI_RD_itt_gaussid,UL_CI_RD_gaussid)
    
    ## Results Gauss id ----
    # RD_mean_gaussid = mean(RD_itt_gaussid)
    
    itt_without_error_gaussid = length(RD_itt_gaussid)
    Coverage_rate_gaussid = length(which(RD_th > LL_CI_RD_itt_gaussid & RD_th < UL_CI_RD_itt_gaussid))/length(RD_itt_gaussid)*100
    Absolute_Bias_itt_gaussid = abs(RD_th - RD_itt_gaussid)
    Bias_itt_gaussid = RD_itt_gaussid - RD_th
    Relative_Bias_itt_gaussid = RD_itt_gaussid/RD_th
    res_gaussid = data_frame(Risk_difference_Estimate_gaussid = RD_itt_gaussid,
                             Number_of_itteration_gaussid = itt_without_error_gaussid,
                             LL_95 = LL_CI_RD_itt_gaussid,
                             UL_95 = UL_CI_RD_itt_gaussid,
                             SE_itt = SE_itt_tot_gaussid,
                             Coverage_rate_gaussid = Coverage_rate_gaussid,
                             Abs_Bias_iteration_gaussid = Absolute_Bias_itt_gaussid,
                             Bias_iteration_gaussid = Bias_itt_gaussid,
                             Relative_bias_iteration_gaussid = Relative_Bias_itt_gaussid,
                             itt_para = itt_para,
                             type_error = sum_fit$optinfo$conv$opt)
    
    res_gaussid[unique(which(is.na(res_gaussid) | abs(res_gaussid$Risk_difference_Estimate_gaussid)>=1 | res_gaussid$type_error != 0,arr.ind = TRUE)[,1]),1:9] = NA
  }else{
    res_gaussid = data_frame(Risk_difference_Estimate_gaussid = NA,
                             Number_of_itteration_gaussid = NA,
                             LL_95 = NA,
                             UL_95 = NA,
                             SE_itt = NA,
                             Coverage_rate_gaussid = NA,
                             Abs_Bias_iteration_gaussid = NA,
                             Bias_iteration_gaussid = NA,
                             Relative_bias_iteration_gaussid = NA,
                             itt_para = itt_para,
                             type_error = "error")
  }
  
  return(res_gaussid)
  
}

# Input : 
# data_itt    : data set generated as a Data.frame
# itt_para    : which itteration of the paralellism is analyzed (it's too know which data we use on the 1000 simulated)
# Scenario    : Scenario corresponding to the the data set analyzed

# Output :
# Results of an iteration for one method by a table with the results : (Risk difference estimated, CI at 95%, Coverage_rage,Bias ....)


fun_sim_binom_glmer <- function(data_itt,itt_para,Scenario = list(pI,pC,k,m,icc,nb_cov,nb_cov_clus,p_cov,p_cov_clus,OR_cov,OR_cov_clus)){
  
  
  ## Scenario parameters ----
  pI = as.numeric(Scenario[1])
  pC = as.numeric(Scenario[2])
  k = as.numeric(Scenario[3])
  m = as.numeric(Scenario[4])
  icc = as.numeric(Scenario[5])
  nb_cov = as.numeric(Scenario[6])
  nb_cov_clus = as.numeric(Scenario[7])
  p_cov = as.vector(unlist(Scenario[8]))
  p_cov_clus = as.vector(unlist(Scenario[9]))
  OR_cov = as.vector(unlist(Scenario[10]))
  OR_cov_clus = as.vector(unlist(Scenario[11]))
  
  RD_th = pI-pC
  nb_cov_tot = nb_cov+nb_cov_clus
  
  
  # initial_value_gee = c(1,0,rep(0.5,(nb_cov_tot)))
  ## Vectors empty needed ----
  
  RD_itt_binlogit = c()
  SE_itt_tot_binlogit = c()
  LL_CI_RD_itt_binlogit = c()
  UL_CI_RD_itt_binlogit = c()
  
  ## Formula for glmer ----
  
  Formula_glmer = fun_formula_for_glmer(nb_cov,nb_cov_clus)
  
  
  
  ## GLMM  ----
  ## Binomial Logit  ----
  
  cap_op2  <- capture.output(suppressMessages(res <- fit_binlogit <- try(glmer(formula = Formula_glmer,
                                                                               data = data_itt,
                                                                               family = binomial()),silent = TRUE)))
  
  sum_fit = summary(fit_binlogit)
  
  if(inherits(res, "try-error")==FALSE ){
    
    ### ALL include ----
    
    # Margin mean
    Trt_data <- NoTrt_data <- data_itt
    Trt_data$Arm <- 1 
    NoTrt_data$Arm <- 0
    
    Pred_Trt = predict(object = fit_binlogit,newdata = Trt_data)
    Pred_NoTrt = predict(object = fit_binlogit,newdata = NoTrt_data)
    
    P_Trt_est_binlogit = mean(invlogit(Pred_Trt))
    P_NoTrt_est_binlogit = mean(invlogit(Pred_NoTrt))
    
    RD_binlogit = P_Trt_est_binlogit - P_NoTrt_est_binlogit
    RD_itt_binlogit = c(RD_itt_binlogit,RD_binlogit)
    
    #Degree of freedom Corrected
    
    df_BW = 2*k - (nb_cov_tot + 1 + 1) # Number of cluster - number of fixed effect (nb_cov_tot + Arm + intercept)
    
    se_tot_binlogit = fun_se_delta_meth_logit_glmer(fit_binlogit,nb_cov,nb_cov_clus,data_itt)
    se_tot_binlogit = as.numeric(se_tot_binlogit)
    SE_itt_tot_binlogit = c(SE_itt_tot_binlogit,se_tot_binlogit)
    
    # CI RD 
    LL_CI_RD_binlogit = RD_binlogit - qt(0.975,(df_BW))*se_tot_binlogit
    LL_CI_RD_itt_binlogit = c( LL_CI_RD_itt_binlogit, LL_CI_RD_binlogit)
    
    UL_CI_RD_binlogit = RD_binlogit + qt(0.975,(df_BW))*se_tot_binlogit
    UL_CI_RD_itt_binlogit =c(UL_CI_RD_itt_binlogit,UL_CI_RD_binlogit)
    
    ## Results Bin logit ----
    # RD_mean_binlogit = mean(RD_itt_binlogit)
    itt_without_error_binlogit = length(RD_itt_binlogit)
    Coverage_rate_binlogit = length(which(RD_th > LL_CI_RD_itt_binlogit & RD_th < UL_CI_RD_itt_binlogit))/length(RD_itt_binlogit)*100
    Absolute_Bias_itt_binlogit = abs(RD_th - RD_itt_binlogit)
    Bias_itt_binlogit = RD_itt_binlogit - RD_th 
    Relative_Bias_itt_binlogit = RD_itt_binlogit/RD_th
    res_binlogit = data_frame(Risk_difference_Estimate_binlogit = RD_itt_binlogit,
                              Number_of_itteration_binlogit = itt_without_error_binlogit,
                              LL_95 = LL_CI_RD_itt_binlogit,
                              UL_95 = UL_CI_RD_itt_binlogit,
                              SE_itt = SE_itt_tot_binlogit,
                              Coverage_rate_binlogit = Coverage_rate_binlogit,
                              Abs_Bias_iteration_binlogit = Absolute_Bias_itt_binlogit,
                              Bias_iteration_binlogit = Bias_itt_binlogit,
                              Relative_bias_iteration_binlogit = Relative_Bias_itt_binlogit,
                              itt_para = itt_para,
                              type_error = sum_fit$optinfo$conv$opt)
    
    res_binlogit[unique(which(is.na(res_binlogit) | abs(res_binlogit$Risk_difference_Estimate_binlogit) >= 1 | res_binlogit$type_error != 0 ,arr.ind = TRUE)[,1]),1:9] = NA
  }else{
    res_binlogit = data_frame(Risk_difference_Estimate_binlogit = NA,
                              Number_of_itteration_binlogit = NA,
                              LL_95 = NA,
                              UL_95 = NA,
                              SE_itt = NA,
                              Coverage_rate_binlogit = NA,
                              Abs_Bias_iteration_binlogit = NA,
                              Bias_iteration_binlogit = NA,
                              Relative_bias_iteration_binlogit = NA,
                              itt_para = itt_para,
                              type_error = "error")
  }
  
  
  
  
  
  return(res_binlogit)
  
}

# Input : 
# data_itt    : data set generated as a Data.frame
# itt_para    : which itteration of the paralellism is analyzed (it's too know which data we use on the 1000 simulated)
# Scenario    : Scenario corresponding to the the data set analyzed

# Output :
# Results of an iteration for one method by a table with the results : (Risk difference estimated, CI at 95%, Coverage_rage,Bias ....)


fun_sim_poiss_glmer <- function(data_itt,itt_para,matcor_type = "exchangeable",Scenario = list(pI,pC,k,m,icc,nb_cov,nb_cov_clus,p_cov,p_cov_clus,OR_cov,OR_cov_clus),Cor.FG=TRUE){
  
  ## Scenario parameters ----
  pI = as.numeric(Scenario[1])
  pC = as.numeric(Scenario[2])
  k = as.numeric(Scenario[3])
  m = as.numeric(Scenario[4])
  icc = as.numeric(Scenario[5])
  nb_cov = as.numeric(Scenario[6])
  nb_cov_clus = as.numeric(Scenario[7])
  p_cov = as.vector(unlist(Scenario[8]))
  p_cov_clus = as.vector(unlist(Scenario[9]))
  OR_cov = as.vector(unlist(Scenario[10]))
  OR_cov_clus = as.vector(unlist(Scenario[11]))
  
  RD_th = pI-pC
  nb_cov_tot = nb_cov+nb_cov_clus
  
  initial_value_gee = c(1,0,rep(0.5,(nb_cov_tot)))
  ## Vectors empty needed ----
  
  RD_itt_poisslog = c()
  SE_itt_tot_poisslog = c()
  LL_CI_RD_itt_poisslog = c()
  UL_CI_RD_itt_poisslog = c()
  
  ## Formula for glmer ----
  
  Formula_glmer = fun_formula_for_glmer(nb_cov,nb_cov_clus)
  
  ## GLMM  ----
  
  ## Poisson Log  ----
  
  
  cap_op2  <- capture.output(suppressMessages(res <- fit_poisslog <- try(glmer(formula = Formula_glmer,
                                                                               data = data_itt,
                                                                               family = poisson()),silent = TRUE)))
  sum_fit = summary(fit_poisslog)
  
  if(inherits(res, "try-error")==FALSE){
    
    # Margin mean
    Trt_data <- NoTrt_data <- data_itt
    Trt_data$Arm <- 1 
    NoTrt_data$Arm <- 0
    
    Pred_Trt = predict(object = fit_poisslog,newdata = Trt_data)
    Pred_NoTrt = predict(object = fit_poisslog,newdata = NoTrt_data)
    
    P_Trt_est_poisslog = mean(exp(Pred_Trt))
    P_NoTrt_est_poisslog = mean(exp(Pred_NoTrt))
    
    RD_poisslog = P_Trt_est_poisslog - P_NoTrt_est_poisslog
    RD_itt_poisslog = c(RD_itt_poisslog,RD_poisslog)
    
    #Degree of freedom Corrected
    
    df_BW = 2*k - (nb_cov_tot + 1 + 1) # Number of cluster - number of fixed effect (nb_cov_tot + Arm + Intercept)
    
    se_tot_poisslog = fun_se_delta_meth_log_glmer(fit_poisslog,nb_cov,nb_cov_clus,data_itt)
    se_tot_poisslog = as.numeric(se_tot_poisslog)
    SE_itt_tot_poisslog = c(SE_itt_tot_poisslog,se_tot_poisslog)
    
    # CI RD 
    LL_CI_RD_poisslog = RD_poisslog - qt(0.975,(2*k-2))*se_tot_poisslog
    LL_CI_RD_itt_poisslog = c( LL_CI_RD_itt_poisslog, LL_CI_RD_poisslog)
    
    UL_CI_RD_poisslog = RD_poisslog + qt(0.975,(2*k-2))*se_tot_poisslog
    UL_CI_RD_itt_poisslog =c(UL_CI_RD_itt_poisslog,UL_CI_RD_poisslog)
    
    ## Results Poiss log ----
    # RD_mean_poisslog = mean(RD_itt_poisslog)
    itt_without_error_poisslog = length(RD_itt_poisslog)
    Coverage_rate_poisslog = length(which(RD_th > LL_CI_RD_itt_poisslog & RD_th < UL_CI_RD_itt_poisslog))/length(RD_itt_poisslog)*100
    Absolute_Bias_itt_poisslog = abs(RD_th - RD_itt_poisslog)
    Bias_itt_poisslog = RD_itt_poisslog - RD_th 
    Relative_Bias_itt_poisslog = RD_itt_poisslog/RD_th
    res_poisslog = data_frame(Risk_difference_Estimate_poisslog = RD_itt_poisslog,
                              Number_of_itteration_poisslog = itt_without_error_poisslog,
                              LL_95 = LL_CI_RD_itt_poisslog,
                              UL_95 = UL_CI_RD_itt_poisslog,
                              SE_itt = SE_itt_tot_poisslog,
                              Coverage_rate_poisslog = Coverage_rate_poisslog,
                              Abs_Bias_iteration_poisslog = Absolute_Bias_itt_poisslog,
                              Bias_iteration_poisslog = Bias_itt_poisslog,
                              Relative_bias_iteration_poisslog = Relative_Bias_itt_poisslog,
                              itt_para = itt_para,
                              type_error = sum_fit$optinfo$conv$opt)
    
    res_poisslog[unique(which(is.na(res_poisslog) | abs(res_poisslog$Risk_difference_Estimate_poisslog)>=1 | res_poisslog$type_error != 0 ,arr.ind = TRUE)[,1]),1:9] = NA
    
  }else{
    res_poisslog = data_frame(Risk_difference_Estimate_poisslog = NA,
                              Number_of_itteration_poisslog = NA,
                              LL_95 = NA,
                              UL_95 = NA,
                              SE_itt = NA,
                              Coverage_rate_poisslog = NA,
                              Abs_Bias_iteration_poisslog = NA,
                              Bias_iteration_poisslog = NA,
                              Relative_bias_iteration_poisslog = NA,
                              itt_para = itt_para,
                              type_error = "error")
    
  }
  
  
  return(res_poisslog)
  
}

# 15) Function use in the parallelism for GLMM  ----

# fun_para_analyse_%method_glmer : Are the function used during the parallelism
# They download the data set corresponding to the scenario and iteration corresponding during parallelism 
# and stock the result of one scenario in a Excel file created in advance but we will see that in R files number 3 and 4 of the whole simulation study.

# Input 
# itt_para       : Is the iteration of used in the parallelism mechanism, to take into account of a specific iteration (In our simulation study 1 to 1000)
# n              : Scenario number used
# Scenario_use   : List of the parameters of the scenario 'n' used
# Data_file      : Path where the data are saved
# G_res_file     : Path where the results will be save for the Gaussian distribution
# B_res_file     : Path where the results will be save for the binomial distribution
# P_res_file     : Path where the results will be save for the Poisson distribution

# Output
# Write and stock the result for one iterations ('itt_para') in a specific Excel file for the scenario 'n'


fun_para_analyse_gauss_glmer <- function(itt_para,n,Scenario_use,Data_file,G_res_file){
  
  data = read.csv(paste(paste(Data_file,"/Scenario_",n,sep=""),sep ="",paste("/Data_itt_",sep = "",itt_para,"_scen_",n,".csv")), sep=";")
  
  a = fun_sim_gauss_glmer(data,itt_para,Scenario = Scenario_use)
  
  write.table(a,file = paste(G_res_file,paste("/Data_output_Gauss_id_Scenario_",sep = "",n,".csv"),sep = ""),append = TRUE,row.names = FALSE,col.names = FALSE,sep = ";",dec = ",")
}

fun_para_analyse_binom_glmer <- function(itt_para,n,Scenario_use,Data_file,B_res_file){
  
  data = read.csv(paste(paste(Data_file,"/Scenario_",n,sep=""),sep ="",paste("/Data_itt_",sep = "",itt_para,"_scen_",n,".csv")), sep=";")
  
  a = fun_sim_binom_glmer(data,itt_para,Scenario = Scenario_use)
  
  write.table(a,file = paste(B_res_file,paste("/Data_output_Binom_logit_Scenario_",sep = "",n,".csv"),sep = ""),append = TRUE,row.names = FALSE,col.names = FALSE,sep = ";",dec = ",")
}

fun_para_analyse_poiss_glmer <- function(itt_para,n,Scenario_use,Data_file,P_res_file){
  
  data = read.csv(paste(paste(Data_file,"/Scenario_",n,sep=""),sep ="",paste("/Data_itt_",sep = "",itt_para,"_scen_",n,".csv")), sep=";")
  
  a = fun_sim_poiss_glmer(data,itt_para,Scenario = Scenario_use)
  
  write.table(a,file = paste(P_res_file,paste("/Data_output_Poiss_log_Scenario_",sep = "",n,".csv"),sep = ""),append = TRUE,row.names = FALSE,col.names = FALSE,sep = ";",dec = ",")
}


save.image("C:/Users/pereiramacedo/Desktop/2023_06_23_Article_1_V1/Papier/Code_github/WS_functions.RData")
