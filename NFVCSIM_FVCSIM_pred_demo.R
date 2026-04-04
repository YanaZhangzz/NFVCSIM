rm(list=ls())
setwd("~/Documents/NFVCSIM/Code_NFVCSIM/NFVCSIM_code")
library(Matrix)
library(poweRlaw)
library(MultiRNG)
library(MASS)
library(mvtnorm)
library(KernSmooth)
library(VGAM)
library(methods)
Rcpp::sourceCpp("~/Documents/NFVCSIM/Code_NFVCSIM/NFVCSIM_code/NFVCSIM_functions_rcpp.cpp")
Rcpp::sourceCpp("~/Documents/NFVCSIM/Code_NFVCSIM/NFVCSIM_code/NFVCSIM_pred_functions_rcpp.cpp")
# Rcpp::sourceCpp("~/Documents/NFVCSIM/Code_NFVCSIM/NFVCSIM_code/FVCSIM_pred_functions_rcpp.cpp")
source("~/Documents/NFVCSIM/Code_NFVCSIM/NFVCSIM_code/NFVC_func.R")
Rcpp::sourceCpp("~/Documents/NFVCSIM/Code_NFVCSIM/NFVCSIM_code/NFVC_pred_func.cpp")
source("~/Documents/NFVCSIM/Code_NFVCSIM/NFVCSIM_code/net_type_functions.R")

n = 200
max.points = 100
Nrep=100
D0 = Matrix(0, 2, 2)
B0 = Matrix(0, 2, 2)
Sige1 = Diagonal(n = 2, x = c(0.4,0.6))
Sige1[1,2] = 0.1
Sige1[2,1] = 0.1
Sige=Sige1
# Sige2 = Diagonal(n = 2, x = c(0.4,0.6))
# Sige2[1,2] = 0
# Sige2[2,1] = 0
# Sige = Sige2
# Sige3 = Diagonal(n = 2, x = c(0.4,0.6))
# Sige3[1,2] = 0.45
# Sige3[2,1] = 0.45
# Sige = Sige3
# D0 = as.matrix(D0)
# B0 = as.matrix(B0)
Sige = as.matrix(Sige)

# beta10=c(0.5,-0.3,0.8)
# beta20=c(0.3,0.5,-0.8)
# beta0 = cbind(beta10,beta20)



# NFVCSIM_pred = list(test_IMSE=numeric(Nrep),
#                     train_IMSE = numeric(Nrep),
#                     test_R2 = numeric(Nrep),
#                     train_R2 = numeric(Nrep))
# FVCSIM_pred = list(test_IMSE=numeric(Nrep),
#                    train_IMSE = numeric(Nrep),
#                    test_R2 = numeric(Nrep),
#                    train_R2 = numeric(Nrep))

NFVC_pred = list(test_IMSE=numeric(Nrep),
                    train_IMSE = numeric(Nrep),
                    test_R2 = numeric(Nrep),
                    train_R2 = numeric(Nrep))




for (i in 1:Nrep){
  
  repeat{
    tryCatch({
      cat(i,"\r")
      
      # eps = matrix(rnorm(n*2*max.points),nrow=n*2)
      # eps = matrix(rt(n*2*max.points,df=3),nrow=n*2)
      eps = matrix(rgamma(n*2*max.points,shape=1.5,scale=2),nrow=n*2)
      eps = eps-mean(eps)
      
      # W = create_W3(10,20)
      W = create_W5(10,20)

      dat = get_predfunc_W_network_rcpp(n,max.points,"3",Sige,eps,W)
      # dat = get_predfunc_network_rcpp(n,max.points,"3",Sige,eps)
      t_points <- seq(0, 1, length.out = max.points)
      train_Tpoints <- matrix(rep(t_points, each = n*0.8), nrow = n*0.8, ncol = max.points)
      
      NFVC_esti_pred = NFVC_pred_func_rcpp(est_entire_timeline,dat,max.points,train_Tpoints,verbose=FALSE)
      
      NFVC_pred$test_IMSE[i] = NFVC_esti_pred$test_IMSE
      NFVC_pred$train_IMSE[i] = NFVC_esti_pred$train_IMSE
      NFVC_pred$test_R2[i] = NFVC_esti_pred$test_R2
      NFVC_pred$train_R2[i] = NFVC_esti_pred$train_R2
      
      
      # NFVCSIM_esti_pred = NFVCSIM_pred_func_rcpp(dat,max.points,D0,B0,Sige,beta0,verbose=FALSE)
      # 
      # NFVCSIM_pred$test_IMSE[i] = NFVCSIM_esti_pred$test_IMSE
      # NFVCSIM_pred$train_IMSE[i] = NFVCSIM_esti_pred$train_IMSE
      # NFVCSIM_pred$test_R2[i] = NFVCSIM_esti_pred$test_R2
      # NFVCSIM_pred$train_R2[i] = NFVCSIM_esti_pred$train_R2
      # # train_coef[i] = NFVCSIM_esti_pred$train_coef
      # 
      # FVCSIM_esti_pred = FVCSIM_pred_func_rcpp(dat,max.points,beta0,verbose=FALSE)
      # 
      # FVCSIM_pred$test_IMSE[i] = FVCSIM_esti_pred$test_IMSE
      # FVCSIM_pred$train_IMSE[i] = FVCSIM_esti_pred$train_IMSE
      # FVCSIM_pred$test_R2[i] = FVCSIM_esti_pred$test_R2
      # FVCSIM_pred$train_R2[i] = FVCSIM_esti_pred$train_R2
      
      
  
  break
    }, error=function(e){
      cat("⚠️ Error in iteration", i,":",conditionMessage(e),"\n")
      cat("Retrying...\n")
    })
    
  }

}

all_result = list(NFVC_results = NFVC_pred)

filename = paste("NFVC-W5-pred-gamma-Sige1-", n, "-", max.points, ".rda", sep = "")
save(all_result, file = filename)
# load("~/Documents/NFVCSIM/Code_NFVCSIM/NFVCSIM_code/NFVCSIM-3-pred-norm-alpha1.8-Sige1-50-25.rda")
NFVC_final_results <- list(
  test_IMSE = round(median(all_result$NFVC_results$test_IMSE),4),
  train_IMSE = round(median(all_result$NFVC_results$train_IMSE),4),
  test_R2 = round(median(all_result$NFVC_results$test_R2),4),
  train_R2 = round(median(all_result$NFVC_results$train_R2),4)
)
NFVC_final_results
# NFVCSIM_final_results <- list(
#   test_IMSE = round(median(all_result$NFVCSIM_results$test_IMSE),4),
#   train_IMSE = round(median(all_result$NFVCSIM_results$train_IMSE),4),
#   test_R2 = round(median(all_result$NFVCSIM_results$test_R2),4),
#   train_R2 = round(median(all_result$NFVCSIM_results$train_R2),4)
# )
# 
# FVCSIM_final_results <- list(
#   test_IMSE = round(median(all_result$FVCSIM_results$test_IMSE),4),
#   train_IMSE = round(median(all_result$FVCSIM_results$train_IMSE),4),
#   test_R2 = round(median(all_result$FVCSIM_results$test_R2),4),
#   train_R2 = round(median(all_result$FVCSIM_results$train_R2),4)
# )
# 
# NFVCSIM_final_results
# FVCSIM_final_results
# 
# 
