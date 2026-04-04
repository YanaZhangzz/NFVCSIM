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
library(np)
source("JOE_estimation.R")
source("JOE_grad_hessian.R")
source("JOE_inference.R")
source("NFVCSIM_functions.R")
source("FVCSIM_functions.R")
source("~/Documents/NFVCSIM/Code_NFVCSIM/NFVCSIM_code/net_type_functions.R")
source("~/Documents/NFVCSIM/Code_NFVCSIM/NFVCSIM_code/NFVC_func.R")

n = 100
max.points = 50

Nrep = 100

d_alpha = rep(list(matrix(0, nrow = max.points, ncol = 8)), length(Nrep))
g = rep(list(matrix(0, nrow = n, ncol = 4)), length(Nrep))
beta = rep(list(matrix(0, nrow = 3, ncol = 2)), length(Nrep))
true_g = rep(list(matrix(0, nrow = n, ncol = 2)), length(Nrep))
FVCSIM_alpha = rep(list(matrix(0,nrow = max.points,ncol=4)),length(Nrep))
FVCSIM_g =  rep(list(matrix(0, nrow = n, ncol = 4)), length(Nrep))
FVCSIM_beta = rep(list(matrix(0, nrow = 3, ncol = 2)), length(Nrep))
NFVC_esti_coef1 = rep(list(matrix(0,nrow=max.points,ncol=3)),length(Nrep))
NFVC_esti_coef2 = rep(list(matrix(0,nrow=max.points,ncol=3)),length(Nrep))
# # ###设置初始值
D0 = Matrix(0, 2, 2)
B0 = Matrix(0, 2, 2)
# Sige1 = Diagonal(n = 2, x = c(0.4,0.6))
# Sige1[1,2] = 0.1
# Sige1[2,1] = 0.1
# Sige = Sige1
# Sige2 = diag(c(0.4,0.6))
# Sige2[1,2] = 0
# Sige2[2,1] = 0
# Sige2 = as.matrix(Sige2)
# Sige = Sige2
Sige3 = Diagonal(n = 2, x = c(0.4,0.6))
Sige3[1,2] = -0.3
Sige3[2,1] = -0.3
Sige3 = as.matrix(Sige3)
Sige=Sige3
D0=as.matrix(D0)
B0=as.matrix(B0)
Sige=as.matrix(Sige)

g0=matrix(0.5,nrow=n,ncol=4)
beta10=c(0.5,-0.3,0.8)
beta20=c(0.3,0.5,-0.8)
beta0 = cbind(beta10,beta20)


for (i in 1:Nrep) {
  
  repeat{
    tryCatch({
      cat(i,"\r")
      
      # eps = matrix(rnorm(n*2*max.points),nrow=n*2)
      # eps = matrix(rt(n*2*max.points,df=3),nrow=n*2)
      eps = matrix(rgamma(n*2*max.points,shape=1.5,scale=2),nrow=n*2)
      eps = eps-mean(eps)

      # W = create_W3(10,10)
      W = create_W5(10,10)
      
      # dat = get_func_W_network(n,max.points,"3",Sige,eps,W)
      dat = get_func_2g_W_network(n,max.points,"3",Sige,eps,W)
      # dat = get_func_network(n,max.points,"3",Sige,eps)
      # W = dat$W
      true_g[[i]] = matrix(dat$g,ncol=2)
    
      # YK=get.YK(dat$Yna,dat$Tpoints)
      
      result = NFVCSIM_iter_rcpp(5,1e-5,dat$Y,dat$X,dat$z,W,D0,B0,Sige,g0,beta0,max.points,verbose=FALSE)
      d_alpha[[i]] = result$d_alpha
      g[[i]] = result$g
      beta[[i]] = result$beta

      FVCSIM_coef = FVCSIM_iter_rcpp(5,1e-5,dat$Y,dat$X,dat$z,g0,beta0,max.points,verbose=FALSE)
      FVCSIM_alpha[[i]] = FVCSIM_coef$alpha
      FVCSIM_g [[i]] = FVCSIM_coef$g
      FVCSIM_beta[[i]] = FVCSIM_coef$beta
      
      # t_points <- seq(0, 1, length.out = max.points)
      # Tpoints <- matrix(rep(t_points, each = n), nrow = n, ncol = max.points)
      # 
      # 
      NFVC_coef1 = est_entire_timeline(dat$Y1,dat$X,Tpoints,dat$W,num.points=max.points)
      NFVC_coef2 = est_entire_timeline(dat$Y2,dat$X,Tpoints,dat$W,num.points=max.points)
      NFVC_esti_coef1[[i]] = NFVC_coef1$est_coef
      NFVC_esti_coef2[[i]] = NFVC_coef2$est_coef
      
      break
    }, error=function(e){
      cat("⚠️ Error in iteration", i,":",conditionMessage(e),"\n")
      cat("Retrying...\n")
    })
    
  }
}
all_result = list(d_alpha=d_alpha,true_g=true_g,g=g,beta=beta,FVCSIM_beta=FVCSIM_beta,FVCSIM_g=FVCSIM_g,FVCSIM_alpha=FVCSIM_alpha)
filename = paste("NFVCSIM-W5-gamma-2g-Sige3(-0.3)-", n, "-", max.points, ".rda", sep = "")
save(all_result, file = filename)


NFVC_d_beta1 = NFVC_result$NFVC_esti_coef1
NFVC_d_beta2 = NFVC_result$NFVC_esti_coef2
NFVC_d_beta = rep(list(matrix(0,nrow=max.points,ncol=6)),length(Nrep))
for(i in 1:Nrep){
  NFVC_d_beta[[i]] = cbind(NFVC_d_beta1[[i]],NFVC_d_beta2[[i]])
}
###esti value
esti_d_alpha = apply(simplify2array(d_alpha),1:2,median)
esti_g = apply(simplify2array(g),1:2,median)
esti_g1 = apply(simplify2array(g),1:2,median)[,1]
esti_g2 = apply(simplify2array(g),1:2,median)[,2]
esti_beta = apply(simplify2array(beta),1:2,median,na.rm=TRUE)
esti_FVCSIM_alpha = apply(simplify2array(FVCSIM_alpha),1:2,median)
esti_FVCSIM_g = apply(simplify2array(FVCSIM_g),1:2,median)
esti_FVCSIM_beta = apply(simplify2array(FVCSIM_beta),1:2,median)
esti_MNFVC_d_alpha = apply(simplify2array(MNFVC_d_beta),1:2,median)
esti_cov11 = apply(simplify2array(CST11),1:2,mean)
esti_cov22 = apply(simplify2array(CST22),1:2,mean)

esti_NFVC = apply(simplify2array(NFVC_d_beta),1:2,median)


###true value
tt = seq(0,1,length=max.points)
true_d_alpha = cbind(0.2*sin(2*tt*pi)+0.3,2*(tt-tt^2),2*(tt^2-tt)+0.5, 0.2*cos(2*tt*pi)+0.3,
                     tt^2, 5*(tt-0.5)^2, (1-tt)^2, tt^0.5)
FVCSIM_true_alpha = cbind(tt^2, 5*(tt-0.5)^2, (1-tt)^2, tt^0.5)
true_d_alpha = cbind(0.2*sin(2*tt*pi)+0.3, tt^2, 5*(tt-0.5)^2, 
                     0.2*cos(2*tt*pi)+0.3,(1-tt)^2, tt^0.5)


IMSE_NFVC_coef =colMeans((true_d_alpha-esti_NFVC)^2)
print(c('d11','b11','b21','d22', 'b12', 'b22'))
r1 = IMSE_NFVC_coef
round(r1, 4)
SUP_NFVC_coef = apply(((true_d_alpha-esti_NFVC)^2),2,max)
round(SUP_NFVC_coef,4)

beta1 = c(2,-1,3)/sqrt(14)
beta2 = c(1,2,-3)/sqrt(14)
true_beta=cbind(beta1,beta2)



###IMSE
IMSE_d_alpha =colMeans((true_d_alpha-esti_d_alpha)^2)
print(c('d11', 'd21','d12','d22', 'alpha11','alpha21','alpha12', 'alpha22'))
IMSE_g1 = mean((true_g1-esti_g1)^2)
IMSE_g2 = mean((true_g2-esti_g2)^2)
MSE_beta = sqrt((true_beta-esti_beta)^2)
round(MSE_beta,4)

esti_beta
true_beta


IMSE_coef =colMeans((true_d_alpha-esti_d_alpha)^2)
print(c('d11', 'd21','d21','d22', 'b11','b21','b12', 'b22'))
r1 = IMSE_coef
round(r1, 4)
SUP_coef = apply(((true_d_alpha-esti_d_alpha)^2),2,max)
round(SUP_coef,4)

IMSE_coef =colMeans((FVCSIM_true_alpha-esti_FVCSIM_alpha)^2)
print(c('d11', 'd21','d21','d22', 'b11','b21','b12', 'b22'))
r2 = IMSE_coef
round(r2, 4)
SUP_coef = apply(((FVCSIM_true_alpha-esti_FVCSIM_alpha)^2),2,max)
round(SUP_coef,4)


MSE_beta = rep(list(rep(0,6)),length(Nrep))
for(i in 1:Nrep){
  MSE_beta[[i]]=(as.vector(beta[[i]])-as.vector(true_beta))^2
}
MSE_b = sqrt(colMeans(do.call(rbind,MSE_beta)))
r3 =round(MSE_b,4)
r3

FVCSIM_MSE_beta = rep(list(rep(0,6)),length(Nrep))
for(i in 1:Nrep){
  FVCSIM_MSE_beta[[i]]=(as.vector(FVCSIM_beta[[i]])-as.vector(true_beta))^2
}
FVCSIM_MSE_b = sqrt(colMeans(do.call(rbind,FVCSIM_MSE_beta)))
r4 =round(FVCSIM_MSE_b,4)
r4

