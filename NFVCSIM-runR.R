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
# set.seed(1235)
source("JOE_estimation.R")
source("JOE_grad_hessian.R")
source("JOE_inference.R")
source("NFVCSIM_functions.R")
source("FVCSIM_functions.R")
Rcpp::sourceCpp("~/Documents/NFVCSIM/Code_NFVCSIM/NFVCSIM_code/NFVCSIM_functions_rcpp.cpp")
Rcpp::sourceCpp("~/Documents/NFVCSIM/Code_NFVCSIM/NFVCSIM_code/FVCSIM_functions_rcpp.cpp")
source("~/Documents/NFVCSIM/Code_NFVCSIM/NFVCSIM_code/net_type_functions.R")
# source("~/Documents/NFVCSIM/Code_NFVCSIM/NFVCSIM_code/NFVC_func.R")

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
# NFVC_esti_coef1 = rep(list(matrix(0,nrow=max.points,ncol=3)),length(Nrep))
# NFVC_esti_coef2 = rep(list(matrix(0,nrow=max.points,ncol=3)),length(Nrep))
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
      
      # dat = get_func_W_network_rcpp(n,max.points,"3",Sige,eps,W)
      dat = get_func_2g_W_network_rcpp(n,max.points,"3",Sige,eps,W)
      # dat = get_func_network_rcpp(n,max.points,"3",Sige,eps)
      # W = dat$W
      true_g[[i]] = matrix(dat$g,ncol=2)
    
      # YK=get.YK(dat$Yna,dat$Tpoints)
      
      result = NFVCSIM_iter_rcpp(5,1e-5,dat$Y,dat$X,dat$z,W,D0,B0,Sige,g0,beta0,max.points,verbose=FALSE)
      d_alpha[[i]] = result$d_alpha
      g[[i]] = result$g
      beta[[i]] = result$beta

      FVCSIM_coef = FVCSIM_iter_rcpp(2,1e-4,dat$Y,dat$X,dat$z,g0,beta0,max.points,verbose=FALSE)
      FVCSIM_alpha[[i]] = FVCSIM_coef$alpha
      FVCSIM_g [[i]] = FVCSIM_coef$g
      FVCSIM_beta[[i]] = FVCSIM_coef$beta
      
      # t_points <- seq(0, 1, length.out = max.points)
      # Tpoints <- matrix(rep(t_points, each = n), nrow = n, ncol = max.points)
      # 
      # 
      # NFVC_coef1 = est_entire_timeline(dat$Y1,dat$X,Tpoints,dat$W,num.points=max.points)
      # NFVC_coef2 = est_entire_timeline(dat$Y2,dat$X,Tpoints,dat$W,num.points=max.points)
      # NFVC_esti_coef1[[i]] = NFVC_coef1$est_coef
      # NFVC_esti_coef2[[i]] = NFVC_coef2$est_coef
      
      break
    }, error=function(e){
      cat("⚠️ Error in iteration", i,":",conditionMessage(e),"\n")
      cat("Retrying...\n")
    })
    
  }
}
all_result = list(d_alpha=d_alpha,true_g=true_g,g=g,beta=beta,FVCSIM_beta=FVCSIM_beta,FVCSIM_g=FVCSIM_g,FVCSIM_alpha=FVCSIM_alpha)
# all_result = list(d_alpha=d_alpha,true_g=true_g,g=g,beta=beta,MNFVC_d_beta=MNFVC_d_beta)
filename = paste("NFVCSIM-W5-gamma-2g-Sige3(-0.3)-", n, "-", max.points, ".rda", sep = "")
save(all_result, file = filename)

# NFVC_result = list(NFVC_esti_coef1 = NFVC_esti_coef1, NFVC_esti_coef2 = NFVC_esti_coef2)
# filename = paste("NFVC-2g-W5-gamma-Sige3-", n, "-", max.points, ".rda", sep = "")
# save(NFVC_result,file = filename)

# esti_cov = est_covt(dat$Y,dat$X,dat$z,tt,dat$W,esti_d_alpha,esti_g,esti_beta)
# true_cov = getCovCST(tt)
# load("~/Documents/NFVCSIM/Code_NFVCSIM/NFVCSIM_code/NFVCSIM-W5-t3-2g-Sige1-100-50.rda")
# d_alpha=all_result$d_alpha
# g=all_result$g
# beta=all_result$beta
# # CST11=all_result$CST11
# # CST22=all_result$CST22
# FVCSIM_alpha = all_result$FVCSIM_alpha
# FVCSIM_g = all_result$FVCSIM_g
# FVCSIM_beta = all_result$FVCSIM_beta
NFVC_d_beta1 = NFVC_result$NFVC_esti_coef1
NFVC_d_beta2 = NFVC_result$NFVC_esti_coef2
NFVC_d_beta = rep(list(matrix(0,nrow=max.points,ncol=6)),length(Nrep))
for(i in 1:Nrep){
  NFVC_d_beta[[i]] = cbind(NFVC_d_beta1[[i]],NFVC_d_beta2[[i]])
}
###esti value
# esti_d_alpha = apply(simplify2array(d_alpha),1:2,median)
# esti_g = apply(simplify2array(g),1:2,median)
# esti_g1 = apply(simplify2array(g),1:2,median)[,1]
# esti_g2 = apply(simplify2array(g),1:2,median)[,2]
# esti_beta = apply(simplify2array(beta),1:2,median,na.rm=TRUE)
# esti_FVCSIM_alpha = apply(simplify2array(FVCSIM_alpha),1:2,median)
# esti_FVCSIM_g = apply(simplify2array(FVCSIM_g),1:2,median)
# esti_FVCSIM_beta = apply(simplify2array(FVCSIM_beta),1:2,median)
# esti_MNFVC_d_alpha = apply(simplify2array(MNFVC_d_beta),1:2,median)
# esti_cov11 = apply(simplify2array(CST11),1:2,mean)
# esti_cov22 = apply(simplify2array(CST22),1:2,mean)

esti_NFVC = apply(simplify2array(NFVC_d_beta),1:2,median)


###true value
tt = seq(0,1,length=max.points)
# true_d_alpha = cbind(0.2*sin(2*tt*pi)+0.3,2*(tt-tt^2),2*(tt^2-tt)+0.5, 0.2*cos(2*tt*pi)+0.3,
#                      tt^2, 5*(tt-0.5)^2, (1-tt)^2, tt^0.5)
# FVCSIM_true_alpha = cbind(tt^2, 5*(tt-0.5)^2, (1-tt)^2, tt^0.5)
true_d_alpha = cbind(0.2*sin(2*tt*pi)+0.3, tt^2, 5*(tt-0.5)^2, 
                     0.2*cos(2*tt*pi)+0.3,(1-tt)^2, tt^0.5)


IMSE_NFVC_coef =colMeans((true_d_alpha-esti_NFVC)^2)
print(c('d11','b11','b21','d22', 'b12', 'b22'))
r1 = IMSE_NFVC_coef
round(r1, 4)
SUP_NFVC_coef = apply(((true_d_alpha-esti_NFVC)^2),2,max)
round(SUP_NFVC_coef,4)

# beta1 = c(2,-1,3)/sqrt(14)
# beta2 = c(1,2,-3)/sqrt(14)
# true_beta=cbind(beta1,beta2)

# true_gg = plot_g(dat$Y,dat$X,dat$z,W,true_d_alpha,true_beta)
# esti_gg = plot_g(dat$Y,dat$X,dat$z,W,esti_d_alpha,esti_beta)
# zz=seq(-2,2,length=20)
# par(mfrow=c(1,2))
# plot(zz,true_gg$p_g1,xlab="z",ylab="g1")
# lines(zz,esti_gg$p_g1)
# plot(zz,true_gg$p_g2,xlab="z",ylab="g2")
# lines(zz,esti_gg$p_g2)

# zb1=z%*%beta[,1]
# true_g1=cos(zb1)
# plot(sort(zb1),true_g1[order(zb1)],type="l")
# lines(esti_g$p_g1)

###IMSE
# IMSE_d_alpha =colMeans((true_d_alpha-esti_d_alpha)^2)
# print(c('d11', 'd21','d12','d22', 'alpha11','alpha21','alpha12', 'alpha22'))
# IMSE_g1 = mean((true_g1-esti_g1)^2)
# IMSE_g2 = mean((true_g2-esti_g2)^2)
MSE_beta = sqrt((true_beta-esti_beta)^2)
round(MSE_beta,4)

esti_beta
true_beta


par(mfrow=c(2,4))
plot(tt,true_d_alpha[,1],type="l",col="black")
lines(tt,esti_d_alpha[,1],col="blue",lty=2)
plot(tt,true_d_alpha[,2],type="l",col="black")
lines(tt,esti_d_alpha[,2],col="blue",lty=2)
plot(tt,true_d_alpha[,3],type="l",col="black")
lines(tt,esti_d_alpha[,3],col="blue",lty=2)
plot(tt,true_d_alpha[,4],type="l",col="black")
lines(tt,esti_d_alpha[,4],col="blue",lty=2)
plot(tt,true_d_alpha[,5],type="l",col="black")
lines(tt,esti_d_alpha[,5],col="blue",lty=2)
lines(tt,esti_FVCSIM_alpha[,1],type="l",col="green",lty=2)
plot(tt,true_d_alpha[,6],type="l",col="black")
lines(tt,esti_d_alpha[,6],col="blue",lty=2)
lines(tt,esti_FVCSIM_alpha[,2],type="l",col="green",lty=2)
plot(tt,true_d_alpha[,7],type="l",col="black")
lines(tt,esti_d_alpha[,7],col="blue",lty=2)
lines(tt,esti_FVCSIM_alpha[,3],type="l",col="green",lty=2)
plot(tt,true_d_alpha[,8],type="l",col="black")
lines(tt,esti_d_alpha[,8],col="blue",lty=2)
lines(tt,esti_FVCSIM_alpha[,4],type="l",col="green",lty=2)


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


################IMSE_SUP#################

IMSE_coef=rep(list(matrix(0, nrow = max.points, ncol = 8)), length(Nrep))
for(i in 1:Nrep){
  IMSE_coef[[i]] =apply(((true_d_alpha-d_alpha[[i]])^2),2,mean)
}
IMSE_d_alpha = apply(do.call(rbind,IMSE_coef),2,mean)
print(c('d11', 'd21','d12','d22', 'alpha11','alpha21','alpha12', 'alpha22'))
r1 = round(IMSE_d_alpha,4)
r1
IMSE_sd = apply(do.call(rbind,IMSE_coef),2,sd)
round(IMSE_sd,4)

# 
# IMSE_coef=rep(list(matrix(0, nrow = max.points, ncol = 8)), length(Nrep))
# for(i in 1:Nrep){
#   IMSE_coef[[i]] =apply(((true_d_alpha-MNFVC_d_beta[[i]])^2),2,mean)
# }
# IMSE_d_alpha = apply(do.call(rbind,IMSE_coef),2,mean)
# print(c('d11', 'd21','d12','d22', 'alpha11','alpha21','alpha12', 'alpha22'))
# r1 = round(IMSE_d_alpha,4)
# r1
# IMSE_sd = apply(do.call(rbind,IMSE_coef),2,sd)
# round(IMSE_sd,4)



SUP_coef=rep(list(matrix(0, nrow = max.points, ncol = 8)), length(Nrep))
for(i in 1:Nrep){
  SUP_coef[[i]] =apply(((true_d_alpha-d_alpha[[i]])^2),2,max)
}
SUP_d_alpha = apply(do.call(rbind,SUP_coef),2,mean)
print(c('d11', 'd21','d12','d22', 'alpha11','alpha21','alpha12', 'alpha22'))
r2 = round(SUP_d_alpha,4)
r2
SUP_sd = apply(do.call(rbind,SUP_coef),2,sd)
round(SUP_sd,4)

FVCSIM_IMSE = rep(list(matrix(0, nrow = max.points, ncol = 4)), length(Nrep))
for(i in 1:Nrep){
  FVCSIM_IMSE[[i]] = apply(((FVCSIM_true_alpha-FVCSIM_alpha[[i]])^2),2,mean)
}
FVCSIM_IMSE_coef = apply(do.call(rbind,FVCSIM_IMSE),2,mean)
round(FVCSIM_IMSE_coef,4)
FVCSIM_IMSE_sd = apply(do.call(rbind,FVCSIM_IMSE),2,sd)
round(FVCSIM_IMSE_sd,4)

FVCSIM_SUP = rep(list(matrix(0, nrow = max.points, ncol = 4)), length(Nrep))
for(i in 1:Nrep){
  FVCSIM_SUP[[i]] = apply(((FVCSIM_true_alpha-FVCSIM_alpha[[i]])^2),2,max)
}
FVCSIM_SUP_coef = apply(do.call(rbind,FVCSIM_SUP),2,median)
round(FVCSIM_SUP_coef,4)
FVCSIM_SUP_sd = apply(do.call(rbind,FVCSIM_SUP),2,sd)
round(FVCSIM_SUP_sd,4)

MSE_beta = rep(list(rep(0,6)),length(Nrep))
for(i in 1:Nrep){
  MSE_beta[[i]]=(as.vector(beta[[i]])-as.vector(true_beta))^2
}
MSE_b = sqrt(colMeans(do.call(rbind,MSE_beta)))
r3 =round(MSE_b,4)
r3
MSE_b_sd = apply(do.call(rbind,MSE_beta),2,sd)
round(MSE_b_sd,4)



FVCSIM_MSE_beta = rep(list(rep(0,6)),length(Nrep))
for(i in 1:Nrep){
  FVCSIM_MSE_beta[[i]]=(as.vector(FVCSIM_beta[[i]])-as.vector(true_beta))^2
}
FVCSIM_MSE_b = sqrt(colMeans(do.call(rbind,FVCSIM_MSE_beta)))
r3 =round(FVCSIM_MSE_b,4)
r3
FVCSIM_MSE_b_sd = apply(do.call(rbind,FVCSIM_MSE_beta),2,sd)
round(FVCSIM_MSE_b_sd,4)


# MSE_beta =sqrt((as.vector(true_beta)-as.vector(esti_beta))^2)
# round(MSE_beta, 4)
# 
# MSE_FVCSIM_beta = sqrt((as.vector(true_beta)-as.vector(esti_FVCSIM_beta))^2)
# round(MSE_FVCSIM_beta,4)
# 
# 

#########################plot
library(ggplot2)
library(patchwork)
library(cowplot)

plot_one_curve <- function(t_list, y_list, labels, title,y_limits = NULL,y_breaks=NULL,show_legend=FALSE) {
  # 全部可能的标签
  all_labels <- c("True", "NFVCSIM", "NFVCM", "FVCSIM")
  
  df <- do.call(rbind, lapply(seq_along(y_list), function(i) {
    data.frame(
      t = t_list[[i]],
      y = y_list[[i]],
      Source = factor(labels[i], levels = all_labels)
    )
  }))
  
  p=ggplot(df, aes(x = t, y = y, color = Source, linetype = Source)) +
    geom_line(size = 0.7, na.rm = TRUE) +
    labs(title = title, x = "", y = "") +
    scale_color_manual(values = c(
      # "True" = "#B3B3B3",
      # "NFVCSIM" = "#E9A3C9",
      # "NFVCM" = "#80B1D3",
      # "FVCSIM" ="#B3DE69"
      "True"="#B3B3B3",
      "NFVCSIM"="#FC4E2A",
      "NFVCM"="#3690C0",
      "FVCSIM"="#B3DE69"
    )) +
    scale_linetype_manual(values = c(
      "True" = "solid", 
      "NFVCSIM" = "longdash", 
      "NFVCM" = "dotdash",
      "FVCSIM" = "twodash"
    )) +
    theme_minimal(base_size = 10) +
    theme(
      legend.title = element_blank(),
      plot.title = element_text(size = 10, hjust = 0.5),
      axis.title = element_text(size = 9),
      axis.text = element_text(size = 8),
      panel.grid=element_blank(),
      panel.border = element_rect(color="black",fill=NA),
      axis.line = element_line(color="black")
    )
  if (!show_legend) {
    p <- p + theme(legend.position = "none")
  }
  
  if (!is.null(y_limits)) {
    p <- p + coord_cartesian(ylim = y_limits)
  }
  if (!is.null(y_breaks)) {
    p <- p + scale_y_continuous(breaks = y_breaks, labels = sprintf("%.2f", y_breaks))
  }
  
  
  return(p)
}

d11 <- plot_one_curve(list(t, t), list(true_d_alpha[,1], esti_d_alpha[,1]), c("True", "NFVCSIM"), expression(hat(d)[11](t)),y_limits = c(0,0.5),show_legend = FALSE)
d21 <- plot_one_curve(list(t, t), list(true_d_alpha[,2], esti_d_alpha[,2]), c("True", "NFVCSIM"), expression(hat(d)[21](t)),y_limits=c(0,0.5),show_legend = FALSE,y_breaks = seq(0, 0.5, by = 0.1))
d12 <- plot_one_curve(list(t, t), list(true_d_alpha[,3], esti_d_alpha[,3]), c("True", "NFVCSIM"), expression(hat(d)[12](t)),y_limits=c(0,0.5),show_legend=FALSE,y_breaks = seq(0, 0.5, by = 0.1))
d22 <- plot_one_curve(list(t, t), list(true_d_alpha[,4], esti_d_alpha[,4]), c("True", "NFVCSIM"), expression(hat(d)[22](t)),y_limits=c(0,0.5),show_legend = FALSE)

d_plot = wrap_plots(d11,d21,d12,d22,ncol=4,guides = "collect")&
  theme(legend.position = "bottom",
        legend.key.size = unit(1.5, "lines"),  # 放大图例键
        legend.text = element_text(size = 14))  # 放大图例文本
d_plot

a11 <- plot_one_curve(list(t, t), list(true_d_alpha[,5], esti_d_alpha[,5]), c("True", "NFVCSIM"), expression(hat(a)[11](t)),y_limits = c(0,1),show_legend = FALSE)
a21 <- plot_one_curve(list(t, t), list(true_d_alpha[,6], esti_d_alpha[,6]), c("True", "NFVCSIM"), expression(hat(a)[21](t)),y_limits=c(0,1),show_legend = FALSE)
a12 <- plot_one_curve(list(t, t), list(true_d_alpha[,7], esti_d_alpha[,7]), c("True", "NFVCSIM"), expression(hat(a)[12](t)),y_limits=c(0,1),show_legend=FALSE)
a22 <- plot_one_curve(list(t, t), list(true_d_alpha[,8], esti_d_alpha[,8]), c("True", "NFVCSIM"), expression(hat(a)[22](t)),y_limits=c(0,1),show_legend = FALSE)

a_plot = wrap_plots(a11,a21,a12,a22,ncol=4,guides = "collect")&
  theme(legend.position = "bottom",
        legend.key.size = unit(1.5, "lines"),  # 放大图例键
        legend.text = element_text(size = 14))  # 放大图例文本
a_plot

# eigenvalue1=eigen(esti_cov11)$values
# eigenvectors1 = eigen(esti_cov11)$vectors
# 
# eigenvalue2=eigen(esti_cov22)$values
# eigenvectors2 = eigen(esti_cov22)$vectors
# 
# eigenvalue1
# eigenvalue2
# eigenvectors1
# eigenvectors2
# 
# true_CST=getCovCST(t)
# true_eigenvalues1 = eigen(true_CST$cov11)$values
# true_eigenvectors1 = eigen(true_CST$cov11)$vectors
# true_eigenvalues2 = eigen(true_CST$cov22)$values
# true_eigenvectors2 = eigen(true_CST$cov22)$vectors
# 
# true_eigenvalues1
# true_eigenvalues2
# true_eigenvectors1
# true_eigenvectors2
# 
# plot(t,true_eigenvectors1[,1])
# lines(t,eigenvectors1[,1])
# plot(t,true_eigenvectors1[,2])
# lines(t,eigenvectors1[,2])
# plot(t,true_eigenvectors2[,1])
# lines(t,eigenvectors2[,1])
# plot(t,true_eigenvectors2[,2])
# lines(t,eigenvectors2[,2])
# 
# EUC_dis=norm(true_CST$cov11-esti_cov11,"F")
# EUC_dis


