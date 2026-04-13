NFVCSIM_pred_func = function(dat)
{
  Y1 = dat$Y1
  Y2 = dat$Y2
  X = dat$X
  z = dat$z
  W = dat$W
  
  n = nrow(Y1)
  
  train_size = floor(4*n/5)
  train_grid = seq(1,train_size)
  test_grid = seq(train_size+1,n)
  
  train_Y1 = Y1[train_grid,]
  train_Y2 = Y2[train_grid,]
  train_Y = rbind(train_Y1,train_Y2)
  train_X = X[train_grid,]
  train_z = z[train_grid,]
  train_W = W[train_grid,train_grid]
  
  test_Y1 = Y1[test_grid,]
  test_Y2 = Y2[test_grid,]
  test_Y = rbind(test_Y1,test_Y2)
  test_X = X[test_grid,]
  test_z = z[test_grid,]
  test_W = W[test_grid,test_grid]
  
  g0_train = matrix(0.5,nrow=n*4/5,ncol=4)
  
  train_coef = NFVCSIM_iter_rcpp(8,1e-5,train_Y,train_X,train_z,train_W,D0,B0,Sige1,g0_train,beta0,max.points,verbose=FALSE)
  
  ans = train_coef$d_alpha
  
  test_g = est_abt_rcpp(test_Y,test_X,test_z,test_W,ans,train_coef$beta)
  test_eta = est_eta(test_Y,test_X,test_z,test_W,ans,test_g,train_coef$beta)
  
  d.mat = train_coef$d_alpha[,1:4]
  alpha.mat = train_coef$d_alpha[,-(1:4)]
  
  test_X_hat = kronecker(diag(2),test_X)
  test_g.mat = matrix(test_g[,1:2],nrow=n*2/5)
  test_eta.mat = rbind(test_eta[,1,],test_eta[,2,])
  
  test_Ymat = matrix(0,nrow=2*n/5,ncol=max.points)
  for(m in 1:max.points){
    Xalpha = test_X_hat%*%alpha.mat[m,]+test_g.mat+test_eta.mat[,m]
    D = matrix(d.mat[m,],nrow=2)
    DkW = kronecker(t(D),test_W)
    DkW2 = DkW%*%DkW
    DkW3 = DkW2%*%DkW
    DkW4 = DkW3%*%DkW
    DkW5 = DkW4%*%DkW
    DkW6 = DkW5%*%DkW
    
    Yvec = Xalpha+DkW%*%Xalpha+DkW2%*%Xalpha+DkW3%*%Xalpha+DkW4%*%Xalpha+DkW5%*%Xalpha+DkW6%*%Xalpha
    test_Ymat[,m] = Yvec
  }
  
  test_IMSE = mean(colMeans(test_Ymat-test_Y)^2)
  
  test_A = test_Y - test_Ymat
  test_SSE = sum(test_A^2)
  test_mean_testY = apply(test_Y,2,mean)
  test_B = test_Y-test_mean_testY
  test_SST = sum(test_B^2)
  test_R2 = 1-test_SSE/test_SST
  
  
  train_g = est_abt_rcpp(train_Y,train_X,train_z,train_W,ans,train_coef$beta)
  train_eta = est_eta(train_Y,train_X,train_z,train_W,ans,train_g,train_coef$beta)
  
  train_X_hat = kronecker(diag(2),train_X)
  train_g.mat = matrix(train_g[,1:2],nrow=n*2*4/5)
  train_eta.mat = rbind(train_eta[,1,],train_eta[,2,])
  
  train_Ymat = matrix(0,nrow=2*n*4/5,ncol=max.points)
  for(m in 1:max.points){
    Xalpha = train_X_hat%*%alpha.mat[m,]+train_g.mat+train_eta.mat[,m]
    D = matrix(d.mat[m,],nrow=2)
    DkW = kronecker(t(D),train_W)
    DkW2 = DkW%*%DkW
    DkW3 = DkW2%*%DkW
    DkW4 = DkW3%*%DkW
    DkW5 = DkW4%*%DkW
    DkW6 = DkW5%*%DkW
    
    Yvec = Xalpha+DkW%*%Xalpha+DkW2%*%Xalpha+DkW3%*%Xalpha+DkW4%*%Xalpha+DkW5%*%Xalpha+DkW6%*%Xalpha
    train_Ymat[,m] = Yvec
  }
  
  train_IMSE = mean(colMeans(train_Ymat-train_Y)^2)
  
  train_A = train_Y - train_Ymat
  train_SSE = sum(train_A^2)
  train_mean_trainY = apply(train_Y,2,mean)
  train_B = train_Y - train_mean_trainY
  train_SST = sum(train_B^2)
  train_R2 = 1-train_SSE/train_SST
  
  list(test_IMSE,train_IMSE,test_R2,train_R2)
  
}


est_et<-function(Y,X,z,t_est,W,ans,bns,beta){
  
  n = nrow(Y)/2
  Y1=Y[1:n,]
  Y2=Y[-(1:n),]
  p = ncol(X)
  t = seq(0,1,length=max.points)
  d11s.mat = matrix(rep(ans[,1],each=n),nrow=n,ncol=max.points)
  d12s.mat = matrix(rep(ans[,3],each=n),nrow=n,ncol=max.points)
  d21s.mat = matrix(rep(ans[,2],each=n),nrow=n,ncol=max.points)
  d22s.mat = matrix(rep(ans[,4],each=n),nrow=n,ncol=max.points)
  alpha1s.seq = ans[,5:(4+p)]
  alpha2s.seq = ans[,(5+p):(p*2+4)]
  Y1t = Y1- X%*%t(alpha1s.seq)-d11s.mat*W%*%Y1-d21s.mat*W%*%Y2
  Y2t = Y2- X%*%t(alpha2s.seq)-d22s.mat*W%*%Y2-d12s.mat*W%*%Y1
  Yt=rbind(Y1t,Y2t)
  
  est_a=bns[,1:2]
  est_b=bns[,3:4]
  
  Yi_1=Y1t-est_a[,1]
  Yi_2=Y2t-est_a[,2]
  
  beta1=beta[,1]
  beta2=beta[,2]
  
  h12 = rule.thumb(Y1t)
  h22 = rule.thumb(Y2t)
  
  k1_t=dnorm((t-t_est)/h12)/h12
  mx1_t=cbind(rep(1,max.points),(t-t_est)/h12)
  k2_t=dnorm((t-t_est)/h22)/h22
  mx2_t=cbind(rep(1,max.points),(t-t_est)/h22)
  
  ##esti eta_1
  sum1=matrix(0,nrow=2,ncol=2)
  for(m in 1:max.points){
    sum1=sum1+k1_t[m]*mx1_t[m,]%*%t(mx1_t[m,])
  }
  sum_inv1=solve(sum1)
  
  sum2=matrix(0,nrow=n,ncol=2)
  for(m in 1:max.points){
    sum2=sum2+k1_t[m]*Yi_1[,m]%*%t(mx1_t[m,])
  }
  
  eta1=sum2%*%sum_inv1
  
  ##esti eta_2
  sum3=matrix(0,nrow=2,ncol=2)
  for(m in 1:max.points){
    sum3=sum3+k2_t[m]*mx2_t[m,]%*%t(mx2_t[m,])
  }
  sum_inv2=solve(sum3)
  
  sum4=matrix(0,nrow=n,ncol=2)
  for(m in 1:max.points){
    sum4=sum4+k2_t[m]*Yi_2[,m]%*%t(mx2_t[m,])
  }
  
  eta2=sum4%*%sum_inv2
  
  return(cbind(eta1[,1],eta2[,1]))

}

est_eta<-function(Y,X,z,W,ans,bns,beta){
  
  n=nrow(Y)/2
  eta=array(0,dim=c(n,2,max.points))
  t=seq(0,1,length=max.points)
  for(s in 1:max.points){
    ts=seq(0,1,length=max.points)
    eta[,,s]=est_et(Y,X,z,ts[s],W,ans,bns,beta)
  }
  return(eta)
}
