####FVCSIM
rule.thumb_FVCSIM = function(x){
  0.79*min(quantile(x,0.75,na.rm=T)-quantile(x,0.25,na.rm=T),sd(x,na.rm=T))*length(x)^(-1/5)
}

FVCSIM_est_alphat=function(Y,X,z,s,bns,beta){
  n = nrow(Y)/2
  Y1 = Y[1:n,]
  Y2 = Y[-(1:n),]
  g1 = bns[,1]
  g2 = bns[,2]
  beta1=beta[,1]
  beta2=beta[,2]
  
  diff1 = Y1[,1]-z%*%beta1
  diff2 = Y2[,1]-z%*%beta2
  
  h11=rule.thumb_FVCSIM(diff1)
  h21=rule.thumb_FVCSIM(diff2)
  # h11=dpill(Y1,z%*%beta1)
  # h21=dpill(Y2,z%*%beta2)
  
  # D_im1 = c(X[i,],X[i,]%*%(s_m_s)/h11)
  # D_im2 = c(X[i,],X[i,]%*%(s_m_s)/h21)
  # 
  sm=seq(0,1,length=max.points)
  Wh1 = matrix(0,nrow=4,ncol=4)
  for(i in 1:n){
    for(m in 1:max.points){
      K1 = dnorm((sm[m]-s)/h11)/h11
      # D_im1 = rbind(matrix(rep(X[i,],max.points),nrow=2),X[i,]%*%(s_m_s)/h11)
      D_im1 = c(X[i,],X[i,]*(sm[m]-s)/h11)
      Wh1 = Wh1+ K1*D_im1%*%t(D_im1)
    }
  }
  
  Wh2 = matrix(0,nrow=4,ncol=4)
  for(i in 1:n){
    for(m in 1:max.points){
      K2 = dnorm((sm[m]-s)/h21)/h21
      D_im2 = c(X[i,],X[i,]*(sm[m]-s)/h21)
      Wh2 = Wh2+ K2*D_im2%*%t(D_im2)
    }
  }
  
  absum1 = matrix(0,nrow=4,ncol=1)
  for(i in 1:n){
    for(m in 1:max.points){
      K1 = dnorm((sm[m]-s)/h11)/h11
      D_im1 = c(X[i,],X[i,]*(sm[m]-s)/h11)
      absum1 = absum1 +K1*D_im1*(Y1[i,m]-g1[i])
    }
  }
  ab1 = solve(Wh1)%*%absum1
  
  absum2 = matrix(0,nrow=4,ncol=1)
  for(i in 1:n){
    for(m in 1:max.points){
      K2 = dnorm((sm[m]-s)/h21)/h21
      D_im2 = c(X[i,],X[i,]*(sm[m]-s)/h21)
      absum2 = absum2 +K2*D_im2*(Y2[i,m]-g2[i])
    }
  }
  ab2 = solve(Wh2)%*%absum2
  
  a = cbind(ab1[1:2,],ab2[1:2,])
  
  return(a)
  
}

FVCSIM_est_alpha<-function(Y,X,z,bns,beta){
  s = seq(from=0,to=1,length.out=max.points)
  alpha=matrix(0,nrow=max.points,ncol=4)
  for(t in 1:max.points){
    alpha[t,] = FVCSIM_est_alphat(Y,X,z,s[t],bns,beta)
  }
  return(alpha)
}

FVCSIM_est_ab<-function(Y,X,z,z_est,ans,beta,verbose=F){
  
  n1 = nrow(Y)/2
  Y1=Y[1:n1,]
  Y2=Y[-(1:n1),]
  
  alpha1s.seq = ans[,1:2]
  alpha2s.seq = ans[,3:4]
  Y1t = Y1- X%*%t(alpha1s.seq)
  Y2t = Y2- X%*%t(alpha2s.seq)
  Yt=rbind(Y1t,Y2t)
  
  beta1 = beta[,1]
  beta2 = beta[,2]
  
  diff1 = Y1t[,1]-z%*%beta1
  diff2 = Y2t[,1]-z%*%beta2
  
  h11=rule.thumb_FVCSIM(diff1)
  h21=rule.thumb_FVCSIM(diff2)
  
  # h11=dpill(Y1t,z%*%beta1)
  # h21=dpill(Y2t,z%*%beta2)
  z_est_rep = matrix(rep(z_est,n),nrow=n,byrow=TRUE)
  mx1_z=cbind(rep(1,n1),(z-z_est_rep)%*%beta1/h11)
  mx2_z=cbind(rep(1,n1),(z-z_est_rep)%*%beta2/h21)
  k1_z=dnorm((z-z_est_rep)%*%beta1/h11)/h11
  k2_z=dnorm((z-z_est_rep)%*%beta2/h21)/h21
  
  ##eati a1,b1
  
  sig_sum1=matrix(0,ncol=2,nrow=2)
  for(i in 1:n1)
  {
    sig_sum1=sig_sum1+max.points*(k1_z[i])*(mx1_z[i,])%*%t(mx1_z[i,])
  }
  sig_inv1=solve(sig_sum1)
  
  re1=c(0,0)
  for(i in 1:n1){
    for(j in 1:max.points){
      re1=re1+k1_z[i]*(mx1_z[i,])*Y1t[i,j]
    }
  }
  res1=sig_inv1%*%re1
  
  ##esti a2,b2
  
  sig_sum2=matrix(0,ncol=2,nrow=2)
  for(i in 1:n1){
    sig_sum2=sig_sum2+max.points*(k2_z[i])*(mx2_z[i,])%*%t(mx2_z[i,])
  }
  sig_inv2 <- solve(sig_sum2)
  
  re2=c(0,0)
  for(i in 1:n1){
    for(j in 1:max.points){
      re2=re2+k2_z[i]*(mx2_z[i,])*Y2t[i,j]
    }
  }
  res2=sig_inv2%*%re2
  
  rest=cbind(res1,res2)
  
  return(c(rest[1,],rest[2,]))
}

FVCSIM_est_abt<-function(Y,X,z,ans,beta){
  n=nrow(Y)/2
  hh=matrix(0,nrow=n,ncol=4)
  for(i in 1:n){
    hh[i,]=FVCSIM_est_ab(Y,X,z,z[i,],ans,beta)
  }
  hh_a=hh[,1:2]
  hh_b=hh[,3:4]
  return(hh)
}

FVCSIM_est_beta<-function(Y,X,z,ans,bns,beta){
  n1 = nrow(Y)/2
  Y1=Y[1:n1,]
  Y2=Y[-(1:n1),]
  
  alpha1s.seq = ans[,1:2]
  alpha2s.seq = ans[,3:4]
  Y1t = Y1- X%*%t(alpha1s.seq)
  Y2t = Y2- X%*%t(alpha2s.seq)
  Yt=rbind(Y1t,Y2t)
  
  beta1=beta[,1]
  beta2=beta[,2]
  
  diff1 = Y1t[,1]-z%*%beta1
  diff2 = Y2t[,1]-z%*%beta2
  
  h11=rule.thumb_FVCSIM(diff1)
  h21=rule.thumb_FVCSIM(diff2)
  # h11=dpill(Y1t,z%*%beta1)
  # h21=dpill(Y2t,z%*%beta2)
  
  est_a=bns[,1:2]
  est_b=bns[,3:4]
  
  ##esti beta1
  Omega_1=matrix(0,nrow=3,ncol=3)
  
  for (it in 1:n1){
    sum_1 = 0
    for (i in 1:n1){
      delta_z = z[i,]-z[it,]
      Kt = as.vector(dnorm(delta_z%*%beta1/h11)/h11)
      sum_1 = sum_1+Kt
    }
  }
  
  for(i in 1:n1){
    for(j in 1:n1){
      delta_z=z[i,]-z[j,]
      k1=as.vector(dnorm(delta_z%*%beta1/h11)/h11)/(sum_1/n1)
      Omega_1=Omega_1+max.points*(est_b[i,1])^2*(k1)*delta_z%*%t(delta_z)/(h11^2)
    }
  }
  Omega_inv_1=solve(Omega_1)
  
  beta_s_1=c(0,0,0)
  
  for(i in 1:n1){
    for(j in 1:n1){
      for(t in 1:max.points){
        delta_z=z[i,]-z[j,]
        k1=as.vector(dnorm(delta_z%*%beta1/h11)/h11)/(sum_1/n1)
        beta_s_1=beta_s_1+est_b[i,1]*(delta_z/h11)*k1*(Y1t[i,t]-est_a[i,1])
      }
    }
  }
  beta_1=Omega_inv_1%*%beta_s_1
  
  ##esti beta2
  
  for (it in 1:n1){
    sum_2 = 0
    for (i in 1:n1){
      delta_z = z[i,]-z[it,]
      Kt = as.vector(dnorm(delta_z%*%beta2/h21)/h21)
      sum_2 = sum_2+Kt
    }
  }
  
  Omega_2=matrix(0,nrow=3,ncol=3)
  for(i in 1:n1){
    for(j in 1:n1){
      delta_z=z[i,]-z[j,]
      k2=as.vector(dnorm(delta_z%*%beta2/h21)/h21)/(sum_2/n1)
      Omega_2=Omega_2+max.points*(est_b[i,2])^2*(k2)*delta_z%*%t(delta_z)/(h21^2)
    }
  }
  Omega_inv_2=solve(Omega_2)
  
  beta_s_2=c(0,0,0)
  
  for(i in 1:n1){
    for(j in 1:n1){
      for(t in 1:max.points){
        delta_z=z[i,]-z[j,]
        k2=as.vector(dnorm(delta_z%*%beta2/h21)/h21)/(sum_2/n1)
        beta_s_2=beta_s_2+est_b[i,2]*(delta_z/h21)*k2*(Y2t[i,t]-est_a[i,2])
      }
    }
  }
  beta_2=Omega_inv_2%*%beta_s_2
  
  
  beta_11=as.vector(beta_1/as.vector(sqrt(t(beta_1)%*%(beta_1))))
  beta_21=as.vector(beta_2/as.vector(sqrt(t(beta_2)%*%(beta_2))))
  return(cbind(beta_11,beta_21))
  # return(sum_1)
  # list(sum_1=sum_1,Y1t=Y1t,h11=h11,diff1=diff1)
  
}

FVCSIM_iteration<-function(max_iter,tol,Y,X,z,g0,beta0){
  
  ans = FVCSIM_est_alpha(Y,X,z=z,g0,beta=beta0)
  bns = FVCSIM_est_abt(Y,X,z=z,ans=ans,beta=beta0)
  cns = FVCSIM_est_beta(Y,X,z=z,ans=ans,bns=bns,beta=beta0)
  bnst = FVCSIM_est_abt(Y,X,z=z,ans=ans,beta=cns)
  
  beta = cns
  iter = 0
  
  while(iter<max_iter){
    
    beta_old = beta
    anst = FVCSIM_est_alpha(Y,X,z=z,bns=bnst,beta=beta_old)
    bnstt = FVCSIM_est_abt(Y,X,z=z,ans=anst,beta=beta_old)
    cnstt = FVCSIM_est_beta(Y,X,z=z,ans=anst,bns=bnstt,beta=beta_old)
    beta = cnstt
    
    if (sum((beta - beta_old)^2,na.rm=TRUE) < tol) {
      break
    }
    
    iter = iter+1
    
  }
  list(alpha=anst, g = bnstt ,beta = cnstt)
  
}

FVCSIM_est_cov<-function(Y,X,z,t,t_est,ans,bns,beta){
  
  n=nrow(Y)/2
  Y1=Y[1:n,]
  Y2=Y[-(1:n),]
  alpha1s.seq = ans[,1:2]
  alpha2s.seq = ans[,3:4]
  Y1t = Y1- X%*%t(alpha1s.seq)
  Y2t = Y2- X%*%t(alpha2s.seq)
  Yt=rbind(Y1t,Y2t)
  
  est_a=bns[,1:2]
  est_b=bns[,3:4]
  
  Yi_1=Y1t-est_a[,1]
  Yi_2=Y2t-est_a[,2]
  
  beta1=beta[,1]
  beta2=beta[,2]
  
  diff1 = Y1t[,1]-z%*%beta1
  diff2 = Y2t[,1]-z%*%beta2
  
  h12=rule.thumb_FVCSIM(diff1)
  h22=rule.thumb_FVCSIM(diff2)
  # h12=dpill(Y1t,z%*%beta1)
  # h22=dpill(Y2t,z%*%beta2)
  # 
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

FVCSIM_est_eta<-function(Y,X,z,t,ans,bns,beta){
  
  n=nrow(Y)/2
  eta=array(0,dim=c(n,2,max.points))
  for(s in 1:max.points){
    ts=seq(0.01,0.99,length=max.points)
    eta[,,s]=FVCSIM_est_cov(Y,X,z,t,ts[s],ans,bns,beta)
  }
  return(eta)
}

FVCSIM_est_covt<-function(Y,X,z,t,ans,bns,beta){
  
  n = nrow(Y)/2
  eta=array(0,dim=c(n,max.points,2))
  for(s in 1:max.points){
    ts=seq(0.01,0.99,length=max.points)
    eta[,s,]=FVCSIM_est_cov(Y,X,z,t,ts[s],ans,bns,beta)
  }
  
  c11=matrix(0,max.points,max.points)
  c12=matrix(0,max.points,max.points)
  c21=matrix(0,max.points,max.points)
  c22=matrix(0,max.points,max.points)
  
  for(s in 1:max.points){
    for(t in 1:max.points){
      suu = matrix(0,2,2)
      for(i in 1:n){
        suu=suu+eta[i,s,]%*%t(eta[i,t,])
      }
      Rhat=suu/n
      c11[s,t]=Rhat[1,1]
      c12[s,t]=Rhat[1,2]
      c21[s,t]=Rhat[2,1]
      c22[s,t]=Rhat[2,2]
    }
  }
  cov1122=list(cov11=c11,cov22=c22)
  return(cov1122)
}
