ifelse<-function(x, val1, val0)
{
  if (x)
  {
    return(val1)
  }else
  {
    return(val0)
  }
}

tr<-function(M)
{
  return(sum(diag(M)))
}


getPowerLawW<-function(N, alpha, normalize = T)           
{
  Nfollowers = rpldis(N, 1, alpha)
  A = sapply(Nfollowers, function(n) {              
    vec = rep(0, N) 
    vec[sample(1:N, min(n,N))] = 1 
    return(vec)
  })
  diag(A) = 0 
  ind = which(rowSums(A)==0)                      
  for (i in ind)
  {
    A[i, sample(setdiff(1:N,i), 3)] = 1             
  }
  if (!normalize)
    return(A)
  W = A/rowSums(A)
  return(W)
}

getBlockW<-function(N, Nblock, normalize = T)                                                 
{
  if (N%%Nblock==0){                                                                 
    isDiagList = rep(list(matrix(1, nrow = N/Nblock, ncol = N/Nblock)), Nblock)      
    mList = rep(list(matrix(rbinom((N/Nblock)^2, size = 1, prob = 0.9*N^{-1}),   
                            nrow = N/Nblock, ncol = N/Nblock)), Nblock)
  }
  else
  {
    isDiagList = rep(list(matrix(1, nrow = floor(N/Nblock), ncol = floor(N/Nblock))), Nblock-1) 
    isDiagList[[length(Nblock)]] = matrix(1, nrow = N%%Nblock, ncol = N%%Nblock)
    
    mList = rep(list(matrix(rbinom(floor(N/Nblock)^2, size = 1, prob = 0.9*N^{-1}),     
                            nrow = floor(N/Nblock), ncol = floor(N/Nblock))), Nblock-1)
    mList[[Nblock]] = matrix(rbinom(floor(N/Nblock)^2, size = 1, prob = 0.9*N^{-1}),   
                             nrow = floor(N/Nblock), ncol = floor(N/Nblock))
  }
  isDiag = bdiag(isDiagList)   
  offDiag = which(isDiag == 0, arr.ind = T)                                         
  mList = lapply(mList, function(M){
    ind = which(rowSums(M)==0)
    if (length(ind)>0)
      M[cbind(ind, sample(1:nrow(M), length(ind)))] = 1
    return(M)
  })
  bA = bdiag(mList)
  bA[offDiag] = rbinom(nrow(offDiag), size = 1, prob = 0.3/N)                        
  bA = as.matrix(bA)
  upperInd = which(upper.tri(bA), arr.ind = T)
  
  bA[upperInd[,2:1]] = bA[upper.tri(bA)]
  diag(bA) = 0
  
  ind = which(rowSums(bA)==0)                                 
  for (i in ind)
  {
    bA[i, sample(setdiff(1:N,i), 3)] = 1                      
  }
  
  if (!normalize)
    return(bA)
  W = bA/rowSums(bA)                                        
  return(W)
}

getDyadW<-function(N, N1 = 10, delta = 1.2, normalize = T)                                             ### simulate Dyad network W; N1: number of mutual pairs 2N1*N, delta: P((0,1)) = P((1,0)) = 0.5*N^{delta}
{
  A = matrix(0, nrow = N, ncol = N)                                                                    ### use A to store network structure
  
  ######################################### mutual follower ###########################################
  ind = which(upper.tri(A), arr.ind = T)                                                               ### record all the index of upper triangular matrix
  indM = ind[sample(1:nrow(ind), N*N1),]                                                               ### sample N*N1 as mutual pairs in the upper triangular matrix
  A[indM] = 1                                                                                          ### the corresponding links are set to be 1
  A[indM[,2:1]] = 1                                                                                    ### the following matrix is set to be symmetric, as a result, mutual pairs are 2N1*N
  
  ######################################### single relationship #######################################
  ind1 = which(A==0&upper.tri(A), arr.ind = T)                                                         ### record all the zero pairs in the upper triangular matrix
  indS = ind1[sample(1:nrow(ind1), N^delta),]                                                          ### choose N^delta index as single relations
  tmp = sample(1:nrow(indS), floor(N^delta/2))                                                         ### randomly choose 0.5*N^delta to be in the lower triangular matrix
  indS[tmp,] = indS[tmp, 2:1]                                                                          ### change the corresponding index to be inverse
  A[indS] = 1                                                                                          ### the single following relation is set to be 
  diag(A) = 0                                                                                          ### aii = 0
  if (!normalize)
    return(A)
  W = A/rowSums(A)                                                                                     ### W is row-normalized
  W = as(W, "dgCMatrix")
  return(W)
}

rule.thumb = function(x){
  0.9*min(quantile(x,0.75,na.rm=T)-quantile(x,0.25,na.rm=T),sd(x,na.rm=T))*length(x)^(-1/5)
}

gaussian.kernel = function(u,h){
  u = u/h
  return(dnorm(u))
}

###generate data
get.func.network = function(n, max.points, net_type,Sige1,eps){
  if (net_type=="1")
  {W = getDyadW(n, N1 = 10, delta = 1.2, normalize = T)}
  
  if (net_type=="2")
  {W = getBlockW(n,10)}
  if (net_type=="3")
  {W = getPowerLawW(n, alpha = 2.5, normalize = T) }
  
  
  t = seq(from = 0, to = 1, length.out = max.points) 
  
  alpha11 = matrix(t^2,nrow=1)
  alpha12 = matrix((1-t)^2,nrow=1)
  alpha21 = matrix(5*(t-0.5)^2,nrow=1)
  alpha22 = matrix(t^0.5,nrow=1)
  alpha1.seq = rbind(alpha11,alpha21)
  alpha2.seq = rbind(alpha12,alpha22)
  alpha = rbind(alpha11,alpha21,alpha12,alpha22)
                
  XSigma = 0.5^abs(outer(1:2,1:2,"-"))
  X = mvrnorm(n,mu=rep(0,nrow(XSigma)),Sigma = XSigma)

  ZSigma = diag(3)
  z = mvrnorm(n,mu=rep(0,nrow(ZSigma)),Sigma = ZSigma)
  
  beta1 = c(2,-1,3)/sqrt(14)
  beta2 = c(1,2,-3)/sqrt(14)
  beta=cbind(beta1,beta2)

  g1 = tanh(z%*%beta1)
  g2 = -atan(z%*%beta2)
  g = rbind(g1,g2)
  
  phi11.seq = sqrt(2)*sin(2*pi*t)
  phi12.seq = sqrt(2)*cos(2*pi*t)
  phi21.seq = sqrt(2)*cos(2*pi*t)
  phi22.seq = sqrt(2)*sin(2*pi*t)
  phi1.seq = rbind(phi11.seq,phi12.seq)
  phi2.seq = rbind(phi21.seq,phi22.seq)
  
  bi11 = rnorm(n,0,sqrt(1.2))
  bi12 = rnorm(n,0,sqrt(0.6))
  bi21 = rnorm(n,0,1)
  bi22 = rnorm(n,0,sqrt(0.5))
  bi1.seq = cbind(bi11,bi12)
  bi2.seq = cbind(bi21,bi22)
  
  etai1 = bi1.seq%*%phi1.seq
  etai2 = bi2.seq%*%phi2.seq
  eta = rbind(etai1,etai2)

  ee = eigen(Sige1)
  In = Diagonal(n = n, x = 1)
  eps1 = kronecker(t(sqrt(ee$values)*t(ee$vectors)), In)%*%eps
  
  d11 = matrix(0.2*sin(2*t*pi)+0.3,nrow = 1)
  d12 = matrix(2*(t^2-t)+0.5,nrow=1)
  d21 = matrix(2*(t-t^2),nrow=1)
  d22 = matrix(0.2*cos(2*t*pi)+0.3,nrow=1)
  d.mat = rbind(d11,d21,d12,d22)
  
  
  Ymat = matrix(0,nrow=2*n,ncol=max.points)
  for(m in 1:max.points){
    D = matrix(d.mat[,m],nrow=2)
    A = matrix(alpha[,m],nrow=2)
    DkW = kronecker(t(D),W)
    DkW2 = DkW%*%DkW
    DkW3 = DkW2%*%DkW
    DkW4 = DkW3%*%DkW
    right.mat = as.vector(X%*%A)+g+eps1[,m]
    Ymat[,m] = as.vector(right.mat+DkW%*%right.mat+DkW2%*%right.mat+DkW3%*%right.mat+DkW4%*%right.mat)
  }
  
  Y1 = Ymat[(1:n),]
  Y2 = Ymat[-(1:n),]
  
  YNA = matrix(NA, 2*n, max.points)
  
  Tpoints = matrix(NA, n, max.points)
  
  for (i in 1:n) {
    total.points = sample(min.points:max.points, 1, prob = rep(1:(max.points-min.points+1)/(max.points-min.points+1)))
    
    ind = sort(sample(1:max.points, size = total.points, replace = F))
    
    Tpoints[i, 1:total.points] = s[ind]
    
    YNA[i, 1:total.points] = Ymat[i,ind]
    YNA[n+i, 1:total.points] = Ymat[n+i,ind]
  }
  
  Y1NA = YNA[(1:N),]
  Y2NA = YNA[-(1:N),]
  
  
  list(Y = Ymat, YNA = YNA,Y1 = Y1, Y2 = Y2, X = X, z = z, g = g,beta=beta,W = W, time_seq = t)
}


###estimate d and alpha
est_d_alpha <- function(Y, X, W,D,B,Sige,g)
{
  n = nrow(Y)/2
  tlsub = seq(from = 0, to = 1, length.out = max.points)
  re <- matrix(0, nrow = max.points, ncol = 8)
  for (ii in 1:(max.points)) {
    Ym = matrix(Y[,ii],ncol=2)-g[,1:2]
    result = MSAR_Lse_rcpp(Ym, X, W, D,B,Sige)
    re[ii,] = as.vector(result$theta)
  }
  
 return(re)
}



est_ab<-function(Y,X,z,z_est,W,ans,beta,verbose=F){
  
  n1 = nrow(Y)/2
  Y1=Y[1:n1,]
  Y2=Y[-(1:n1),]
  
  d11s.mat = matrix(rep(ans[,1],each=n1),nrow=n1,ncol=max.points)
  d12s.mat = matrix(rep(ans[,3],each=n1),nrow=n1,ncol=max.points)
  d21s.mat = matrix(rep(ans[,2],each=n1),nrow=n1,ncol=max.points)
  d22s.mat = matrix(rep(ans[,4],each=n1),nrow=n1,ncol=max.points)
  alpha1s.seq = ans[,5:6]
  alpha2s.seq = ans[,7:8]
  Y1t = Y1- X%*%t(alpha1s.seq)-d11s.mat*W%*%Y1-d21s.mat*W%*%Y2
  Y2t = Y2- X%*%t(alpha2s.seq)-d22s.mat*W%*%Y2-d12s.mat*W%*%Y1

  Yt=rbind(Y1t,Y2t)
  
  beta1=beta[,1]
  beta2=beta[,2]
  
  # h11=dpill(z%*%beta1,Y1t)
  # h21=dpill(z%*%beta2,Y2t)
  h11 = rule.thumb(Y1t)
  h21 = rule.thumb(Y2t)
  
  # zbeta1 = as.numeric(z %*% beta1)
  # zbeta2 = as.numeric(z %*% beta2)
  # 
  # h11=dpill(zbeta1,Y1t) #0.99
  # h21=dpill(zbeta2,Y2t) #0.78
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
  
  return(c(rest[1,],rest[2,])) ####(a1,a2,b1,b2)
}

est_abt<-function(Y,X,z,W,ans,beta){
  n = nrow(Y)/2
  hh=matrix(0,nrow=n,ncol=4)
  for(i in 1:n){
    hh[i,]=est_ab(Y,X,z,z[i,],W,ans,beta)
  }
  # hh_a=hh[,1:2]  ##esti_g
  # hh_b=hh[,3:4]
  return(hh)
}


est_beta<-function(Y,X,z,W,ans,bns,beta){
  n1 = nrow(Y)/2
  Y1=Y[1:n1,]
  Y2=Y[-(1:n1),]
  d11s.mat = matrix(rep(ans[,1],each=n1),nrow=n1,ncol=max.points)
  d12s.mat = matrix(rep(ans[,3],each=n1),nrow=n1,ncol=max.points)
  d21s.mat = matrix(rep(ans[,2],each=n1),nrow=n1,ncol=max.points)
  d22s.mat = matrix(rep(ans[,4],each=n1),nrow=n1,ncol=max.points)
  alpha1s.seq = ans[,5:6]
  alpha2s.seq = ans[,7:8]
  Y1t = Y1- X%*%t(alpha1s.seq)-d11s.mat*W%*%Y1-d21s.mat*W%*%Y2#2n*max.pointsά
  Y2t = Y2- X%*%t(alpha2s.seq)-d22s.mat*W%*%Y2-d12s.mat*W%*%Y1#2n*max.pointsά
  Yt=rbind(Y1t,Y2t)
  
  beta1=beta[,1]
  beta2=beta[,2]
  
  # h11=dpill(Y1t,z%*%beta1)
  # h21=dpill(Y2t,z%*%beta2)
  h11 = rule.thumb(Y1t)
  h21 = rule.thumb(Y2t)
  
  est_a=bns[,1:2]
  est_b=bns[,3:4]
  
  ##esti beta1
  Omega_1=matrix(0,nrow=3,ncol=3)
  sum_1 = numeric(n1)
  for(i in 1:n1){
    for(j in 1:n1){
      if (i != j){
        delta_z=z[i,]-z[j,]
        sum_1[i] = sum_1[i]+dnorm((t(delta_z)%*%beta1)/h11)/h11
      }
    }
  }
  for(i in 1:n1){
    for(j in 1:n1){
      if (i != j){
        delta_z=z[i,]-z[j,]
        k1=as.vector(dnorm(delta_z%*%beta1/h11)/h11)/(sum_1[i]/n1)
        Omega_1=Omega_1+max.points*(est_b[j,1])^2*(k1)*delta_z%*%t(delta_z)/(h11^2)
      }
    }
  }
  Omega_inv_1=solve(Omega_1)
  
  beta_s_1=c(0,0,0)
  
  for(i in 1:n1){
    for(j in 1:n1){
      delta_z=z[i,]-z[j,]
      k1=as.vector(dnorm(delta_z%*%beta1/h11)/h11)/(sum_1[i]/(n1-1))
      for(t in 1:max.points){
        beta_s_1=beta_s_1+est_b[j,1]*(delta_z/h11)*k1*(Y1t[i,t]-est_a[j,1])
      }
    }
  }
  beta_1=Omega_inv_1%*%beta_s_1
  
  ##esti beta2
  sum_2 = numeric(n1)
  for(i in 1:n1){
    for(j in 1:n1){
      if(i != j){
        delta_z = z[i,]-z[j,]
        sum_2[i] = sum_2[i] + dnorm((t(delta_z)%*%beta2)/h21)/h21
      }
    }
  }
  
  Omega_2=matrix(0,nrow=3,ncol=3)
  for(i in 1:n1){
    for(j in 1:n1){
      if(i != j){
        delta_z=z[i,]-z[j,]
        k2=as.vector(dnorm(delta_z%*%beta2/h21)/h21)/(sum_2[i]/(n1-1))
        Omega_2=Omega_2+max.points*(est_b[i,2])^2*(k2)*delta_z%*%t(delta_z)/(h21^2)
      }
    }
  }
  Omega_inv_2=solve(Omega_2)
  
  beta_s_2=c(0,0,0)
  
  for(i in 1:n1){
    for(j in 1:n1){
      delta_z=z[i,]-z[j,]
      k2=as.vector(dnorm(delta_z%*%beta2/h21)/h21)/(sum_2[i]/(n1-1))
      for(t in 1:max.points){
        beta_s_2=beta_s_2+est_b[j,2]*(delta_z/h21)*k2*(Y2t[i,t]-est_a[j,2])
      }
    }
  }
  beta_2=Omega_inv_2%*%beta_s_2
  
  
  beta_11=as.vector(beta_1/as.vector(sqrt(t(beta_1)%*%(beta_1))))
  beta_21=as.vector(beta_2/as.vector(sqrt(t(beta_2)%*%(beta_2))))
  return(cbind(beta_11,beta_21))
  
}


###iteration
iteration<-function(max_iter,tol,Y,X,z,W,D,B,Sige,g,beta,verbose=F){
  
  ans = est_d_alpha(Y,X,W,D,B,Sige,g)
  bns = est_abt(Y,X,z,W,ans=ans,beta)
  cns = est_beta(Y,X,z,W,ans=ans,bns=bns,beta)
  bnst = est_abt(Y,X,z,W,ans,beta=cns)
  anst = est_d_alpha(Y,X,W,D,B,Sige=Sige,g=bnst)
  
  beta = cns
  iter = 0
  
  while(iter<max_iter){
    
    if(verbose)
      cat("","\n")
    
    beta_old = beta
    bnstt = est_abt(Y,X,z,W,ans=anst,beta=beta_old)
    cnstt = est_beta(Y,X,z,W,ans=anst,bns=bnstt,beta=beta_old)
    # beta = 0.6*beta_old+0.4*cnstt
    beta=cnstt
    
    # if (sum((beta - beta_old)^2,na.rm=TRUE) < tol) {
    #   break
    # }
    # 
    
    if (sum((beta[,1] - beta_old[,1])^2) < tol && sum((beta[,2] - beta_old[,2])^2) < tol) break

    iter = iter+1
    
  }
  
  list(d_alpha = anst, g = bnstt ,beta = beta)
  
}


plot_g=function(Y,X,z,W,d_alpha,beta){
  
  n = nrow(Y)/2
  Y1 = Y[1:n,]
  Y2 = Y[-(1:n),]
  
  alpha1s.seq=d_alpha[,5:6]
  alpha2s.seq=d_alpha[,7:8]
  
  d11s.mat = matrix(rep(d_alpha[,1],each=n),nrow=n,ncol=max.points)
  d12s.mat = matrix(rep(d_alpha[,3],each=n),nrow=n,ncol=max.points)
  d21s.mat = matrix(rep(d_alpha[,2],each=n),nrow=n,ncol=max.points)
  d22s.mat = matrix(rep(d_alpha[,4],each=n),nrow=n,ncol=max.points)
  
  beta1=beta[,1]
  beta2=beta[,2]
  
  Y1t = Y1-X%*%t(alpha1s.seq)-d11s.mat*W%*%Y1-d21s.mat*W%*%Y2
  Y2t = Y2-X%*%t(alpha2s.seq)-d22s.mat*W%*%Y2-d12s.mat*W%*%Y1
  
  rest = matrix(0,nrow=20,ncol=2)
  zz = seq(-2,2,length=20)
  for(k in 1:20){
    zzb=zz[k]
    zb1 = z%*%beta1
    zb2 = z%*%beta2
    # h11=dpill(zb1,Y1t)
    # h21=dpill(zb2,Y2t)
    h11 = rule.thumb(Y1t)
    h21 = rule.thumb(Y2t)
    mx1_z=cbind(rep(1,n),(zb1-zzb)/h11)
    mx2_z=cbind(rep(1,n),(zb2-zzb)/h21)
    k1_z=dnorm((zb1-zzb)/h11)/h11
    k2_z=dnorm((zb2-zzb)/h21)/h21
    
    sig_sum1=matrix(0,ncol=2,nrow=2)
    for(i in 1:n){
      sig_sum1=sig_sum1+max.points*(k1_z[i])*(mx1_z[i,])%*%t(mx1_z[i,])
    }
    sig_inv1=solve(sig_sum1)
    
    re1=c(0,0)
    for(i in 1:n){
      for(j in 1:max.points){
        re1=re1+k1_z[i]*(mx1_z[i,])*Y1t[i,j]
      }
    }
    res1=sig_inv1%*%re1
    
    
    sig_sum2=matrix(0,ncol=2,nrow=2)
    for(i in 1:n){
      sig_sum2=sig_sum2+max.points*(k2_z[i])*(mx2_z[i,])%*%t(mx2_z[i,])
    }
    sig_inv2 <- solve(sig_sum2)
    
    re2=c(0,0)
    for(i in 1:n){
      for(j in 1:max.points){
        re2=re2+k2_z[i]*(mx2_z[i,])*Y2t[i,j]
      }
    }
    res2=sig_inv2%*%re2
    
    rest[k,]=cbind(res1[1],res2[1])
    
  }
  list(p_g1=rest[,1],p_g2=rest[,2])
}


#####esti_cov

est_et<-function(Y,X,z,t,t_est,W,ans,bns,beta){
  
  n = nrow(Y)/2
  Y1=Y[1:n,]
  Y2=Y[-(1:n),]
  d11s.mat = matrix(rep(ans[,1],each=n),nrow=n,ncol=max.points)
  d12s.mat = matrix(rep(ans[,3],each=n),nrow=n,ncol=max.points)
  d21s.mat = matrix(rep(ans[,2],each=n),nrow=n,ncol=max.points)
  d22s.mat = matrix(rep(ans[,4],each=n),nrow=n,ncol=max.points)
  alpha1s.seq = ans[,5:6]
  alpha2s.seq = ans[,7:8]
  Y1t = Y1- X%*%t(alpha1s.seq)-d11s.mat*W%*%Y1-d21s.mat*W%*%Y2
  Y2t = Y2- X%*%t(alpha2s.seq)-d22s.mat*W%*%Y2-d12s.mat*W%*%Y1
  Yt=rbind(Y1t,Y2t)
  
  est_a=bns[,1:2]
  est_b=bns[,3:4]
  
  Yi_1=Y1t-est_a[,1]
  Yi_2=Y2t-est_a[,2]
  
  beta1=beta[,1]
  beta2=beta[,2]
  
  h12=dpill(z%*%beta1,Y1t)
  h22=dpill(z%*%beta2,Y2t)/3
  # h12 = rule.thumb(Y1t)
  # h22 = rule.thumb(Y2t)
  
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
    eta[,,s]=est_et(Y,X,z,t,ts[s],W,ans,bns,beta)
  }
  return(eta)
}

est_covt<-function(Y,X,z,W,ans,bns,beta){
  
  n = nrow(Y)/2
  eta=array(0,dim=c(n,max.points,2))
  t = seq(0,1,length=max.points)
  for(s in 1:max.points){
    ts=seq(0,1,length=max.points)
    eta[,s,]=est_et(Y,X,z,t,ts[s],W,ans,bns,beta)
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

getCovCST = function(t)
{
  phi11.seq = matrix(sin(2*t*pi) * sqrt(2), ncol = 1)
  phi12.seq = matrix(cos(2*t*pi) * sqrt(2), ncol = 1)
  phi1.seq = rbind(phi11.seq,phi12.seq)
  phi21.seq = matrix(cos(2*t*pi) * sqrt(2), ncol = 1)
  phi22.seq = matrix(sin(2*t*pi) * sqrt(2), ncol = 1)
  phi2.seq = rbind(phi21.seq,phi22.seq)
  
  cov11 = 1.2*phi11.seq%*%t(phi11.seq)+0.6*phi12.seq%*%t(phi12.seq)
  cov22 = phi21.seq%*%t(phi21.seq)+0.5*phi22.seq%*%t(phi22.seq)
  cov = rbind(cov11,cov22)
  
  list(cov11=cov11,cov22=cov22)
}


