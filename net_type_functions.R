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




plot_g=function(Y,X,z,W,d_alpha,beta){
  
  n = nrow(Y)/2
  Y1 = Y[1:n,]
  Y2 = Y[-(1:n),]
  p = ncol(X)
  
  alpha1s.seq=d_alpha[,5:(p+4)]
  alpha2s.seq=d_alpha[,(p+5):(2*p+4)]
  
  d11s.mat = matrix(rep(d_alpha[,1],each=n),nrow=n,ncol=max.points)
  d12s.mat = matrix(rep(d_alpha[,3],each=n),nrow=n,ncol=max.points)
  d21s.mat = matrix(rep(d_alpha[,2],each=n),nrow=n,ncol=max.points)
  d22s.mat = matrix(rep(d_alpha[,4],each=n),nrow=n,ncol=max.points)
  
  beta1=beta[,1]
  beta2=beta[,2]
  
  Y1t = Y1-X%*%t(alpha1s.seq)-d11s.mat*W%*%Y1-d21s.mat*W%*%Y2
  Y2t = Y2-X%*%t(alpha2s.seq)-d22s.mat*W%*%Y2-d12s.mat*W%*%Y1
  
  rest = matrix(0,nrow=20,ncol=2)
  zz = seq(-1.5,1.5,length=20)
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

plot_g_real=function(Y,X,z,W,d_alpha,beta){
  
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
  zz1 = seq(min(z%*%beta1),max(z%*%beta1),length=20)
  zz2 = seq(min(z%*%beta2),max(z%*%beta2),length=20)
  for(k in 1:20){
    zzb1=zz1[k]
    zzb2=zz2[k]
    zb1 = z%*%beta1
    zb2 = z%*%beta2
    # h11=dpill(zb1,Y1t)
    # h21=dpill(zb2,Y2t)
    h11 = rule.thumb(Y1t)
    h21 = rule.thumb(Y2t)
    mx1_z=cbind(rep(1,n),(zb1-zzb1)/h11)
    mx2_z=cbind(rep(1,n),(zb2-zzb2)/h21)
    k1_z=dnorm((zb1-zzb1)/h11)/h11
    k2_z=dnorm((zb2-zzb2)/h21)/h21
    
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


rule.thumb = function(x){
  0.9*min(quantile(x,0.75,na.rm=T)-quantile(x,0.25,na.rm=T),sd(x,na.rm=T))*length(x)^(-1/5)
}

gaussian.kernel = function(u,h){
  u = u/h
  return(dnorm(u))
}

get.YK <- function(Y,Tpoints)
{
  # dimensions
  N = nrow(Y)/2
  
  #Y---2n*maxpoints
  Y1s=Y[1:N,]
  Y2s=Y[-(1:N),]
  
  
  tlsub = seq(0,1,length=max.points)
  re = matrix(0,nrow=2*n,ncol=max.points)
  for(ii in 1:max.points){
    
    t_est = tlsub[ii]
    T_diff = Tpoints - t_est#N*numpointsά
    h = rule.thumb(T_diff[1,])/3
    Ker_mat = gaussian.kernel(T_diff, h)
    
    T_diff.hat = rbind(T_diff,T_diff)
    Ker_mat.hat = gaussian.kernel(T_diff.hat, h)
    
    
    ns = rowSums(!is.na(Tpoints))
    
    ns.hat = rowSums(!is.na(rbind(Tpoints,Tpoints)))
    
    fs = rowSums(Ker_mat, na.rm = T)/ns 
    
    fs.hat = rowSums(Ker_mat.hat, na.rm = T)/ns.hat
    
    Y_Ker_mat_1 = Y1s*Ker_mat/ns/fs
    Y2_Ker_mat_1 = Y1s^2*Ker_mat/ns/fs
    Y_Ker_vec_1 = rowSums(Y_Ker_mat_1, na.rm = T) 
    Y2_Ker_vec_1 = rowSums(Y2_Ker_mat_1, na.rm = T) 
    
    Y_Ker_mat_2 = Y2s*Ker_mat/ns/fs
    Y2_Ker_mat_2 = Y2s^2*Ker_mat/ns/fs
    Y_Ker_vec_2 = rowSums(Y_Ker_mat_2, na.rm = T) 
    Y2_Ker_vec_2 = rowSums(Y2_Ker_mat_2, na.rm = T) 
    
    Y_Ker_mat.hat = rbind(Y_Ker_mat_1,Y_Ker_mat_2)
    Y2_Ker_mat.hat = rbind(Y2_Ker_mat_1,Y2_Ker_mat_2)
    #Y_Ker_mat = Y_Ker_mat.hat/ns/fs
    #Y2_Ker_mat = Y2_Ker_mat.hat/ns/fs
    Y_Ker_vec = rowSums(Y_Ker_mat.hat, na.rm = T)
    Y2_Ker_vec = rowSums(Y2_Ker_mat.hat, na.rm = T)
    
    zeta = apply(Y*Ker_mat.hat, 1, function(s){return(mean(s, na.rm = T))})
    zetat = t((zeta/fs.hat)%*%t(matrix(1, nrow = 2*N, ncol = 1)))
    
    Ys = zeta/fs.hat
    
    # Ym = matrix(Ys,ncol=2)
    
    re[,ii] = Ys
  }
  
  return(re)
}


create_W3 <- function(nrow,ncol, torus = FALSE, style = "W") {
  n <- nrow * ncol
  
  # 从 Rook 邻接开始
  nb_rook <- spdep::cell2nb(nrow, ncol, type = "rook", torus = torus)
  
  # 过滤掉上下连接，只保留左右邻居
  nb_left_right <- nb_rook
  for(i in 1:length(nb_left_right)){
    neighbors <- nb_left_right[[i]]
    if(length(neighbors) > 0){
      row_i <- ceiling(i / ncol)
      keep_neighbors <- sapply(neighbors, function(j) {
        row_j <- ceiling(j / ncol)
        row_i == row_j  # 在同一行
      })
      nb_left_right[[i]] <- neighbors[keep_neighbors]
    }
  }
  
  W3 <- spdep::nb2mat(nb_left_right, style = style, zero.policy = TRUE)
  rownames(W3) <- 1:n
  colnames(W3) <- 1:n
  
  return(W3)
}


create_W5 <- function(nrow,ncol, torus = FALSE, style = "W") {
  n <- nrow * ncol
  
  # 直接创建 Queen 邻接
  nb_queen <- spdep::cell2nb(nrow, ncol, type = "queen", torus = torus)
  W5 <- spdep::nb2mat(nb_queen, style = style, zero.policy = FALSE)
  rownames(W5) <- 1:n
  colnames(W5) <- 1:n
  
  return(W5)
}


