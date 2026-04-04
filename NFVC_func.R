est_rho_beta <- function(t_est, Y, X, Tpoints, W, Ytrue = NULL, verbose = F)
{
  # dimensions
  N = nrow(Y)
  
  # kernel matrices
  T_diff = Tpoints - t_est
  h = rule.thumb(T_diff[1,])/3
  Ker_mat = gaussian.kernel(T_diff, h)
  ns = rowSums(!is.na(Tpoints))
  fs = rowSums(Ker_mat, na.rm = T)/ns # density estimation
  #fs = sum(Ker_mat, na.rm = T)/sum(ns)
  
  Y_Ker_mat = Y*Ker_mat/ns/fs
  Y2_Ker_mat = Y^2*Ker_mat/ns/fs
  Y_Ker_vec = rowSums(Y_Ker_mat, na.rm = T) # Y_Ker_vec is kernel approximation to Y(t)
  Y2_Ker_vec = rowSums(Y2_Ker_mat, na.rm = T) # Y2_Ker_vec is kernel approximation to Y^2
  
  zeta = apply(Y*Ker_mat, 1, function(s){return(mean(s, na.rm = T))})
  zetat = t((zeta/fs)%*%t(matrix(1, nrow = N, ncol = 1)))
  
  # nu_mat is kernel approximation to Y%*%t(Y)
  nu_mat = Y_Ker_vec%*%t(Y_Ker_vec)
  #diag(nu_mat) = Y2_Ker_vec
  
  # some basic matrices
  I_N = Diagonal(n = N, x = 1)
  tWW = t(W)%*%W
  d_tWW = diag(tWW)
  tW_W = t(W) + W
  
  rho = 0.5
  del = 1
  iter = 0
  while (max(abs(del))>10^{-4} &  iter < 20)
  {
    #### Update beta
    St = diag(N)-rho*W
    Mt = t(St)%*%St
    Dt = diag(diag(Mt)^{-1})
    Ut = Dt%*%Mt ### Ui is the i-th column of Ut
    A = Dt%*%t(St)%*%X
    beta = solve(t(A)%*%A)%*%(t(A)%*%rowSums(Ut*zetat))   
    
    if (verbose)
      cat(del, "\n")
    
    #### Update rho
    # S matrix and related derivatives
    S = I_N - rho*W
    tSS = t(S)%*%S
    tSS_rho = -tW_W + 2*rho*tWW # 1st order derivative of tSS
    tSS_rho2 = 2*tWW # 2nd order derivative of tSS
    
    # D matrix and related derivatives
    D = 1/diag(tSS)
    D1 = -2*rho*D^2*diag(tWW)
    D2 = -2*D^2*d_tWW + 8*rho^2*D^3*d_tWW^2
    
    # U_mat and u_mat u1_mat
    U_mat = D*tSS
    U1_mat = D1*tSS - D*tW_W + 2*rho*D*tWW
    u_mat = t(U1_mat)%*%U_mat
    U2_mat = D2*tSS + 2*D1*tSS_rho+D*tSS_rho2
    u1_mat = t(U2_mat)%*%U_mat + t(U1_mat)%*%U1_mat
    
    # xi vectors and v_vec
    Xbeta = X%*%beta
    xi = D*t(S)%*%Xbeta
    xi_rho = (D1*t(S) - D*t(W))%*%Xbeta
    xi_rho2 = (D2*t(S) - 2*D1*t(W))%*%Xbeta
    xi_beta = D*t(S)%*%X
    v_vec_rho = t(U1_mat)%*%xi + t(U_mat)%*%xi_rho
    v1_vec_rho = t(U2_mat)%*%xi + 2*t(U1_mat)%*%xi_rho + t(U_mat)%*%xi_rho2
    
    # first order derivatives
    # nu_mat = Ytrue%*%t(Ytrue)
    # Y_Ker_vec = Ytrue
    G1_rho = 2/N*sum(nu_mat*u_mat)
    G2_rho = -2/N*sum(Y_Ker_vec*v_vec_rho)
    G3_rho = 2/N*sum(xi_rho*xi)
    G_rho = G1_rho + G2_rho + G3_rho
    
    # # use F matrix
    # F_vec = D*t(S)%*%(S%*%Ytrue - Xbeta)
    # F_rho = (D1*t(S) - D*t(W))%*%(S%*%Ytrue - Xbeta)-D*t(S)%*%W%*%Ytrue
    # mean(F_vec*F_rho)*2
    
    # second order derivative
    
    H1_rho = 2/N*sum(nu_mat*u1_mat)
    H2_rho = -2/N*sum(Y_Ker_vec*v1_vec_rho)
    H3_rho = 2/N*sum(xi_rho*xi_rho + xi_rho2*xi)
    H_rho = H1_rho + H2_rho + H3_rho
    del = H_rho^{-1}*G_rho
    rho = rho - del
    iter = iter + 1
    if(rho<0 | rho>1) {rho = runif(1)}
    
  }
  
  return(c(rho, beta))
  
}


est_entire_timeline <- function(Y, X, Tpoints, W, num.points)
{
  tlsub = seq(from = 0, to = 1, length.out = num.points)
  re <- matrix(0, nrow = num.points, ncol = 3)
  for (ii in 1:num.points) {
    #cat(ii, "\r")
    re[ii,] = est_rho_beta(tlsub[ii], Y, X, Tpoints, W)
  }
  list(est_coef = re, timeline = tlsub)
}

obtain.kernel = function(Y, Tpoints, t)
{
  n = dim(Tpoints)[1]
  max.points = dim(Tpoints)[2]
  Ys = Tpoints[1,] 
  h = 1*density(Ys[Ys<999])$bw
  ### Store kernel
  K = matrix(999, n, max.points)
  for (i in 1:n) {
    ti = Tpoints[i,]
    total.points = sum(ti<999)
    K[i, 1:total.points] = gaussian.kernel.mat(t-ti[1:total.points], h)
  }
  list(K = K, h = h)
}

### Define kernel functions
gaussian.kernel.mat = function(u, h)
{
  ### Input t is a vector
  u = u/h
  return(exp(-0.5*u^2)/sqrt(2*pi))
}


rule.thumb = function(x)
{
  1.06*sd(x, na.rm = T)*length(x)^{-1/5}
}

gaussian.kernel = function(u, h)
{
  ### Input t is a vector
  u = u/h
  return(exp(-0.5*u^2)/sqrt(2*pi))
}
