library(genscore)
library(mvtnorm)
library(Matrix)
library(MASS)
library(glasso)

get_data_aug = function(n, beta, beta2, burnin, skip){
  
  p = dim(beta)[1]
  X = matrix(0, nrow = n, ncol = p)
  temp = rnorm(p)
  
  # burn in period
  for(t in 1:burnin){
    for(l in 1:p){
      vec_condition = temp[-l]
      vec_beta = beta[,l][-l]
      vec_beta2 = beta2[,l][-l]
      A = 2*t(vec_condition^2)%*%vec_beta + beta2[l,l]
      B = 2*t(vec_condition)%*%vec_beta2 + beta[l,l]
      mu  = -B/(2*A); sig = sqrt(-1/(2*A))
      temp[l] = rnorm(1, mean = mu, sd = sig)
    }
    if(t%%10000==0){
    }
  }
  
  #samples obtained by skipping 500 samples
  k = 1
  samples = 1
  while(samples <= n){
    for(l in 1:p){
      vec_condition = temp[-l]
      vec_beta = beta[,l][-l]
      vec_beta2 = beta2[,l][-l]
      A = 2*t(vec_condition^2)%*%vec_beta + beta2[l,l]
      B = 2*t(vec_condition)%*%vec_beta2 + beta[l,l]
      mu  = -B/(2*A); sig = sqrt(-1/(2*A))
      temp[l] = rnorm(1, mean = mu, sd = sig)
    }
    if(k%%skip == 1){
      X[samples,] <- temp
      samples<-samples+1
    }
    k = k+1
  }
  return(X)
}
get_element_aug = function(data){
  # this function is used to get gamma and g for single normal conditional graph
  if(is.null(dim(data))){data = matrix(data, nrow = 1, ncol = length(data))}
  
  n = dim(data)[1]
  p = dim(data)[2]
  s = 2*(p + choose(p,2))
  
  pos_mat = matrix(1, p,p)
  pos_mat[upper.tri(pos_mat)] = 0
  diag(pos_mat) = 0
  tmat_lower = which(pos_mat !=0, arr.ind = T)
  tmat = tmat_lower[,c(2,1)]
  
  tmp_deri = function(x, j){
    # a tmp function to get the derivatives of the vector(x1,..xp, x1^2.. xp^2,x_jx_k, x^2_jx^2_k), j is the corresponding coordinate
    p = length(x)
    deri_upper = rep(0,p)
    deri_upper[j] = 1
    
    deri_mid = rep(0,p)
    deri_mid[j] = 2*x[j]
    
    tx_left0 = x[tmat[,1]]
    tx_right0 = x[tmat[,2]]
    
    deri_left0 = tx_left0
    deri_left0[which(tmat[,1]==j)] = 1
    deri_left0[which(tmat[,1]!=j)] = 0
    deri_left0 = deri_left0*tx_right0
    
    deri_right0 = tx_right0
    deri_right0[which(tmat[,2]==j)] = 1
    deri_right0[which(tmat[,2]!=j)] = 0
    deri_right0 = deri_right0*tx_left0
    
    deri_lower0 = deri_left0 + deri_right0
    
    tx_left = x[tmat[,1]]^2
    tx_right = x[tmat[,2]]^2
    
    deri_left = tx_left
    deri_left[which(tmat[,1]==j)] = 2*x[j]
    deri_left[which(tmat[,1]!=j)] = 0
    deri_left = deri_left*tx_right
    
    deri_right = tx_right
    deri_right[which(tmat[,2] == j)] = 2*x[j]
    deri_right[which(tmat[,2] != j)] = 0
    deri_right = deri_right*tx_left
    
    deri_lower = deri_left + deri_right
    
    deri_full = c(deri_upper, deri_mid, 2*deri_lower0, 2*deri_lower)
    return(deri_full)
  }
  
  tmp_deri_ii = function(x, deri){
    # a tmp function to get the derivatives of the vector(x1,..xp, x1^2.. xp^2,x_jx_k, x^2_jx^2_k), j is the corresponding coordinate
    # deri is the vector we got in tmp_deri
    p = length(x)
    tmp_sub = deri[1:p]
    j = which(tmp_sub!=0)
    deri_ii = deri/x[j]
    deri_ii[1:p] = rep(0,p)
    deri_ii[(2*p+1):(2*p+choose(p,2))] = 0
    return(deri_ii)
  }
  
  # now calculate the empirical Gamma and G for the graph
  
  Gamma = matrix(0, s, s)
  G = rep(0, s)
  
  for(i in 1:n){
    x_i = data[i,]
    Gamma_i = matrix(0, s, s)
    G_i = rep(0, s)
    for(j in 1:p){
      deri_i = tmp_deri(x_i,j)
      deri_ii = tmp_deri_ii(x_i, deri_i)
      
      Gamma_i = Gamma_i + deri_i %*% t(deri_i)
      G_i = G_i + deri_ii
      print(paste('i = ', i, '; j = ', j, sep = ''))
    }
    Gamma = Gamma + Gamma_i
    G = G + G_i
  }
  
  return(list(Gamma = Gamma/n, G = G/n))
}
ADMM_gen = function(rho, lambda, alpha, threshold, maxit, elem_x, elem_y){
  
  gamma_x = elem_x$Gamma
  gamma_y = elem_y$Gamma
  g_x = elem_x$G
  g_y = elem_y$G
  s = dim(gamma_x)[1]
  
  Inv_mat1 = solve(gamma_x + gamma_y + diag(rho, s))
  Inv_mat2 = solve(gamma_y + diag(rho,s))
  
  theta_x = rep(0, s)
  delta = rep(0, s)
  Z1 = rep(0, s)
  Z2 = rep(0, s)
  U1 = rep(0, s)
  U2 = rep(0, s)
  
  error = 1 + threshold
  i = 0
  while((error > threshold)&(i<=maxit)){
    
    theta_0 = c(theta_x, delta)
    theta_x = Inv_mat1%*%(gamma_y%*%delta - g_x - g_y + rho*(Z1 - U1))
    delta = Inv_mat2%*%(gamma_y%*%theta_x + g_y + rho*(Z2 - U2))
    
    A1 = theta_x + U1
    A2 = delta + U2
    
    Z1 = sign(A1)*pmax(0, abs(A1)-lambda/rho)
    Z1[1:(2*p)] = A1[1:(2*p)]
    Z2 = sign(A2)*pmax(0, abs(A2)-lambda*alpha/rho)
    
    U1 = U1 + (theta_x - Z1)
    U2 = U2 + (delta - Z2)
    
    theta_1 = c(theta_x, delta)
    dif = theta_1 - theta_0
    error = norm(dif, type = '2')/norm(theta_0, type = '2')
    i = i+1
    print(paste(lambda, '+',alpha, '+',i, '+', error, sep = ' '))
  }
  return(list(theta_x = Z1, delta = Z2))
}
ADMM_single  = function(rho, lambda, threshold, maxit, elem){
  gamma = elem$Gamma
  g = elem$G
  
  s = dim(gamma)[1]
  Inv_mat = solve(gamma + diag(rho,s))
  
  theta = rep(0, s)
  Z1 = rep(0, s)
  U1 = rep(0, s)
  error = 1 + threshold
  while((error > threshold)&(i<=maxit)){
    theta_0 = theta
    
    theta = Inv_mat %*%(-g + rho*Z1 - U1)
    
    A1 = theta + U1/rho
    Z1 = sign(A1)*pmax(0, abs(A1)-lambda/rho)
    Z1[1:(2*p)] = A1[1:(2*p)]
    
    U1 = U1 + rho*(theta - Z1)
    
    theta_1 = theta
    dif = theta_1 - theta_0
    error = norm(dif, type = '2')/norm(theta_0, type = '2')
    i = i+1
    print(paste(lambda, '+',alpha, '+',i, '+', error, sep = ' '))
  }
  return(Z1)
  
  
}
ROC_score = function(Est, SuppDelta, p){
  # Para:
  #   Est: The object returned by JGL function
  #   SuppDelta: Support Matrix for Delta
  #   p: Number of nodes
  # Retunrs:
  #   A vector of c(FPR, TPR)
  Edgeshat <- matrix(0, nrow=p, ncol=p)
  for (i in 1:p){
    for (j in 1:p){
      diff <- abs(Est[i,j])
      if (diff>0){
        Edgeshat[i,j] <- 1
      }
    }
  }
  
  v.estimate <- c()
  v.true <- c()
  for (i in 1:(p - 1)){
    v.estimate <- c(v.estimate, Edgeshat[i, (i + 1):p])
    v.true <- c(v.true, SuppDelta[i, (i + 1):p])
  }
  
  n <- length(v.estimate)
  TP <- rep(0, n)
  TN <- rep(0, n)
  FP <- rep(0, n)
  FN <- rep(0, n)
  
  TP <- as.double(sum(v.estimate & v.true))
  TN <- as.double(sum(!v.estimate & !v.true))
  FP <- as.double(sum(v.estimate & !v.true))
  FN <- as.double(sum(!v.estimate & v.true))
  
  TPR <- TP / (TP + FN)
  FPR <- FP / (FP + TN)
  
  return(c(FPR, TPR))
}
ROC_Nai = function(thetax_e, thetay_e, SuppDelta,p,tol){
  ThetaX <- thetax_e
  ThetaY <- thetay_e
  Edgeshat <- matrix(0, nrow=p, ncol=p)
  
  for (i in 1:p){
    for (j in 1:p){
      diff <- abs(ThetaX[i, j] - ThetaY[i, j])
      if (diff>tol){
        Edgeshat[i,j] <- 1
      }
    }
  }
  
  v.estimate <- c()
  v.true <- c()
  for (i in 1:(p - 1)){
    v.estimate <- c(v.estimate, Edgeshat[i, (i + 1):p])
    v.true <- c(v.true, SuppDelta[i, (i + 1):p])
  }
  
  n <- length(v.estimate)
  TP <- rep(0, n)
  TN <- rep(0, n)
  FP <- rep(0, n)
  FN <- rep(0, n)
  
  TP <- as.double(sum(v.estimate & v.true))
  TN <- as.double(sum(!v.estimate & !v.true))
  FP <- as.double(sum(v.estimate & !v.true))
  FN <- as.double(sum(!v.estimate & v.true))
  
  TPR <- TP / (TP + FN)
  FPR <- FP / (FP + TN)
  
  return(c(FPR, TPR))
}
get_element_fast = function(data, trace = T){
  # this function is used to get gamma and g for single normal conditional graph (x1,..xp, x1^2.. xp^2,x_jx_k, x^2_jx^2_k)
  if(is.null(dim(data))){data = matrix(data, nrow = 1, ncol = length(data))}
  
  n = dim(data)[1]
  p = dim(data)[2]
  s = 2*(p + choose(p,2))
  
  pos_mat = matrix(1, p,p)
  pos_mat[upper.tri(pos_mat)] = 0
  diag(pos_mat) = 0
  tmat_lower = which(pos_mat !=0, arr.ind = T)
  tmat = tmat_lower[,c(2,1)]
  
  tmp_deri = function(data,j){
    
    deri_top = matrix(0,p,n)
    deri_top[j,] = 1
    
    deri_mid = matrix(0,p,n)
    deri_mid[j,] = 2*data[,j]
    
    deri_low = t(data[,tmat[,1]]*data[,tmat[,2]])
    if(n==1){deri_low = t(deri_low)}
    deri_low[which((tmat[,1]!=j)&(tmat[,2]!=j)),] = 0
    deri_low = sweep(deri_low, 2, data[,j], '/')
    deri_low = 2*deri_low
    
    deri_butt = t(data[,tmat[,1]]^2*data[,tmat[,2]]^2)
    if(n==1){deri_butt = t(deri_butt)}
    deri_butt[which((tmat[,1]!=j)&(tmat[,2]!=j)),] = 0
    deri_butt = sweep(deri_butt, 2, data[,j], '/')
    deri_butt = 2*2*deri_butt
    
    deri_full = rbind(deri_top, deri_mid, deri_low, deri_butt)
    deri_ii_full = rbind(matrix(0,p,n), sweep(deri_mid, 2, data[,j],'/'), matrix(0,dim(deri_low)[1],n), sweep(deri_butt, 2, data[,j], '/'))
    
    Gammaj = deri_full%*%t(deri_full)
    Gj = rowMeans(deri_ii_full)
    return(list(Gamma = Gammaj, G = Gj))
  }
  
  Gamma = 0
  G = 0
  for(j in 1:p){
    deri_res = tmp_deri(data = data,j = j)
    Gamma = Gamma + deri_res$Gamma
    G = G + deri_res$G
    if(trace == T){print(j)}
  }
  
  Gamma = Gamma/n
  G = G
  
  return(list(Gamma = Gamma, G = G))
}
get_element_block = function(data, trace = T,weights,boot){
  if(is.null(dim(data))){data = matrix(data, nrow = 1, ncol = length(data))}
  n = dim(data)[1]
  p = dim(data)[2]
  
  
  
  Gamma = matrix(0, 2*(p^2),2*(p^2))
  G = rep(0, 2*p^2)
  deri_list=list()
  for(j in 1:p){
    if(trace ==T){
      print(j)}
    pos = ((2*p*(j-1))+1) : ((2*p*(j-1)) + 2*p)
    
    deri_raw = NULL
    for(i in 1:n){
      datai = data[i,]
      dataij = datai[j]
      
      linmat = (datai)%*%t(datai)
      linmat_j = linmat[,j]
      linmat_deri_j = 2*(linmat_j/dataij)
      
      quamat = (datai^2)%*%t(datai^2)
      quamat_j = quamat[,j]
      quamat_deri_j = 2*2*quamat_j/(dataij)
      quamat_deri_j[j] = 1
      
      if(boot==T){
        deri_j = c(linmat_deri_j,quamat_deri_j)*sqrt(as.complex(weights[i]))
      }else{
        deri_j = c(linmat_deri_j,quamat_deri_j)
      }
      
      deri_raw = cbind(deri_raw, deri_j)
      
      lin_deri_ii = rep(0,p)
      lin_deri_ii[j] = 2
      qua_deri_ii = 2*2*quamat_j/(dataij^2)
      qua_deri_ii[j] = 0
      
      deri_j_ii = c(lin_deri_ii, qua_deri_ii)*weights[i]
      G[pos] = G[pos]+deri_j_ii/n
    }
    deri_list[[j]] = deri_raw
    Gamma[pos,pos] = deri_raw%*%t(deri_raw)/n
  }
  return(list(Gamma = Gamma,G = G, deri= deri_list))
}
ADMM_gen_block = function(rho, lambda, alpha, threshold, maxit, elem_x, elem_y){
  
  tmp_mat = diag(1,p)
  diag_set_raw = NULL
  for(j in 1:p){
    diag_set_raw = c(diag_set_raw, c(tmp_mat[j,], tmp_mat[j,]))
  }
  diag_set = which(diag_set_raw!=0)
  
  
  gamma_x = elem_x$Gamma
  gamma_y = elem_y$Gamma
  g_x = elem_x$G
  g_y = elem_y$G
  s = dim(gamma_x)[1]
  
  Inv_mat1 = solve(gamma_x + gamma_y + diag(rho, s))
  Inv_mat2 = solve(gamma_y + diag(rho,s))
  
  theta_x = rep(0, s)
  delta = rep(0, s)
  Z1 = rep(0, s)
  Z2 = rep(0, s)
  U1 = rep(0, s)
  U2 = rep(0, s)
  
  error = 1 + threshold
  i = 0
  while((error > threshold)&(i<=maxit)){
    
    theta_0 = c(theta_x, delta)
    theta_x = Inv_mat1%*%(gamma_y%*%delta - g_x - g_y + rho*(Z1 - U1))
    delta = Inv_mat2%*%(gamma_y%*%theta_x + g_y + rho*(Z2 - U2))
    
    A1 = theta_x + U1
    A2 = delta + U2
    
    Z1 = sign(A1)*pmax(0, abs(A1)-lambda/rho)
    Z1[diag_set] = A1[diag_set]
    Z2 = sign(A2)*pmax(0, abs(A2)-lambda*alpha/rho)
    
    U1 = U1 + (theta_x - Z1)
    U2 = U2 + (delta - Z2)
    
    theta_1 = c(theta_x, delta)
    dif = theta_1 - theta_0
    error = norm(dif, type = '2')/norm(theta_0, type = '2')
    i = i+1
    print(paste(lambda, '+',alpha, '+',i, '+', error, sep = ' '))
  }
  return(list(theta_x = Z1, delta = Z2))
}
setwd('~')
setwd('Research/dsm')

args = commandArgs(TRUE)

n = as.numeric(args[[1]])
p = as.numeric(args[[2]])
alpha = as.numeric(args[[3]])
seed = as.numeric(args[[4]])


domain = make_domain(type="R", p=p)
rho = 1
dm = 1

#lambda_seq = seq(from = 1e-3, to = 1, length = 20)
#lambda_seq1 = seq(from = 1e-4, to = 0.1, length = 20)
lambda_seq = seq(from = 1e-3, to = 3, length = 100)
lambda_seq1 = seq(from = 1e-4, to = 0.1, length = 100)
attributes = list(n = n, p = p, alpha = alpha, seed = seed)

beta = matrix(0,p,p)
for(j in 2:p){
  beta[j,j-1] = -0.04
  beta[j-1, j] = -0.04
  if(j>=3){
    beta[j, j-2] = -0.04
    beta[j-2, j] = -0.04
  }
}
diag(beta) = 8/50

beta2 = matrix(0,p,p)
for(j in 6:p){
  beta2[j-5,j] = -0.04
  beta2[j,j-5] = -0.04
}
diag(beta2) = -1

Delta_true = matrix(0,p,p)
for(j in 1:10){
  Delta_true[1,3*j] = 0.5
  Delta_true[j*3,1] = 0.5
}
beta_y = beta - Delta_true
Delta_v = Delta_true[lower.tri(Delta_true)]

set.seed(seed)
x_data = get_data_aug(2*n, beta, beta2, 1000, 10)
y_data = get_data_aug(2*n, beta_y, beta2, 1000, 10)
eles_x = get_element_block(x_data[1:n,],T,rep(1,n),F)
eles_y = get_element_block(y_data[1:n,],T,rep(1,n),F)
#eles_x_2 = get_element_aug(x_data[(n+1):(2*n),])
#eles_y_2 = get_element_aug(y_data[(n+1):(2*n),])

#data_list = list(x_data_all = x_data, y_data_all = y_data, beta = beta, beta_y = beta_y, eles_x = eles_x, eles_y = eles_y,
#                 eles_x_2 = eles_x_2, eles_y_2 = eles_y_2)
#check_dir_data = paste('Data_AugNCSimuP1/n', n, 'p',p,'/',sep = '')
#if(!dir.exists(check_dir_data)){dir.create(check_dir_data)}
#saveRDS(data_list, paste(check_dir_data, 'seed', seed, '.rds', sep = ''))




#ROC for differential graph
dif_supp = (Delta_true!=0)
ROC_res = matrix(nrow = length(lambda_seq), ncol = 2)
ROC_glasso = matrix(nrow = length(lambda_seq), ncol = 2)
ROC_gen = matrix(nrow = length(lambda_seq), ncol = 2)
res_delta = list()
glasso_delta = list()
gen_delta = list()
Loss_res = matrix(nrow = length(lambda_seq), ncol = 2)
Loss_glasso = matrix(nrow = length(lambda_seq), ncol = 2)
Loss_gen = matrix(nrow = length(lambda_seq), ncol = 2)

for(i in 1:length(lambda_seq)){
  
  res = ADMM_gen_block(rho, lambda_seq[i], alpha, 1e-3, 1000, eles_x, eles_y)
  delta_res = res$delta
  delta_sm = matrix(0,p,p)
  delta_sm2 = matrix(0,p,p)
  for(j in 1:p){
    pos = ((2*p*(j-1))+1) : ((2*p*(j-1)) + 2*p)
    delta_sm[,j] = delta_res[pos][1:p]
    delta_sm2[,j] = delta_res[pos][-(1:p)]
  }
  #delta_sm = delta_sm*(abs(delta_sm)<=t(abs(delta_sm))) + t(delta_sm)*(t(abs(delta_sm))<=abs(delta_sm))
  #delta_sm2 = delta_sm2*(abs(delta_sm2)<=t(abs(delta_sm2))) + t(delta_sm2)*(t(abs(delta_sm2))<=abs(delta_sm2))
  delta_sm = (delta_sm + t(delta_sm))/2
  delta_sm2 =(delta_sm2 + t(delta_sm2))/2
  sm_est = (delta_sm!=0)|(delta_sm2!=0)
  
  
  glx = glasso(cov(x_data[1:n,]), rho = lambda_seq1[i], penalize.diagonal = F)
  gly = glasso(cov(y_data[1:n,]), rho = lambda_seq1[i], penalize.diagonal = F)
  
  elts_gauss_x = get_elts(NULL, x_data[1:n,], setting="gaussian", centered=T, domain=domain,
                          diagonal_multiplier=dm,  use_C=TRUE)
  elts_gauss_y = get_elts(NULL, y_data[1:n,], setting="gaussian", centered=T, domain=domain,
                          diagonal_multiplier=dm,  use_C=TRUE)
  est_gen_x = estimate(x_data, setting="gaussian", elts=elts_gauss_x, centered=T, 
                       domain=domain, symmetric="symmetric", lambda1s = lambda_seq1[i],
                       verbose=TRUE, tol=1e-3, maxit=1000, BIC_refit=T, diagonal_multiplier=dm,
                       return_raw=TRUE)
  est_gen_y = estimate(y_data, setting="gaussian", elts=elts_gauss_y, centered=T, 
                       domain=domain, symmetric="symmetric", lambda1s = lambda_seq1[i],
                       verbose=TRUE, tol=1e-3, maxit=1000, BIC_refit=T, diagonal_multiplier=dm,
                       return_raw=TRUE)
  
  ROC_res[i,] = ROC_score(sm_est, dif_supp, p)
  ROC_glasso[i,] = ROC_Nai(glx$wi, gly$wi, dif_supp, p, 1e-3)
  ROC_gen[i,] = ROC_Nai(est_gen_x$raw_estimates[[1]], est_gen_y$raw_estimates[[1]], dif_supp, p, 1e-3)
  
  res_delta[[i]] = delta_sm
  glasso_delta[[i]] = glx$wi - gly$wi
  gen_delta[[i]] = est_gen_x$raw_estimates[[1]] - est_gen_y$raw_estimates[[1]]
  
  loss_res = norm(delta_sm - Delta_true, type = 'F')/sqrt(p)
  Loss_res[i,] = c(sum(delta_sm!=0)/2, loss_res)
  
  delta_gl = glx$wi - gly$wi
  diag(delta_gl) = 0
  loss_gl = norm(delta_gl - Delta_true, type = 'F')/sqrt(p)
  Loss_glasso[i,] = c(sum(delta_gl!=0)/2, loss_gl)
  
  delta_gen = est_gen_x$raw_estimates[[1]] - est_gen_y$raw_estimates[[1]]
  diag(delta_gen) = 0
  loss_gen = norm(delta_gen - Delta_true, type = 'F')/sqrt(p)
  Loss_gen[i,] = c(sum(delta_gen!=0)/2, loss_gen)
  
}
plot(ROC_res, type = 'l', col = 'blue', xlab = 'FPR', ylab = 'TPR', ylim = c(0,1),main = paste('Delta Graph for ','n = ', n, ', p = ', p, sep = ''))
lines(ROC_gen, col = 'red', lty = 2)
lines(ROC_glasso, col = 'green', lty = 3)

res_all = list(ROC_res = ROC_res, ROC_glasso = ROC_glasso, ROC_gen = ROC_gen,
               Loss_res = Loss_res, Loss_glasso = Loss_glasso, Loss_gen = Loss_gen, res_delta = res_delta, glasso_delta = glasso_delta, gen_delta = gen_delta,
               attributes = attributes)
check_dir_res = paste('Res_AugNCSimuP1_block/n', n, 'p',p,'_syl/',sep = '')
if(!dir.exists(check_dir_res)){dir.create(check_dir_res)}
saveRDS(res_all, paste(check_dir_res, 'seed', seed,'alpha',alpha, '.rds', sep = ''))
