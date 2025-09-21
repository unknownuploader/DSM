setwd('~')

args = commandArgs(TRUE)

n = as.numeric(args[[1]])
p = as.numeric(args[[2]])
alpha = as.numeric(args[[3]])
seed = as.numeric(args[[4]])
m = as.numeric(args[[5]])

library(mvtnorm)
library(Matrix)
library(MASS)
library(parallel)
library(doParallel)
library(clime)
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
get_element_fast = function(data, trace = T, get_deri){
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
    return(list(Gamma = Gammaj, G = Gj, deri_i = deri_full, deri_ii = deri_ii_full))
  }
  
  Gamma = 0
  G = 0
  deri_list = list()
  deri_ii_list = list()
  for(j in 1:p){
    deri_res = tmp_deri(data = data,j = j)
    Gamma = Gamma + deri_res$Gamma
    G = G + deri_res$G
    if(get_deri == T){
      deri_list[[j]] = deri_res$deri_i
      deri_ii_list[[j]] = deri_res$deri_ii
    }
    if(trace == T){print(j)}
  }
  
  Gamma = Gamma/n
  G = G
  if(get_deri == T){
    return(list(Gamma = Gamma, G = G, deri = deri_list, deri_ii = deri_ii_list))
  }else{
    return(list(Gamma = Gamma, G = G))}
}
get_element_single = function(data,i,dir){
  dir_single = paste0(dir,'/',i,'.rds')
  eles = get_element_fast(data[i,],trace = F)
  saveRDS(eles, dir_single)
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
dynamic_p4_I = function(n,p,alpha,seed,m,basedir,l,l2){
  setwd('~')
  mc_function_dir = function(c, crit_pos,M,log_dir){
    write.table(paste('Running ',c,'th',' Monte Carlo...'), file = log_dir,eol = '\n', col.names = F, row.names = F, append = T)
    #print(c)
    Gauweights = rnorm(n,0,1)
    x_ele = get_element_block(x_data[1:n,],F,Gauweights,T)
    y_ele = get_element_block(y_data[1:n,],F,Gauweights,T)
    
    Gamma_boot = rbind(cbind(x_ele$Gamma + y_ele$Gamma,-y_ele$Gamma),cbind(-y_ele$Gamma, y_ele$Gamma))
    G_boot = c(x_ele$G+y_ele$G, -y_ele$G)
    
    
    leading = -M%*%(Gamma_boot%*%thetahat_long + G_boot)*sqrt(n)
    leading = as.numeric(leading)
    leading_sub = leading[crit_pos]
    return(max(abs(leading_sub)))
  }
  # check directories
  wd = getwd()
  data_dir = paste0(wd,'/' ,basedir, '/Data_AugNCSimuP4_I_block_full')
  if(!dir.exists(data_dir)){dir.create(data_dir)}
  base_dir_data = paste0(data_dir,'/n', n, 'p',p,'m',m,'/',sep = '')
  if(!dir.exists(base_dir_data)){dir.create(base_dir_data)}
  res_dir = paste0(wd,'/' ,basedir, '/Res_AugNCSimuP4_block_full')
  if(!dir.exists(res_dir)){dir.create(res_dir)}
  check_dir_res = paste(res_dir,'/n', n, 'p',p,'_syl_new/',sep = '')
  if(!dir.exists(check_dir_res)){dir.create(check_dir_res)}
  log_dir = paste0(base_dir_data, 'seed',seed,'_log_new.txt')
  
  
  #basic settings
  rho = 1
  lambda = l*sqrt(log(p)/n)
  lambda2 = l2*sqrt(log(p)/n)
  k = 5
  attributes = list(n = n, p = p, alpha = alpha, seed = seed, m = m,l = l)
  
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
  for(j in 1:m){
    Delta_true[1,3*j] = 0.5
    Delta_true[3*j,1] = 0.5
  }
  beta_y = beta - Delta_true
  
  betax_vec = c(diag(beta), diag(beta2), beta2[lower.tri(beta2)], beta[lower.tri(beta)])
  betay_vec = c(diag(beta_y), diag(beta2), beta2[lower.tri(beta2)], beta_y[lower.tri(beta_y)])
  Delta_v = betax_vec - betay_vec
  
  write.table('Getting data....', file = log_dir,eol = '\n', col.names = F, row.names = F)
  
  set.seed(seed)
  x_data = get_data_aug(2*n, beta, beta2, 1000, 10)
  y_data = get_data_aug(2*n, beta_y, beta2, 1000, 10)
  write.table('Calculating elements....', file = log_dir,eol = '\n', col.names = F, row.names = F, append = T)
  
  eles_x_aug_block = get_element_block(x_data[1:n,],F,rep(1,n),F)
  eles_y_aug_block = get_element_block(y_data[1:n,],F,rep(1,n),F)
  eles_x_aug_block1 = get_element_block(x_data[(n+1):(2*n),], F, rep(1, n),F)
  eles_y_aug_block1 = get_element_block(y_data[(n+1):(2*n),], F, rep(1,n),F)
  
  write.table('Finished calculating elements....', file = log_dir,eol = '\n', col.names = F, row.names = F, append = T)
  
  sm_est_aug = ADMM_gen_block(rho, lambda, alpha, 1e-3, 1000, eles_x_aug_block, eles_y_aug_block)
  write.table('Finished estimating....', file = log_dir,eol = '\n', col.names = F, row.names = F, append = T)
  Gamma_hat_aug = rbind(cbind(eles_x_aug_block$Gamma+eles_y_aug_block$Gamma, -eles_y_aug_block$Gamma),cbind(-eles_y_aug_block$Gamma, eles_y_aug_block$Gamma))
  Ghat_aug = c(eles_x_aug_block$G + eles_y_aug_block$G, -eles_y_aug_block$G)
  
  thetahat_long = c(sm_est_aug$theta_x, sm_est_aug$delta)
  
  #thetahat_long = c(thetax_long, delta_long)
  write.table('Finished extension....', file = log_dir,eol = '\n', col.names = F, row.names = F, append = T)
  
  kx = clime_block((eles_x_aug_block1$Gamma),lambda2,p = p)
  ky = clime_block(eles_y_aug_block1$Gamma, lambda2,p = p)
  Mhat_aug = rbind(cbind(kx,kx), cbind(kx, kx+ky))
  M = Mhat_aug[-(1:(2*p^2)),]
  write.table('Finished CLIME....', file = log_dir,eol = '\n', col.names = F, row.names = F, append = T)
  theta_d_aug = thetahat_long - Mhat_aug%*%(Gamma_hat_aug%*%thetahat_long + Ghat_aug)
  delta_d_aug = theta_d_aug[-(1:(2*p^2))]
  
  t = 1
  tmp_mat = diag(1,p)
  position_set_raw = NULL
  for(j in 1:p){
    position_set_raw = c(position_set_raw, c(tmp_mat[j,], tmp_mat[j,]))
  }
  position_set = which(position_set_raw!=0)
  position_all = 1:(2*p^2)
  inf_graph_aug = matrix(0,p,p)
  repeat{
    inf_graph_aug_tmp = inf_graph_aug
    crit_position = setdiff(position_all, position_set)
    TB = c()
    
    for(c in 1:100){
      TB = c(TB, mc_function_dir(c, crit_position,log_dir = log_dir,M = M))
    }
    write.table('MC finished!', file = log_dir,eol = '\n', col.names = F, row.names = F, append = T)
    
    crit_aug = quantile(TB,0.95)
    
    
    delta_matlin = matrix(0,p,p)
    delta_matqua = matrix(0,p,p)
    for(j in 1:p){
      pos = ((2*p*(j-1))+1) : ((2*p*(j-1)) + 2*p)
      vec_sub = delta_d_aug[pos]
      delta_matlin[,j] = vec_sub[1:p]
      delta_matqua[,j] = vec_sub[(p+1):(2*p)]
    }
    
    rej_setlin = which(sqrt(n)*abs(delta_matlin)>=crit_aug, arr.ind = T)
    rej_setqua = which(sqrt(n)*abs(delta_matqua)>=crit_aug, arr.ind = T)
    
    inf_graph_aug[rej_setlin] = 1
    inf_graph_aug[rej_setqua] = 1
    
    inf_graph_aug = inf_graph_aug + t(inf_graph_aug)
    inf_graph_aug = (inf_graph_aug!=0)
    diag(inf_graph_aug) = 1
    
    for(j in 1:p){
      pos = ((2*p*(j-1))+1) : ((2*p*(j-1)) + 2*p)
      position_set_raw[pos][1:p] = inf_graph_aug[,j]
      position_set_raw[pos][(p+1):(2*p)] = inf_graph_aug[,j]
    }
    position_set = which(position_set_raw!=0)
    
    max_deg_aug = max(colSums(inf_graph_aug))-1
    indicator = max(abs(inf_graph_aug - inf_graph_aug_tmp))
    t = t+1
    
    if((max_deg_aug>k)|(indicator==0) ){break}
    
  }
  
  typeI_rej = as.numeric(max_deg_aug > k)
  write.table(max_deg_aug, file = log_dir,eol = '\n', col.names = F, row.names = F, append = T)
  result = list(inf_graph = inf_graph_aug, typeI_rej = typeI_rej, m = m, crit = crit_aug, attributes = attributes)
  saveRDS(result, paste(check_dir_res,'alpha',alpha, 'seed', seed,'m',m ,'_I.rds', sep = ''))
  write.table('Finished!', file = log_dir,eol = '\n', col.names = F, row.names = F, append = T)
  #return(result)
}
dynamic_p4_II = function(n,p,alpha,seed,m,basedir,l,l2){
  setwd('~')
  mc_function_dir = function(c, crit_pos,M,log_dir){
    write.table(paste('Running ',c,'th',' Monte Carlo...'), file = log_dir,eol = '\n', col.names = F, row.names = F, append = T)
    Gauweights = rnorm(n,0,1)
    x_ele = get_element_block(x_data[1:n,],F,Gauweights,T)
    y_ele = get_element_block(y_data[1:n,],F,Gauweights,T)
    
    Gamma_boot = rbind(cbind(x_ele$Gamma + y_ele$Gamma,-y_ele$Gamma),cbind(-y_ele$Gamma, y_ele$Gamma))
    G_boot = c(x_ele$G+y_ele$G, -y_ele$G)
    
    
    leading = -M%*%(Gamma_boot%*%thetahat_long + G_boot)*sqrt(n)
    leading = as.numeric(leading)
    leading_sub = leading[crit_pos]
    return(max(abs(leading_sub)))
  }
  # check directories
  wd = getwd()
  data_dir = paste0(wd,'/' ,basedir, '/Data_AugNCSimuP4_II_block_full')
  if(!dir.exists(data_dir)){dir.create(data_dir)}
  base_dir_data = paste0(data_dir,'/n', n, 'p',p,'m',m,'/',sep = '')
  if(!dir.exists(base_dir_data)){dir.create(base_dir_data)}
  res_dir = paste0(wd,'/' ,basedir, '/Res_AugNCSimuP4_II_block_full')
  if(!dir.exists(res_dir)){dir.create(res_dir)}
  check_dir_res = paste(res_dir,'/n', n, 'p',p,'_syl_new/',sep = '')
  if(!dir.exists(check_dir_res)){dir.create(check_dir_res)}
  log_dir = paste0(base_dir_data, 'seed',seed,'_log_new.txt')
  
  
  #basic settings
  rho = 1
  dm = 1
  lambda = l*sqrt(log(p)/n)
  lambda2 = l2*sqrt(log(p)/n)
  k = 5
  attributes = list(n = n, p = p, alpha = alpha, seed = seed, m = m,l = l)
  
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
  for(j in 1:m){
    Delta_true[1,j*3] = 0.5
    Delta_true[j*3,1] = 0.5
  }
  beta_y = beta - Delta_true
  
  betax_vec = c(diag(beta), diag(beta2), beta2[lower.tri(beta2)], beta[lower.tri(beta)])
  betay_vec = c(diag(beta_y), diag(beta2), beta2[lower.tri(beta2)], beta_y[lower.tri(beta_y)])
  Delta_v = betax_vec - betay_vec
  
  write.table('Getting data....', file = log_dir,eol = '\n', col.names = F, row.names = F)
  
  set.seed(seed)
  x_data = get_data_aug(2*n, beta, beta2, 1000, 10)
  y_data = get_data_aug(2*n, beta_y, beta2, 1000, 10)
  write.table('Calculating elements....', file = log_dir,eol = '\n', col.names = F, row.names = F, append = T)
  
  eles_x_aug_block = get_element_block(x_data[1:n,],T,rep(1,n),F)
  eles_y_aug_block = get_element_block(y_data[1:n,],T,rep(1,n),F)
  eles_x_aug_block1 = get_element_block(x_data[(n+1):(2*n),], T, rep(1, n),F)
  eles_y_aug_block1 = get_element_block(y_data[(n+1):(2*n),], T, rep(1,n),F)
  
  #eles_x_2 = get_element_fast(x_data[(n+1):(2*n),],F,F)
  #eles_y_2 = get_element_fast(y_data[(n+1):(2*n),],F,F)
  write.table('Finished calculating elements....', file = log_dir,eol = '\n', col.names = F, row.names = F, append = T)
  
  sm_est_aug = ADMM_gen_block(rho, lambda, alpha, 1e-3, 1000, eles_x_aug_block, eles_y_aug_block)
  
  Gamma_hat_aug = rbind(cbind(eles_x_aug_block$Gamma+eles_y_aug_block$Gamma, -eles_y_aug_block$Gamma),cbind(-eles_y_aug_block$Gamma, eles_y_aug_block$Gamma))
  Ghat_aug = c(eles_x_aug_block$G + eles_y_aug_block$G, -eles_y_aug_block$G)
  
  thetahat_long = c(sm_est_aug$theta_x, sm_est_aug$delta)
  
  
  kx = clime_block((eles_x_aug_block1$Gamma),lambda2,p=p)
  ky = clime_block(eles_y_aug_block1$Gamma, lambda2,p=p)
  Mhat_aug = rbind(cbind(kx,kx), cbind(kx, kx+ky))
  M = Mhat_aug[-(1:(2*p^2)),]
  
  theta_d_aug = thetahat_long - Mhat_aug%*%(Gamma_hat_aug%*%thetahat_long + Ghat_aug)
  delta_d_aug = theta_d_aug[-(1:(2*p^2))]
  
  t = 1
  tmp_mat = diag(1,p)
  position_set_raw = NULL
  for(j in 1:p){
    position_set_raw = c(position_set_raw, c(tmp_mat[j,], tmp_mat[j,]))
  }
  position_set = which(position_set_raw!=0)
  position_all = 1:(2*p^2)
  inf_graph_aug = matrix(0,p,p)
  repeat{
    inf_graph_aug_tmp=inf_graph_aug
    crit_position = setdiff(position_all, position_set)
    TB = c()
    
    for(c in 1:100){
      TB = c(TB, mc_function_dir(c, crit_position,log_dir = log_dir,M = M))
    }
    write.table('MC finished!', file = log_dir,eol = '\n', col.names = F, row.names = F, append = T)
    
    crit_aug = quantile(TB,0.95)
    
    delta_matlin = matrix(0,p,p)
    delta_matqua = matrix(0,p,p)
    for(j in 1:p){
      pos = ((2*p*(j-1))+1) : ((2*p*(j-1)) + 2*p)
      vec_sub = delta_d_aug[pos]
      delta_matlin[,j] = vec_sub[1:p]
      delta_matqua[,j] = vec_sub[(p+1):(2*p)]
    }
    
    rej_setlin = which(sqrt(n)*abs(delta_matlin)>=crit_aug, arr.ind = T)
    rej_setqua = which(sqrt(n)*abs(delta_matqua)>=crit_aug, arr.ind = T)
    
    inf_graph_aug[rej_setlin] = 1
    inf_graph_aug[rej_setqua] = 1
    
    inf_graph_aug = inf_graph_aug + t(inf_graph_aug)
    inf_graph_aug = (inf_graph_aug!=0)
    diag(inf_graph_aug) = 1
    
    for(j in 1:p){
      pos = ((2*p*(j-1))+1) : ((2*p*(j-1)) + 2*p)
      position_set_raw[pos][1:p] = inf_graph_aug[,j]
      position_set_raw[pos][(p+1):(2*p)] = inf_graph_aug[,j]
    }
    position_set = which(position_set_raw!=0)
    
    max_deg_aug = max(colSums(inf_graph_aug))-1
    indicator = max(abs(inf_graph_aug - inf_graph_aug_tmp))
    t = t+1
    
    if((max_deg_aug>k)|(indicator==0) ){break}
    
  }
  
  typeII_rej = as.numeric(max_deg_aug <= k)
  write.table(max_deg_aug, file = log_dir,eol = '\n', col.names = F, row.names = F, append = T)
  result = list(inf_graph = inf_graph_aug, typeII_rej = typeII_rej, m = m, crit = crit_aug, attributes = attributes)
  saveRDS(result, paste(check_dir_res,'alpha',alpha, 'seed', seed,'m',m ,'_II.rds', sep = ''))
  write.table('Finished!', file = log_dir,eol = '\n', col.names = F, row.names = F, append = T)
  #return(result)
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
clime_block = function(matrix,lambda,p){
  ret_mat = matrix(0, 2*p^2, 2*p^2)
  for(j in 1:p){
    print(j)
    pos = ((2*p*(j-1))+1) : ((2*p*(j-1)) + 2*p)
    submat = matrix[pos,pos]
    clime_obj = clime(submat, lambda = lambda, sigma = T, standardize = F,perturb=F)
    submat_res = clime_obj$Omegalist[[1]]
    ret_mat[pos,pos] = submat_res
  }
  return(ret_mat)
}
basedir = 'Research/dsm'

res = dynamic_p4_II(n = n, p = p,alpha = alpha, seed = seed, m = m, basedir = basedir, l = 5, l2 = 0.5)












