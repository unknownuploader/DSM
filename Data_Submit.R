#library(genscore)
library(mvtnorm)
library(Matrix)
library(MASS)
library(clime)
library(maotai)
library(igraph)
library(EnsDb.Hsapiens.v79)
#library(fastclime)
get_elements = function(x_data, y_data){
  nobs = dim(x_data)[1]
  gamma_x = t(x_data)%*%x_data/nobs
  gamma_y = t(y_data)%*%y_data/nobs
  return(list(gamma_x = gamma_x, gamma_y = gamma_y))
}
ADMM_sym1 = function(rho, lambda, alpha, threshold, maxit, elem){
  
  
  #initialize
  gamma_x = elem$gamma_x
  gamma_y = elem$gamma_y
  theta_x = matrix(0,p,p)
  delta = matrix(0,p,p)
  Z1 = matrix(0,p,p)
  Z2 = matrix(0,p,p)
  U1 = matrix(0,p,p)
  U2 = matrix(0,p,p)
  error = 1 + threshold
  i = 0
  while((error > threshold)&(i<=maxit)){
    
    theta_0 = rbind(theta_x, delta)
    
    C1 = diag(2,p) + (gamma_y%*%delta)/2 + (delta%*%gamma_y)/2 + rho*(Z1 - U1)
    A1 = (gamma_x + gamma_y)/2
    B1 = A1 + diag(1,p)
    theta_x = sylvester(A1, B1, C1)
    C2 = diag(-1,p) + (gamma_y%*%theta_x)/2 + (theta_x%*%gamma_y)/2 + rho*(Z2 - U2)
    A2 = gamma_y/2
    B2 = A2 + diag(1,p)
    delta = sylvester(A2, B2, C2)
    
    A1 = theta_x + U1
    A2 = delta + U2
    
    Z1 = sign(A1)*pmax(0, abs(A1)-lambda/rho)
    diag(Z1) = diag(A1)
    Z2 = sign(A2)*pmax(0, abs(A2)-lambda*alpha/rho)
    
    U1 = U1 + theta_x - Z1
    U2 = U2 + delta - Z2
    
    theta_1 = rbind(theta_x, delta)
    dif = theta_1 - theta_0
    error = norm(dif, type = 'F')/norm(theta_0, type = 'F')
    i = i+1
    print(paste(lambda, '+',alpha, '+',i, '+', error, sep = ' '))
  }
  return(list(theta_x = Z1, delta = Z2))
}
get_element = function(data){
  # this function is used to get gamma and g for single normal conditional graph
  if(is.null(dim(data))){data = matrix(data, nrow = 1, ncol = length(data))}
  
  n = dim(data)[1]
  p = dim(data)[2]
  s = 2*p + choose(p,2)
  
  pos_mat = matrix(1, p,p)
  pos_mat[upper.tri(pos_mat)] = 0
  diag(pos_mat) = 0
  tmat_lower = which(pos_mat !=0, arr.ind = T)
  tmat = tmat_lower[,c(2,1)]
  
  tmp_deri = function(x, j){
    # a tmp function to get the derivatives of the vector(x1,..xp, x1^2.. xp^2, x^2_jx^2_k), j is the corresponding coordinate
    p = length(x)
    deri_upper = rep(0,p)
    deri_upper[j] = 1
    
    deri_mid = rep(0,p)
    deri_mid[j] = 2*x[j]
    
    
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
    
    deri_full = c(deri_upper, deri_mid, 2*deri_lower)
    return(deri_full)
  }
  
  tmp_deri_ii = function(x, deri){
    # a tmp function to get the derivatives of the vector(x1,..xp, x1^2.. xp^2, x^2_jx^2_k), j is the corresponding coordinate
    # deri is the vector we got in tmp_deri
    p = length(x)
    tmp_sub = deri[1:p]
    j = which(tmp_sub!=0)
    deri_ii = deri/x[j]
    deri_ii[1:p] = rep(0,p)
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
tr = function(mat){
  a = sum(diag(mat))
  return(a)
}
matdif = function(A,B){
  pos_vec = c()
  if(is.null(A)){
    return(B)
  }else
  {for(i in 1:(dim(A)[1])){
    for(j in 1:(dim(B)[1])){
      if (max(abs(A[i,] - B[j,])) ==0){
        pos_vec = c(pos_vec,j)
      }
    }
  }
    pos_vec = unique(pos_vec)
    if(is.null(pos_vec)){return(B)}else{
      return(B[-pos_vec,])}}
}
clime_block = function(matrix,lambda,p){
  ret_mat = matrix(0, 2*p^2, 2*p^2)
  for(j in 1:p){
    print(j)
    pos = ((2*p*(j-1))+1) : ((2*p*(j-1)) + 2*p)
    submat = matrix[pos,pos]
    clime_obj = clime(submat, lambda = lambda, sigma = T, standardize = F)
    submat_res = clime_obj$Omegalist[[1]]
    ret_mat[pos,pos] = submat_res
  }
  return(ret_mat)
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
load('RealData.Rdata')
load('breast40.Rdata')

rho = 1
dat_X <- exp_breast_sex_dif[,2:(expset_breast_men_dim+1)]
dat_Y <- exp_breast_sex_dif[,(2+expset_breast_men_dim):(1+expset_breast_men_dim+expset_breast_women_dim)]
dat_X = log(dat_X+1)
dat_Y = log(dat_Y+1)
dat_X_final = t((dat_X - rowMeans(dat_X)))
dat_Y_final = t((dat_Y - rowMeans(dat_Y)))
dat_X_final[which(dat_X_final==0)] = 1e-5
dat_Y_final[which(dat_Y_final==0)] = 1e-5

x_data = dat_X_final[1:100,]
x_data1 = dat_X_final[-(1:100),]
y_data = dat_Y_final[1:100,]
y_data1 = dat_Y_final[-(1:100),]

n = dim(x_data)[1]
p = dim(x_data)[2]
alpha = 0.5


#Gaussian
{
  n = dim(x_data)[1]
  p = dim(x_data)[2]
  alpha = 0.5
  lambda = 0.4
  
  eles = get_elements(x_data, y_data)
  eles1 = get_elements(x_data1, y_data1)
  
  sm_est = ADMM_sym1(rho, lambda, alpha, 1e-2, 1000, eles)
  theta_hat = rbind(sm_est$theta_x, sm_est$delta)
  
  # debiasing
  Gamma_hat = matrix(0, 2*p, 2*p)
  Gamma_hat1 = matrix(0, 2*p, 2*p)
  
  Gamma_hat[1:p, 1:p] = eles$gamma_x + eles$gamma_y
  Gamma_hat[1:p, (p+1):(2*p)] =  -eles$gamma_y
  Gamma_hat[(p+1):(2*p), 1:p] = -eles$gamma_y
  Gamma_hat[(p+1):(2*p), (p+1):(2*p)] = eles$gamma_y
  
  Gamma_hat1[1:p, 1:p] = eles1$gamma_x + eles1$gamma_y
  Gamma_hat1[1:p, (p+1):(2*p)] =  -eles1$gamma_y
  Gamma_hat1[(p+1):(2*p), 1:p] = -eles1$gamma_y
  Gamma_hat1[(p+1):(2*p), (p+1):(2*p)] = eles1$gamma_y
  
  Ghat = rbind(diag(-2,p), diag(1,p))
  
  
  clime_objx = clime(x_data1, lambda = 0.2, standardize = F)
  clime_objy = clime(y_data1, lambda = 0.2, standardize = F)
  kx = clime_objx$Omegalist[[1]]
  ky = clime_objy$Omegalist[[1]]
  
  Mhat = rbind(cbind(kx,kx), cbind(kx, kx+ky))
  
  theta_d = theta_hat - Mhat%*%(Gamma_hat%*%theta_hat + Ghat)
  delta_d = theta_d[(p+1):(2*p),]
  
  t = 1
  edges_set = NULL
  edge_all = which(diag(-1,p)>=0, arr.ind = T)
  edges_set_list = list()
  repeat{
    edge_crit = matdif(edges_set, edge_all)
    TB = NULL
    gamma_list = list()
    for(i in 1:n){
      xi = x_data[i,]
      yi = y_data[i,]
      gamma_xi = xi%*%t(xi)
      gamma_yi = yi%*%t(yi)
      
      gamma_i = rbind(cbind(gamma_xi + gamma_yi, -gamma_yi), cbind(-gamma_yi, gamma_yi))
      gamma_list[[i]] = gamma_i
    }
    for(c in 1:500){ #loop for Monte Carlo
      
      leading = 0
      for(i in 1:n){# loop for number of observations
        
        gamma_i = gamma_list[[i]]
        
        leading_i = -Mhat%*%(gamma_i%*%theta_hat + Ghat)*rnorm(1)/sqrt(n)
        leading = leading + leading_i
        #print(paste('c', c, 'i', i))
      }
      print(c)
      leading_sub = leading[(p+1):(p*2),]
      diag(leading_sub) = 0
      TB = c(TB, max(abs(leading_sub[edge_crit])))
    }
    crit = quantile(TB, 0.95)
    
    delta_d = theta_d[(p+1):(2*p),]
    delta_d[edges_set] = 0
    rej_set = which(sqrt(n)*abs(delta_d) >= crit, arr.ind = T)
    
    edges_set = rbind(edges_set, rej_set)
    edges_set_list[[t]] = edges_set
    inf_graph = matrix(0,p,p)
    inf_graph[edges_set] = 1
    diag(inf_graph) = 0
    max_deg = max(colSums(inf_graph))
    rej_set = as.matrix(rej_set)
    t = t+1
    if((t>=20)|(dim(rej_set)[1] == 0)){break}
  }
  a = inf_graph
  tmp_gaussian = graph_from_adjacency_matrix(
    a,
    mode = c( "undirected"),
    weighted = NULL,
    diag = F,
    add.colnames = NULL,
    add.rownames = NA
  )
  ensg = sapply(exp_breast_sex_dif$Name,function(x){chartr(old=".", new=",", x=x)})
  ensg1 = sapply(ensg, function(x){strsplit(x,split=",")[[1]][1]})
  geneIDs = ensembldb::select(EnsDb.Hsapiens.v79, keys= ensg1, keytype = "GENEID", columns = c("SYMBOL","GENEID"))$SYMBOL
  
  
  ly = readRDS('ly.rds')
  
  color_vec_gau = rep('white', 40)
  hubs_gau = which(degree(tmp_gaussian)>6)
  color_vec_gau[hubs_gau] = 'red'
  name_vec_gau = rep(NA,40)
  name_vec_gau[hubs_gau] = geneIDs[hubs_gau]
  size_vec_gau = rep(6,40)
  size_vec_gau[hubs_gau] = 9
  deg_vec_gau = rep(0,40)
  deg_vec_gau[hubs_gau] = c(0.5*pi)
  dist_vec_gau = rep(0,40)
  dist_vec_gau[hubs_gau] = c(1.5)
  plot(tmp_gaussian,vertex.color=color_vec_gau, edge.width = 1.3, vertex.size = size_vec_gau, layout=ly, vertex.label.cex = 0.8,
       vertex.label.dist = dist_vec_gau, vertex.label = name_vec_gau,vertex.label.color = 'black', vertex.label.degree = deg_vec_gau)
  box()
}

#Intermediate Steps to Simplify Gaussian Multiplier Bootstrap steps
#for(i in 1:n){
#  xi = get_element_block(x_data[i,],T,1,F)
#  yi = get_element_block(y_data[i,],T,1,F)
  
#  saveRDS(xi, paste0('Xaug_block/X_',i,'.rds'))
#  saveRDS(yi, paste0('Yaug_block/Y_',i,'.rds'))
#}


#for(i in 1:n){
#  xi = get_element_block(x_data[i,],T,1,F)
#  yi = get_element_block(y_data[i,],T,1,F)
  
#  saveRDS(xi, paste0('Xaug_block_log/X_',i,'.rds'))
#  saveRDS(yi, paste0('Yaug_block_log/Y_',i,'.rds'))
#}



#AugNC_fb
{
  s = 2*(p + p*(p-1)/2)
  lambda =  1.4
  eles_x_aug_block = get_element_block(x_data,T,rep(1,n),F)
  eles_y_aug_block = get_element_block(y_data,T,rep(1,n),F)
  eles_x_aug_block1 = get_element_block(x_data1, T, rep(1, dim(x_data1)[1]),F)
  eles_y_aug_block1 = get_element_block(y_data1, T, rep(1,dim(y_data1)[1]),F)
  
  Gamma_hat_aug = rbind(cbind(eles_x_aug_block$Gamma+eles_y_aug_block$Gamma, -eles_y_aug_block$Gamma),cbind(-eles_y_aug_block$Gamma, eles_y_aug_block$Gamma))
  Ghat_aug = c(eles_x_aug_block$G + eles_y_aug_block$G, -eles_y_aug_block$G)
  
  sm_est_aug = ADMM_gen_block(rho, lambda, alpha, 1e-2, 1000, eles_x_aug_block, eles_y_aug_block)
  thetahat_long = c(sm_est_aug$theta_x, sm_est_aug$delta)
  
  
  kx = clime_block((eles_x_aug_block1$Gamma), 0.7,p)
  ky = clime_block(eles_y_aug_block1$Gamma, 0.7,p)
  Mhat_aug = rbind(cbind(kx,kx), cbind(kx, kx+ky))
  
  theta_d_aug = thetahat_long - Mhat_aug%*%(Gamma_hat_aug%*%thetahat_long + Ghat_aug)
  delta_d_aug = theta_d_aug[-(1:(2*p^2))]
  
  M = Mhat_aug[-(1:(2*p^2)),]
  leading_all = NULL
  for(i in 1:n){# loop for number of observations
    eles_xi = readRDS(paste0('Xaug_block_log/X_',i,'.rds'))
    eles_yi = readRDS(paste0('Yaug_block_log/Y_',i,'.rds'))
    
    gamma_xi = eles_xi$Gamma
    gamma_yi = eles_yi$Gamma
    Gi = c(eles_xi$G + eles_yi$G, -eles_yi$G)
    
    gamma_i = rbind(cbind(gamma_xi + gamma_yi, -gamma_yi), cbind(-gamma_yi, gamma_yi))
    
    leading_i = -M%*%(gamma_i%*%thetahat_long + Gi)/sqrt(n)
    leading_all = cbind(leading_all, leading_i)
    print(i)
  }
  
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
    tb=  NULL
    boot = 500
    for(i in 1:boot){
      leading_remain = leading_all[crit_position,]
      a = leading_remain%*%rnorm(100,0,1)
      tb = c(tb, max(abs(a)))
    }
    crit_aug = quantile((tb), 0.95)
    
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
    
    diag(inf_graph_aug) = 1
    
    for(j in 1:p){
      pos = ((2*p*(j-1))+1) : ((2*p*(j-1)) + 2*p)
      position_set_raw[pos][1:p] = inf_graph_aug[,j]
      position_set_raw[pos][(p+1):(2*p)] = inf_graph_aug[,j]
    }
    position_set = which(position_set_raw!=0)
    print(paste0('t:',t))
  
    max_deg_aug = max(colSums(inf_graph_aug))-1
    indicator = max(abs(inf_graph_aug - inf_graph_aug_tmp))
    t = t+1
    if((t>20)|(indicator==0) ){break}
    
  }
  #image(inf_graph_aug)
  
  tmp_inf_aug = graph_from_adjacency_matrix(
    inf_graph_aug,
    mode = c( "undirected"),
    weighted = NULL,
    diag = F,
    add.colnames = NULL,
    add.rownames = NA
  )
  
  ly = readRDS('ly.rds')
  hubs_aug = which(degree(tmp_inf_aug)>6)
  color_vec = rep('white', 40)
  color_vec[hubs_aug] = 'red'
  size_vec = rep(6,40)
  size_vec[hubs_aug] = 9
  name_vec = NULL
  name_vec[hubs_aug] = geneIDs[hubs_aug]
  deg_vec = rep(0,40)
  deg_vec[hubs_aug] = c(0.5*pi,1.75*pi,2*pi,1.95*pi)
  dist_vec = rep(0,40)
  dist_vec[hubs_aug] = c(1.5,1.7,2.8,2.8)
  plot(tmp_inf_aug,vertex.color=color_vec,layout = ly, edge.width = 1.3, vertex.size = size_vec, vertex.label.cex = 0.8,
       vertex.label.dist = dist_vec, vertex.label = name_vec, vertex.label.color = 'black',vertex.label.degree=deg_vec)
  box()
  plot(tmp_inf_aug,vertex.color=color_vec,layout = ly, edge.width = 1.3, vertex.size = size_vec, vertex.label.cex = 0.8,
       vertex.label = name_vec, vertex.label.color = 'black')
}
