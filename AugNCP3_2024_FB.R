#library(genscore)
library(mvtnorm)
library(Matrix)
library(MASS)
#library(glasso)
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
mat_sym = function(mat){
  s = mat*(abs(mat)<=t(abs(mat))) + t(mat)*(t(abs(mat))<abs(mat))
  s
}

setwd('~')
setwd('Research/dsm')

args = commandArgs(TRUE)

n = as.numeric(args[[1]])
p = as.numeric(args[[2]])
alpha = as.numeric(args[[3]])
seed = as.numeric(args[[4]])
m = 10

rho = 1
dm = 1
lambda = 5*sqrt(log(p)/n)
attributes = list(n = n, p = p, alpha = alpha, seed = seed)
lambda2 = 0.5*sqrt(log(p)/n)
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

betax_vec = c(diag(beta), diag(beta2), beta2[lower.tri(beta2)], beta[lower.tri(beta)])
betay_vec = c(diag(beta_y), diag(beta2), beta2[lower.tri(beta2)], beta_y[lower.tri(beta_y)])
Delta_v = betax_vec - betay_vec

#Delta_v = Delta_true[lower.tri(Delta_true)]

set.seed(seed)
x_data = get_data_aug(2*n, beta, beta2, 1000, 10)
y_data = get_data_aug(2*n, beta_y, beta2, 1000, 10)

eles_x_aug_block = get_element_block(x_data[1:n,],F,rep(1,n),F)
eles_y_aug_block = get_element_block(y_data[1:n,],F,rep(1,n),F)
eles_x_aug_block1 = get_element_block(x_data[(n+1):(2*n),], F, rep(1, n),F)
eles_y_aug_block1 = get_element_block(y_data[(n+1):(2*n),], F, rep(1,n),F)


s = dim(eles_x_aug_block$Gamma)[1]

sm_est_aug = ADMM_gen_block(rho, lambda, alpha, 1e-3, 1000, eles_x_aug_block, eles_y_aug_block)
Gamma_hat_aug = rbind(cbind(eles_x_aug_block$Gamma+eles_y_aug_block$Gamma, -eles_y_aug_block$Gamma),cbind(-eles_y_aug_block$Gamma, eles_y_aug_block$Gamma))
Ghat_aug = c(eles_x_aug_block$G + eles_y_aug_block$G, -eles_y_aug_block$G)

#xmatlin=matrix(0,p,p)
#xmatqua = matrix(0,p,p)
#dmatlin = matrix(0,p,p)
#dmatqua = matrix(0,p,p)
#for(j in 1:p){
#  pos = ((2*p*(j-1))+1) : ((2*p*(j-1)) + 2*p)
#  xmatlin[,j] = sm_est_aug$theta_x[pos][1:p]
#  xmatqua[,j] = sm_est_aug$theta_x[pos][-(1:p)]
#  dmatlin[,j] = sm_est_aug$delta[pos][1:p]
#  dmatqua[,j] = sm_est_aug$delta[pos][-(1:p)]
#}
#xmatlin = mat_sym(xmatlin)
#xmatqua = mat_sym(xmatqua)
#dmatlin = mat_sym(dmatlin)
#dmatqua = mat_sym(dmatqua)

#xlong = NULL
#dlong = NULL
#for(j in 1:p){
#  xlong = c(xlong, xmatlin[,j], xmatqua[,j])
#  dlong = c(dlong, dmatlin[,j], dmatqua[,j])
#}
#thetahat_long = c(xlong, dlong)
thetahat_long = c(sm_est_aug$theta_x,sm_est_aug$delta)

kx = clime_block((eles_x_aug_block1$Gamma),lambda2,p = p)
ky = clime_block(eles_y_aug_block1$Gamma, lambda2,p = p)
Mhat_aug = rbind(cbind(kx,kx), cbind(kx, kx+ky))
M = Mhat_aug[-(1:(2*p^2)),]

theta_d_aug = thetahat_long - Mhat_aug%*%(Gamma_hat_aug%*%thetahat_long + Ghat_aug)
delta_d_aug = theta_d_aug[-(1:(2*p^2))]


#ThetaX = c(diag(beta), diag(beta2), beta2[lower.tri(beta2)],beta[lower.tri(beta)])
#Delta_true = c(rep(0, 2*p+length(beta2[lower.tri(beta2)] )), Delta_v)
#theta_true = c(ThetaX, Delta_true)

var_e = rep(0, 2*p^2)
var_e = rep(0, 2*p^2)
thetax_long = sm_est_aug$theta_x
delta_long = sm_est_aug$delta
for(i in 1:n){
  x_i = x_data[i,]
  y_i = y_data[i,]
  
  eles_xi = get_element_block(x_i,F,1,F)
  eles_yi = get_element_block(y_i,F,1,F)
  
  z1 = kx%*%(eles_xi$Gamma%*%thetax_long + eles_xi$G)
  z2 = ky%*%(eles_yi$Gamma%*%(thetax_long - delta_long) + eles_yi$G)
  var_e = var_e + z1^2 + z2^2
  print(i)
}
var_e = var_e/n


#tmp_mat = diag(1,p)
#diag_set_raw = NULL
#for(j in 1:p){
#  diag_set_raw = c(diag_set_raw, c(tmp_mat[j,], tmp_mat[j,]))
#}
#diag_set = which(diag_set_raw!=0)


delta_d = delta_d_aug
var_d = var_e

tmp_delta = NULL
for(i in 1:p){
  tmp_delta = c(tmp_delta,rep(0,p), Delta_true[,i])
}

Dif_Supp = (tmp_delta!=0)

clb = delta_d + sqrt(var_d)*qnorm(0.025)/sqrt(n)
cub = delta_d + sqrt(var_d)*qnorm(0.975)/sqrt(n)

supp = which(Dif_Supp!=0, arr.ind = T)
supp_c = which(Dif_Supp==0, arr.ind = T)

mean(sqrt(n)*abs(delta_d[supp_c])/sqrt(var_d[supp_c])>qnorm(0.975))
typeI = sum((tmp_delta[supp_c]<clb[supp_c]) | (tmp_delta[supp_c] > cub[supp_c]))/length(supp_c)
typeII = sum((0>=clb[supp])*(0<=cub[supp]))/length(supp)
tmp = cub - clb
mean(tmp[supp])
mean(tmp[supp_c])
results = list(supp =supp,  supp_c = supp_c,typeI = typeI, typeII = typeII, ci = list(clb, cub))

check_dir_res = paste('Res_AugNCSimuP3/n', n, 'p',p,'_fb/',sep = '')
if(!dir.exists(check_dir_res)){dir.create(check_dir_res)}
saveRDS(results, paste(check_dir_res,'alpha',alpha, 'seed', seed ,'.rds', sep = ''))