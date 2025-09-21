#Figure 1 in the draft
setwd('~')
setwd('Research/dsm')
n=  500
p = 30
alpha = 0.5
ROC_all  = matrix(0, 100, 2)
ROC_all_gau = matrix(0, 50, 2)
ROC_all_gl = matrix(0, 100, 2)
Loss_all = matrix(0, 100,2)
Loss_all_gau  = matrix(0, 50,2)
Loss_all_gl  = matrix(0, 100,2)
for(s in 1:300){
  seed = s
  check_dir_res = paste('Res_AugNCSimuP1/n', n, 'p',p,'_syl/',sep = '')
  res = readRDS(paste(check_dir_res, 'seed', seed,'alpha',alpha, '.rds', sep = ''))
  
  check_dir_res1 = paste('Res_AugNCGaussianP1/n', n, 'p',p,'/',sep = '')
  res1 = readRDS(paste(check_dir_res1, 'seed', seed, '.rds', sep = ''))
  
  ROC_all = ROC_all + res$ROC_res
  ROC_all_gau  = ROC_all_gau + res1$ROC_gau
  ROC_all_gl = ROC_all_gl + res$ROC_glasso
  print(s)
}
ROC_all = ROC_all/300
ROC_all_gau = ROC_all_gau/300
ROC_all_gl = ROC_all_gl/300


a = list(ROC_all = ROC_all, ROC_all_gau = ROC_all_gau, ROC_all_gl = ROC_all_gl)
saveRDS(a, 'AugNCP1_n500.rds')

a = readRDS('AugNCP1_n500.rds')
plot(a$ROC_all_gau, type = 'b', col = 'blue',lty = 2, lwd = 2, xlab = 'FPR', ylab = 'TPR')
lines(a$ROC_all, type = 'b', col = 'blue', lwd = 2)
lines(a$ROC_all_gl, type = 'l', col = 'blue', lwd = 2)




#Table 1 in the draft
setwd('~')
setwd('Research/dsm')
n = 5000
p = 30
alpha = 0.5
mat_res = matrix(0, 10, 300)
typeI = NULL
typeII = NULL
len_supp = NULL
len_suppc = NULL
e_mat = NULL
for(s in 1:295){
  seed = s
  check_dir_res = paste('Res_AugNCSimuP3/n', n, 'p',p,'_syl/',sep = '')
  res = readRDS(paste(check_dir_res,'alpha',alpha, 'seed', seed ,'.rds', sep = ''))
  supp = res$supp
  Delta_v = res$Delta_v[-(1:(2*p))]
  supp_c = which(Delta_v==0, arr.ind = T)
  supp = which(Delta_v!=0, arr.ind = T)
  clb = res$ci[[1]]
  cub = res$ci[[2]]
  delta = res$delta_d[supp]
  e_mat = cbind(e_mat, delta - .5)
  
  len_supp = c(len_supp, mean(cub[supp] - clb[supp]))
  len_suppc = c(len_suppc, mean(cub[supp_c] - clb[supp_c]))
  typeI = c(typeI, res$typeI)
  typeII = c(typeII, res$typeII)
  
  
  print(s)
}

mean(typeI)
1-mean(typeII)
mean(len_supp)
mean(len_suppc)



#Table 2 in the draft
#AugNCSimuP4
setwd('~')
setwd('Research/dsm')
n = 3000
p = 30
m = 5
alpha = .5
e = 0
for(i in 1:300){
  seed = i
  check_dir_res = paste('Res_AugNCSimuP4_block_full/n', n, 'p',p,'_syl_new/',sep = '')
  res = readRDS(paste(check_dir_res,'alpha',alpha, 'seed', seed,'m',m ,'_I.rds', sep = ''))
  e = e+ res$typeI_rej
  print(i)
}
e/300


#AugNCSimuP4
setwd('~')
setwd('Research/dsm')
n = 1500
p = 30
m = 6
alpha = .5
e = 0
for(i in 1:300){
  seed = i
  check_dir_res = paste('Res_AugNCSimuP4_II_block_full/n', n, 'p',p,'_syl_new/',sep = '')
  res = readRDS(paste(check_dir_res,'alpha',alpha, 'seed', seed,'m',m ,'_II.rds', sep = ''))
  e = e+ res$typeII_rej
  print(i)
}
1-e/300
