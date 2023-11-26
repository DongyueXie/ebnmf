datax = read.csv('/home/dxie/Downloads/SRSF3.Counts.csv.gz')
rownames(datax)=datax[,1]
datax[1:5,1:5]
datax = datax[,-1]
dim(datax)
datax = as.matrix(datax)
sum(datax==0)/prod(dim(datax))

library(ebnmf)

set.seed(12345)
K = 5
fit = ebnmf(datax,K,ebpm.fn = c(ebpm::ebpm_gamma,smashrgen::ebpm_pois_sgp),
            smooth_F = T,tol=1e-5,maxiter = 50,warm_start = T,
            smooth_control = list(maxiter=10,m=100),printevery = 1)
par(mfrow=c(5,1))
for(k in 1:K){
  plot(fit$EF[,k],type='l')
}

fit_cold = ebnmf(datax,K,ebpm.fn = c(ebpm::ebpm_gamma,smashrgen::ebpm_pois_sgp),
            smooth_F = T,tol=1e-5,maxiter = 50,warm_start = F,
            smooth_control = list(maxiter=10,m=100),printevery = 1)

set.seed(12345)
fit_tm = fastTopics::fit_poisson_nmf(Matrix::Matrix(datax,sparse=T),K)
plot(fit_tm$F[,1],type='l')
plot(fit_tm$F[,2],type='l')
plot(fit_tm$F[,3],type='l')
plot(fit_tm$F[,4],type='l')
plot(fit_tm$F[,5],type='l')


fit2 = ebnmf(datax,K,ebpm.fn = c(ebpm::ebpm_gamma,smashrgen::ebpm_pois_sgp),
            smooth_F = T,tol=1e-8,maxiter = 50,warm_start = T,over_dispersion = F,
            smooth_control = list(maxiter=10,m=100,opt_method='L-BFGS-B'),printevery = 1,init = 'kmeans')

for(k in 1:K){
  plot(fit2$EF[,k],type='l')
}

fit2_5 = ebnmf(datax,K,ebpm.fn = c(ebpm::ebpm_gamma,smashrgen::ebpm_pois_sgp),
             smooth_F = T,tol=1e-5,maxiter = 50,warm_start = T,over_dispersion = T,
             smooth_control = list(maxiter=10,m=100),printevery = 1,init = 'kmeans')

fit3 = ebnmf(datax,K,ebpm.fn = c(ebpm::ebpm_gamma,smashrgen::ebpm_BMSM),
             smooth_F = T,tol=1e-5,maxiter = 50,warm_start = T,over_dispersion = F,
             smooth_control = list(maxiter=10,m=100,opt_method='L-BFGS-B'),printevery = 1,init = 'kmeans')
#################
K=5
X = datax
n = dim(X)[1]
p = dim(X)[2]
n_points = n*p

x_rs = rowSums(X)
x_cs = colSums(X)

X = Matrix::Matrix(X,sparse = T)
x = Matrix::summary(X)
non0_idx = cbind(x$i,x$j)
calc_ebnmf_obj_sparse(x,n,p,K,iter38$res,non0_idx,iter38$alpha)
calc_ebnmf_obj_sparse(x,n,p,K,iter39$res,non0_idx,iter39$alpha)

k = 2
Ez = calc_EZ(x, iter38$alpha[,k])
l_scale = iter38$res$q_alpha$row$mean*sum(iter38$res$q_alpha$col$mean*iter38$res$qf$Ef[,k])
l_seq = Ez$rs
fit = ebpm_gamma(l_seq,l_scale)

f_seq = Ez$cs
f_scale = iter38$res$q_alpha$col$mean*sum(iter38$res$q_alpha$row$mean*fit$posterior$mean)
temp = pois_sgp(f_seq,sc=f_scale,m=100,X_ind=iter38$res$gf[[k]]$X_ind,kernel_param = iter38$res$gf[[k]]$kernel_param,
                mu=iter38$res$gf[[k]]$mu,
                post_mean = iter38$res$qf_aug[[k]]$mean_log_ind,
                V = iter38$res$qf_aug[[k]]$v_log_ind,verbose=T)



