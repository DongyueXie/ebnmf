set.seed(123)
n = 120
p = 300
K= 3
L = matrix(0, nrow=n, ncol=K)
FF = matrix(0, nrow=K, ncol=p)
L[1:(n/3),1] = 1
L[((n/3)+1):(2*n/3),2] = 1
L[((2*n/3)+1):n,3] = 1
L = L + matrix(runif(n*K,0,0.5),nrow=n)
FF[1,1:(p/3)] = 1+10
FF[2,((p/3)+1):(2*p/3)] = 1+10
FF[3,((2*p/3)+1):p] = 1+10
lambda = L %*% FF
sigma2=0.5
lambda_noised = lambda*exp(matrix(rnorm(n*p,0,sqrt(sigma2)),nrow=n))
X = matrix(rpois(n=length(lambda_noised),lambda_noised),nrow=n)

library(ebnmf)
res= ebnmf(X,3,ebpm.fn = c(ebpm::ebpm_gamma,smashrgen::BMSM),smooth_F = T,tol=1e-5,warm_start = F)
res= ebnmf(X,3,ebpm.fn = c(ebpm::ebpm_gamma,smashrgen::ebpm_pois_sgp),smooth_F = T,tol=1e-5,maxiter = 1000,warm_start = T,
           smooth_control = list(maxiter=3))

plot(res$elbo_trace)
plot(res$EF[,1],type='l')
plot(res$EF[,2],type='l')
plot(res$EF[,3],type='l')

plot(res$EL[,1])
plot(res$EL[,2])
plot(res$EL[,3])
