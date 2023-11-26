#'@title Smoothed Poisson Topic Model
#'@description This function fits Poisson Topic Model with smooth Loading or Factors
#'@param X count matrix
#'@param K number of factors/ranks
#'@param init initialization methods, default is 'fasttopics'; or provide init as a list with L_init, and F_init.
#'@param maxiter,maxiter_init maximum iterations
#'@param tol stop criteria
#'@param ebpm.fn specify functions to use for solving the poisson subproblems
#'@param fix_F if TRUE, F will not be updated and will be fixed at the input value in init.
#'@param smooth_F whether smooth l or f, must match the functions in ebpm.fn
#'@param smooth_control a list. ebnmf_smooth_control_default() gives default settings.
#'@param convergence_criteria 'mKLabs', or 'ELBO'
#'@return EL,EF: posterior of loadings and factors
#'@examples
#'set.seed(123)
#'n = 120
#'p = 256
#'K= 3
#'L = matrix(0, nrow=n, ncol=K)
#'FF = matrix(0, nrow=K, ncol=p)
#'L[1:(n/3),1] = 1
#'L[((n/3)+1):(2*n/3),2] = 1
#'L[((2*n/3)+1):n,3] = 1
#'L = L + matrix(runif(n*K,0,0.5),nrow=n)
#'FF[1,1:(p/3)] = 1+10
#'FF[2,((p/3)+1):(2*p/3)] = 1+10
#'FF[3,((2*p/3)+1):p] = 1+10
#'lambda = L %*% FF
#'X = matrix(rpois(n=length(lambda),lambda),nrow=n)
#'image(X)
#'@import ebpm
#'@import Matrix
#'@import vebpm
#'@importFrom smashrgen ebps
#'@importFrom smashrgen BMSM
#'@importFrom Rfast rowsums
#'@export

ebnmf = function(X,K,
                  init = 'fasttopics',
                 over_dispersion = FALSE,
                  maxiter=50,
                  maxiter_init = 100,
                  tol=1e-3,
                  ebpm.fn=c(ebpm::ebpm_gamma,smashrgen::ebpm_pois_sgp),
                  fix_F = FALSE,
                  smooth_F = TRUE,
                  smooth_control=list(maxiter=5,m=100),
                  warm_start=TRUE,
                  printevery=10,
                  verbose=TRUE,
                  convergence_criteria = 'ELBO'){

  # remove columns that are all 0, and are at the start or end of the matrices
  while(sum(X[,1])==0){
    cat('Removed first column that are all 0')
    cat('\n')
    X = X[,-1]
  }
  while(sum(X[,ncol(X)])==0){
    cat('Removed last column that are all 0')
    cat('\n')
    X = X[,-ncol(X)]
  }

  start_time = Sys.time()
  #browser()
  n = dim(X)[1]
  p = dim(X)[2]
  n_points = n*p

  x_rs = rowSums(X)
  x_cs = colSums(X)

  X = Matrix::Matrix(X,sparse = T)
  x = Matrix::summary(X)
  non0_idx = cbind(x$i,x$j)

  # smooth_control = modifyList(ebnmf_smooth_control_default(),smooth_control,keep.null = TRUE)
  if(length(ebpm.fn)==1){
    ebpm.fn.l = ebpm.fn
    ebpm.fn.f = ebpm.fn
  }
  if(length(ebpm.fn)==2){
    ebpm.fn.l = ebpm.fn[[1]]
    ebpm.fn.f = ebpm.fn[[2]]
  }

  if(verbose){
    cat('initializing loadings and factors...')
    cat('\n')
  }

  res = ebnmf_init(X,K,init,maxiter_init)

  if(smooth_F){
    res$qf_aug = vector("list", K)
  }
  alpha = res$ql$Elogl[x$i,] + res$qf$Elogf[x$j,]
  exp_offset = rowMaxs(alpha)
  alpha = alpha - outer(exp_offset,rep(1,K),FUN='*')
  alpha = exp(alpha)
  alpha = alpha/rowsums(alpha)

  obj = c()
  obj[1] = -Inf

  if(verbose){
    cat('running iterations')
    cat('\n')
  }

  for(iter in 1:maxiter){

    for(k in 1:K){
      Ez = calc_EZ(x, alpha[,k])
      res = ebnmf_update_rank1(Ez$rs,Ez$cs,k,ebpm.fn.l,ebpm.fn.f,res,fix_F,smooth_F,smooth_control,warm_start,over_dispersion)
    }

    if(over_dispersion){
      alpha_row_scale = c(colSums(tcrossprod(res$qf$Ef,res$ql$El)*res$q_alpha$col$mean))
      fit_alpha_row = ebpm_gamma1(x_rs,alpha_row_scale)
      res$q_alpha$row$mean = fit_alpha_row$posterior$mean
      res$q_alpha$row$mean_log = fit_alpha_row$posterior$mean_log
      res$g_alpha$row = fit_alpha_row$fitted_g
      res$H_alpha[1] = calc_H(x_rs,alpha_row_scale,fit_alpha_row$log_likelihood,fit_alpha_row$posterior$mean,fit_alpha_row$posterior$mean_log)

      alpha_col_scale = c(colSums(tcrossprod(res$ql$El,res$qf$Ef)*res$q_alpha$row$mean))
      fit_alpha_col = ebpm_gamma1(x_cs,alpha_col_scale)
      res$q_alpha$col$mean = fit_alpha_col$posterior$mean
      res$q_alpha$col$mean_log = fit_alpha_col$posterior$mean_log
      res$g_alpha$col = fit_alpha_row$fitted_g
      res$H_alpha[2] = calc_H(x_cs,alpha_col_scale,fit_alpha_col$log_likelihood,fit_alpha_col$posterior$mean,fit_alpha_col$posterior$mean_log)
    }

    alpha = res$ql$Elogl[x$i,] + res$qf$Elogf[x$j,]
    exp_offset = rowMaxs(alpha)
    alpha = alpha - outer(exp_offset,rep(1,K),FUN='*')
    alpha = exp(alpha)
    alpha = alpha/rowsums(alpha)

    saveRDS(list(res=res,alpha=alpha),file=paste("debug/iter",iter,".rds",sep=''))

    if(convergence_criteria == 'mKLabs'){
      obj[iter+1] = mKL(x$x,tcrossprod(res$ql$El,res$qf$Ef)[non0_idx])
      if(verbose){
        if(iter%%printevery==0){
          cat(sprintf('At iter %d, mKL(X,LF) = %f',iter,obj[iter+1]))
          cat('\n')
        }
      }
      if(abs(obj[iter+1]-obj[iter])<=tol){
        break
      }
    }

    if(convergence_criteria=='ELBO'){
      obj[iter+1] = calc_ebnmf_obj_sparse(x,n,p,K,res,non0_idx,alpha)
      if(verbose){
        if(iter%%printevery==0){
          print(sprintf('At iter %d, ELBO: %f',iter,obj[iter+1]))
        }
      }
      if((obj[iter+1]-obj[iter])/n_points<tol){
        break
      }
    }
  }
  if(iter==maxiter & verbose){
    message('Reached maximum iterations')
  }

  # calc elbo(approximated)
  if(verbose){
    cat('wrapping-up')
    cat('\n')
  }
  ldf = poisson_to_multinom(res$qf$Ef,res$ql$El)
  EL = ldf$L
  EF = ldf$FF

  fit = list(EL = EL,
             EF = EF,
             elbo=calc_ebnmf_obj_sparse(x,n,p,K,res,non0_idx,alpha),
             d=ldf$s,
             elbo_trace=obj,
             res = res,
             run_time = difftime(Sys.time(),start_time,units='auto'))
  return(fit)
}


calc_ebnmf_obj_sparse = function(x,n,p,K,res,non0_idx,alpha){
  val = sum(x$x*alpha*(res$ql$Elogl[x$i,]+res$qf$Elogf[x$j,] + res$q_alpha$row$mean_log[x$i] + res$q_alpha$col$mean_log[x$j]-log(alpha)))
  obj = val - sum(tcrossprod(res$q_alpha$row$mean*res$ql$El,res$q_alpha$col$mean*res$qf$Ef)) - sum(lfactorial(x$x)) + sum(res$Hl)+sum(res$Hf) + sum(res$H_alpha)
  return(obj)
}

# calc_ebnmf_obj = function(x,n,p,K,res,non0_idx){
#   val = 0
#   qz = calc_qz(n,p,K,res$ql,res$qf)
#   for(k in 1:K){
#     val = val + qz[,,k]*(matrix(res$ql$Elogl[,k],nrow=n,ncol=p,byrow=F)+matrix(res$qf$Elogf[,k],nrow=n,ncol=p,byrow=T)-log(qz[,,k]))
#   }
#   E1 = sum(x$x*val[non0_idx]) - sum(tcrossprod(res$ql$El,res$qf$Ef))
#
#   return(E1+sum(res$Hl)+sum(res$Hf)-sum(lfactorial(x$x)))
# }

calc_H = function(x,s,loglik,pm,pmlog){
  if(is.null(loglik)){
    H = 0
  }else{
    H = loglik - sum(x*log(s)+x*pmlog-pm*s-lfactorial(x))
  }
  H
}

#'@title rank 1 update of the model

ebnmf_update_rank1 = function(l_seq,f_seq,k,ebpm.fn.l,ebpm.fn.f,res,fix_F,smooth_F,ebpm_control,warm_start,over_dispersion){

  # update l
  #l_scale = sum(res$qf$Ef[,k])
  #browser()
  l_scale = res$q_alpha$row$mean*sum(res$q_alpha$col$mean*res$qf$Ef[,k])
  fit = ebpm.fn.l(l_seq,l_scale)
  res$ql$El[,k] = fit$posterior$mean
  res$ql$Elogl[,k] = fit$posterior$mean_log
  res$Hl[k] = calc_H(l_seq,l_scale,fit$log_likelihood,fit$posterior$mean,fit$posterior$mean_log)
  res$gl[[k]] = fit$fitted_g

  if(!fix_F){
    # update f
    # f_scale = sum(res$ql$El[,k])
    f_scale = res$q_alpha$col$mean*sum(res$q_alpha$row$mean*res$ql$El[,k])
    if(smooth_F){
      if(warm_start){
        fit = ebpm.fn.f(f_seq,f_scale,g_init=res$gf[[k]],q_init=res$qf_aug[[k]],control=ebpm_control)
      }else{
        fit = ebpm.fn.f(f_seq,f_scale,control=ebpm_control)
      }
      res$qf$Ef[,k] = fit$posterior$mean
      res$qf$Elogf[,k] = fit$posterior$mean_log
      res$gf[[k]] = fit$fitted_g
      res$qf_aug[[k]] = fit$posterior
    }else{
      fit = ebpm.fn.f(f_seq,f_scale)
      res$qf$Ef[,k] = fit$posterior$mean
      res$qf$Elogf[,k] = fit$posterior$mean_log
    }

    res$Hf[k] = calc_H(f_seq,f_scale,fit$log_likelihood,fit$posterior$mean,fit$posterior$mean_log)

  }

  return(res)

}




#' #'@title Default parameters of ebpm
#' #'@export
#' ebpm_control_default = function(){
#'   list(pi0 = 'estimate',
#'        g_init = NULL,
#'        fix_g = FALSE,
#'        control =  NULL)
#' }


#'@title Default parameters of smooth split
#'@export
ebnmf_smooth_control_default = function(){
  list(wave_trans='ndwt',
       ndwt_method = "ti.thresh",
       filter.number = 1,
       family = 'DaubExPhase',
       ebnm_params=list(),
       maxiter=1,
       maxiter_vga = 10,
       make_power_of_2='extend',
       vga_tol=1e-3,
       tol = 1e-2,
       warmstart=TRUE,
       convergence_criteria = 'nugabs',
       m_init_method_for_init = 'vga')
}





