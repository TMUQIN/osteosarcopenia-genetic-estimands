#!/usr/bin/env Rscript
# Independent R implementation: actual column-major lower triangle; Richardson Jacobian.
args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args)==2)
input_dir <- normalizePath(args[[1]], mustWork=TRUE)
out <- normalizePath(args[[2]], mustWork=TRUE)
rows <- list()
metric <- function(x,w,v,standardize=FALSE,partial=FALSE) {
  S <- matrix(0,8,8)
  S[lower.tri(S,diag=TRUE)] <- x
  S[upper.tri(S)] <- t(S)[upper.tri(S)]
  if (standardize) {
    scales <- rep(1,8); scales[2:5] <- sqrt(diag(S)[2:5])
    S <- S/outer(scales,scales)
  }
  a <- c(0,w,1-w,0,0,0,0,0)
  b <- c(0,0,0,v,1-v,0,0,0)
  vf <- drop(t(a)%*%S%*%a)
  vm <- drop(t(b)%*%S%*%b)
  gamma <- drop(t(a)%*%S%*%b)/vm
  u <- a-gamma*b
  vu <- drop(t(u)%*%S%*%u)
  cf <- sum(a*S[,1]); cu <- sum(u*S[,1]); va <- S[1,1]
  rf <- cf/sqrt(vf*va); ru <- cu/sqrt(vu*va)
  if (partial) {
    # Independently residualize both ALM and F on M using Schur complements.
    ma <- sum(b*S[,1])
    return(c(partial_rg_F_ALM_given_M=cu/sqrt(vu*(va-ma^2/vm))))
  }
  c(gamma=gamma,Var_g_F=vf,Var_g_U=vu,rg_F_ALM=rf,rg_U_ALM=ru,
    Delta_rg=ru-rf,Delta_Cov=cu-cf,variance_retention=vu/vf,covariance_retention=cu/cf)
}
for(side in c('right','left')) {
  obj <- readRDS(file.path(input_dir,paste0(side,'_joint_8trait_ldsc.rds')))
  x <- obj$S[lower.tri(obj$S,diag=TRUE)]
  for(w in c(0,.25,.5,.75,1)) for(v in c(0,.25,.5,.75,1)) {
    f <- function(z)metric(z,w,v)
    a <- f(x); J <- numDeriv::jacobian(f,x,method='Richardson')
    cov <- J%*%obj$V%*%t(J)
    stopifnot(all(diag(cov)>0))
    rows[[length(rows)+1]] <- data.frame(side=side,analysis='grid',w=w,v=v,
      metric=names(a),estimate=as.numeric(a),SE=sqrt(diag(cov)))
  }
  for(kind in c('genetic_SD','partial_rg')) {
    f <- function(z)metric(z,.5,.5,standardize=kind=='genetic_SD',partial=kind=='partial_rg')
    a <- f(x); J <- numDeriv::jacobian(f,x,method='Richardson')
    cov <- J%*%obj$V%*%t(J)
    rows[[length(rows)+1]] <- data.frame(side=side,analysis=kind,w=.5,v=.5,
      metric=names(a),estimate=as.numeric(a),SE=sqrt(diag(cov)))
  }
}
write.table(do.call(rbind,rows),file.path(out,'R2_INDEPENDENT_R_CHECK.tsv'),sep='\t',
  row.names=FALSE,quote=FALSE)
writeLines(c(R.version.string,paste('numDeriv',packageVersion('numDeriv')),
  paste('GenomicSEM',packageVersion('GenomicSEM'))),file.path(out,'R2_R_SESSION_INFO.txt'))
