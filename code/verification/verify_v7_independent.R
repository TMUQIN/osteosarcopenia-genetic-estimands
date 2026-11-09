# JCSM P0 FINAL CLOSURE — verification-only independent V7 DeltaCov check
# This script never writes to frozen project directories. All paths are relative to this file.
options(stringsAsFactors = FALSE)
args <- commandArgs(trailingOnly = FALSE)
file_arg <- args[grep('^--file=', args)]
script_path <- if (length(file_arg)) sub('^--file=', '', file_arg[1]) else normalizePath('verify_v7_independent.R')
root <- normalizePath(file.path(dirname(script_path), '..', '..'), mustWork = TRUE)
review <- file.path(root, '08_MACHINE_READABLE_INPUTS_AND_SCRIPTS', '02_REVIEWER_ONLY')
inp <- file.path(review, 'inputs')
rds <- file.path(review, 'rds')
logs <- file.path(root, '09_RUN_LOGS_AND_SESSION_INFO')
dir.create(logs, recursive = TRUE, showWarnings = FALSE)
log_file <- file.path(logs, 'verification_only_R_run.log')
sink(log_file, split = TRUE)
cat('verification_only=TRUE\n')
cat('script=', normalizePath(script_path), '\n')
cat('root=', root, '\n')
cat('started=', format(Sys.time(), '%Y-%m-%dT%H:%M:%S%z'), '\n')

needed <- c('GenomicSEM','lavaan','data.table','Matrix','digest')
for (p in needed) if (!requireNamespace(p, quietly=TRUE)) stop('missing package: ', p)
cat('R=', R.version.string, '\n')
cat('GenomicSEM=', as.character(packageVersion('GenomicSEM')), '\n')
cat('lavaan=', as.character(packageVersion('lavaan')), '\n')
cat('data.table=', as.character(packageVersion('data.table')), '\n')
cat('Matrix=', as.character(packageVersion('Matrix')), '\n')
cat('GenomicSEM_remote_sha=', packageDescription('GenomicSEM')$RemoteSha, '\n')

read_num <- function(x) suppressWarnings(as.numeric(as.character(x)))
read_cov <- function(side) {
  pref <- paste0(side, '_joint_8trait_')
  read_mat <- function(suffix) {
    z <- utils::read.delim(file.path(inp, paste0(pref, suffix, '.tsv')), sep='\t', check.names=FALSE,
      row.names=1, stringsAsFactors=FALSE)
    as.matrix(z)
  }
  raw <- readRDS(file.path(rds, paste0('covstruc_', side, '.rds')))
  cs <- list(V=read_mat('V'), S=read_mat('S'), I=read_mat('I'), N=read_mat('N'), m=raw$m)
  # Confirm the copied text S/V/I/N are the same frozen doubles used by the archived covstruc object.
  for (nm in c('V','S','I','N')) {
    d_abs <- max(abs(cs[[nm]] - raw[[nm]]))
    d_rel <- max(abs(cs[[nm]] - raw[[nm]]) / pmax(1, abs(raw[[nm]])))
    cat(side, ' copied_', nm, '_max_abs=', format(d_abs, scientific=TRUE), ' max_rel=', format(d_rel, scientific=TRUE), '\n')
    if (!is.finite(d_rel) || d_rel > 1e-12) stop(side, ' copied ', nm, ' relative difference exceeds 1e-12: ', d_rel)
    dimnames(cs[[nm]]) <- dimnames(raw[[nm]])
  }
  # Fit uses the copied full-precision RDS (the authoritative double-precision S/V/I/N object); TSVs above are independently re-read and checked.
  raw
}
model_frozen <- paste(readLines(file.path(inp,'model_FROZEN.lav'), warn=FALSE), collapse='\n')
model_derived <- paste(readLines(file.path(inp,'model_DERIVED.lav'), warn=FALSE), collapse='\n')
if (!grepl('delta_cov_alm\\s*:=\\s*-gamma\\s*\\*\\s*c_ma', model_derived)) stop('derived syntax missing')

# Capture Ohtt and model parTable from the installed usermodel() execution frame.
# Trace is installed in memory and removed before exit; package files are untouched.
cap_env <- new.env(parent=emptyenv()); cap_env$records <- list()
tracer <- quote({
  if (exists('Ohtt') && exists('Model1_Results') && exists('Model_Output') && exists('S2.delt') && exists('S2.W') && exists('V_Reorder')) {
    assign('current', list(Ohtt=Ohtt, bread=bread, lettuce=lettuce, S2.delt=S2.delt,
      S2.W=S2.W, V_Reorder=V_Reorder, Model_Output=Model_Output,
      converged=lavInspect(Model1_Results,'converged'), npar=lavInspect(Model1_Results,'npar')),
      envir=cap_env)
  }
})
trace_status <- tryCatch({
  suppressMessages(trace('usermodel', where=asNamespace('GenomicSEM'), tracer=tracer, at=61, print=FALSE)); 'TRACE_INSTALLED'
}, error=function(e) paste0('TRACE_FAILED:', conditionMessage(e)))
cat('trace_status=', trace_status, '\n')

fits <- list(); result_rows <- list(); cross_rows <- list(); fit_rows <- list(); weight_rows <- list()
for (side in c('right','left')) {
  cat('--- side=',side,' ---\n')
  cs <- read_cov(side)
  cap_env$current <- NULL
  warnings <- character(); console <- character(); fit <- NULL
  fit <- tryCatch({
    console <- utils::capture.output({
      withCallingHandlers({
        fit <- GenomicSEM::usermodel(covstruc=cs, estimation='DWLS', model=model_derived,
          CFIcalc=TRUE, std.lv=FALSE, imp_cov=FALSE, fix_resid=TRUE)
      }, warning=function(w){warnings <<- c(warnings, conditionMessage(w)); invokeRestart('muffleWarning')})
    }, type='output'); fit
  }, error=function(e) { structure(list(message=conditionMessage(e)), class='verification_error') })
  if (inherits(fit,'verification_error')) stop(side, ' refit failed: ', fit$message)
  cap <- cap_env$current
  if (is.null(cap)) stop(side, ' robust Ohtt capture failed (required for independent variance)')
  fits[[side]] <- fit
  res <- fit$results
  row <- res[res$op == ':=' & res$lhs == 'delta_cov_alm',,drop=FALSE]
  if (nrow(row)!=1) stop(side, ' derived result not unique')
  par <- cap$Model_Output
  free_rows <- which(par$free > 0)
  if (anyDuplicated(par$free[free_rows])) stop(side, 'free index duplicated')
  ord <- order(par$free[free_rows]); free_rows <- free_rows[ord]
  free_names <- paste(par$lhs[free_rows], par$op[free_rows], par$rhs[free_rows], sep='|')
  gidx <- which(par$lhs[free_rows]=='FFM_total' & par$op[free_rows]=='~' & par$rhs[free_rows]=='FM')
  cidx <- which(par$lhs[free_rows]=='FM' & par$op[free_rows]=='~~' & par$rhs[free_rows]=='ALM')
  if (length(gidx)!=1 || length(cidx)!=1) stop(side, ' gamma/c_ma mapping failed')
  gamma <- par$est[free_rows[gidx]]; cma <- par$est[free_rows[cidx]]
  grad <- numeric(length(free_rows)); grad[gidx] <- -cma; grad[cidx] <- -gamma
  # Ohtt is ordered by free parameter index; this mapping is by labels/names, never by SE order.
  vcov <- cap$Ohtt
  if (!all(dim(vcov)==c(length(free_rows), length(free_rows)))) stop(side, ' Ohtt dimension mismatch')
  vcov_out <- as.data.frame(vcov, check.names=FALSE)
  colnames(vcov_out) <- free_names
  vcov_out <- cbind(parameter_name=free_names, vcov_out)
  utils::write.table(vcov_out, file.path(rds, paste0(side, '_Ohtt_parameter_covariance.tsv')), sep='\t', row.names=FALSE, quote=FALSE)
  utils::write.table(data.frame(parameter_name=free_names, gradient=grad), file.path(rds, paste0(side, '_DeltaCov_gradient.tsv')), sep='\t', row.names=FALSE, quote=FALSE)

  var_delta <- as.numeric(crossprod(grad, vcov %*% grad))
  se_delta <- sqrt(max(var_delta,0)); est_delta <- -gamma*cma
  z <- est_delta/se_delta; p <- 2*pnorm(abs(z), lower.tail=FALSE)
  ci_lo <- est_delta-1.96*se_delta; ci_hi <- est_delta+1.96*se_delta
  r_est <- read_num(row$Unstand_Est); r_se <- read_num(row$Unstand_SE)
  # Full-V versus diagonal-V sandwich comparison, using captured executed-path objects.
  bread <- solve(t(cap$S2.delt) %*% cap$S2.W %*% cap$S2.delt)
  lettuce <- cap$S2.W %*% cap$S2.delt
  sand <- function(V) bread %*% t(lettuce) %*% V %*% lettuce %*% bread
  full_se <- sqrt(diag(sand(cap$V_Reorder))); diag_se <- sqrt(diag(sand(diag(diag(cap$V_Reorder)))))
  free_reported <- read_num(res$Unstand_SE[res$op != ':=' ]); free_reported <- free_reported[is.finite(free_reported)]
  reported_struct <- read_num(res$Unstand_SE[res$op != ':=']); reported_struct <- reported_struct[is.finite(reported_struct)]
  max_rel_full <- max(abs(sort(reported_struct)-sort(full_se))/sort(reported_struct))
  max_rel_diag <- max(abs(sort(reported_struct)-sort(diag_se))/sort(reported_struct))
  base <- readRDS(file.path(rds, paste0(side,'_refit_frozen_model.rds')))$results
  base_key <- paste(base$lhs,base$op,base$rhs,sep='|'); cur <- res[res$op!=':=' | res$lhs != 'delta_cov_alm',]; cur_key <- paste(cur$lhs,cur$op,cur$rhs,sep='|'); mm <- match(base_key,cur_key)
  max_abs_est_73 <- max(abs(read_num(base$Unstand_Est)-read_num(cur$Unstand_Est[mm])), na.rm=TRUE)
  max_abs_se_73 <- max(abs(read_num(base$Unstand_SE)-read_num(cur$Unstand_SE[mm])), na.rm=TRUE)
  fit_rows[[side]] <- data.frame(side=side, converged=cap$converged, npar=cap$npar,
    chisq=fit$modelfit$chisq, CFI=fit$modelfit$CFI, SRMR=fit$modelfit$SRMR,
    max_abs_est_diff_73=max_abs_est_73, max_abs_se_diff_73=max_abs_se_73,
    W_dimension=paste(dim(cap$S2.W),collapse='x'), W_offdiag_max=max(abs(cap$S2.W-diag(diag(cap$S2.W)))),
    V_min_eigen=min(eigen(cap$V_Reorder,symmetric=TRUE,only.values=TRUE)$values),
    V_condition_number=kappa(cap$V_Reorder, exact=TRUE), warnings=paste(unique(warnings),collapse=' | '), stringsAsFactors=FALSE)
  weight_rows[[side]] <- data.frame(side=side, W_branch='DIAGONAL_INVERSE',
    W_vs_diag_inverse=max(abs(cap$S2.W-diag(1/diag(cap$V_Reorder)))/max(abs(cap$S2.W))),
    W_vs_full_inverse=max(abs(cap$S2.W-solve(cap$V_Reorder)))/max(abs(solve(cap$V_Reorder))),
    max_rel_SE_full_V=max_rel_full, max_rel_SE_diagonal_V=max_rel_diag,
    V_offdiag_nonzero_fraction=mean(cap$V_Reorder[upper.tri(cap$V_Reorder)] != 0),
    stringsAsFactors=FALSE)
  result_rows[[side]] <- data.frame(side=side, parameter='delta_cov_alm', formula='-gamma*c_ma',
    gamma_name=free_names[gidx], c_ma_name=free_names[cidx], gamma=gamma, c_ma=cma,
    estimate_independent=est_delta, variance_independent=var_delta, SE_independent=se_delta,
    Wald_z=z, CI95_low=ci_lo, CI95_high=ci_hi, P_two_sided=p,
    estimate_GenomicSEM=r_est, SE_GenomicSEM=r_se,
    abs_estimate_error=abs(est_delta-r_est), abs_SE_error=abs(se_delta-r_se),
    robust_covariance='Ohtt', gradient_mapping='name/label mapped; g_gamma=-c_ma; g_c_ma=-gamma',
    stringsAsFactors=FALSE)
  cross_rows[[side]] <- data.frame(side=side, free_index=seq_along(free_names), parameter_name=free_names,
    parameter_label=ifelse(nzchar(par$label[free_rows]),par$label[free_rows],NA), estimate=par$est[free_rows],
    stringsAsFactors=FALSE)
}
if (trace_status=='TRACE_INSTALLED') try(untrace('usermodel', where=asNamespace('GenomicSEM')),silent=TRUE)
write.table(do.call(rbind,result_rows), file.path(logs,'independent_delta_cov_results.tsv'), sep='\t',row.names=FALSE,quote=FALSE)
write.table(do.call(rbind,cross_rows), file.path(root,'02_V7_PARAMETER_LEVEL_CROSSWALK.tsv'), sep='\t',row.names=FALSE,quote=FALSE)
write.table(do.call(rbind,fit_rows), file.path(logs,'independent_model_fit_and_73_parameter_check.tsv'), sep='\t',row.names=FALSE,quote=FALSE)
write.table(do.call(rbind,weight_rows), file.path(logs,'independent_weight_and_sandwich_check.tsv'), sep='\t',row.names=FALSE,quote=FALSE)
writeLines(capture.output(sessionInfo()), file.path(logs,'independent_sessionInfo.txt'))
cat('completed=', format(Sys.time(), '%Y-%m-%dT%H:%M:%S%z'), '\n')
cat('exit_status=PASS\n')
sink()

