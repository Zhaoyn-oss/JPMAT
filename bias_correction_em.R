bias_correction_em <- function(beta, var_hat, tol = 1e-4, max_iter = 100) {
  
  ## ---------- Preprocess ----------
  #neither_na <- !(is.na(beta) | is.na(var_hat))
  #beta <- beta[neither_na]
  #nu0 <- var_hat[neither_na]
  nu0 <- var_hat
  n_tax <- length(beta)
  
  if(any(nu0[!is.na(nu0)] <= 0)) stop("Variance must be > 0")
  
  ## ---------- Stabilize variance ----------
  # floor
  eps <- quantile(nu0, 0.05, na.rm=TRUE)
  nu0 <- pmax(nu0, eps, na.rm=TRUE)
  
  ## ---------- Initial values ----------
  pi0 <- 0.75; pi1 <- 0.125; pi2 <- 0.125
  
  delta <- mean(beta, na.rm = TRUE)
  l1 <- quantile(beta, 0.1, na.rm = TRUE) - delta
  l2 <- quantile(beta, 0.9, na.rm = TRUE) - delta
  
  kappa1 <- var(beta, na.rm = TRUE) * 0.5
  kappa2 <- var(beta, na.rm = TRUE) * 0.5
  
  ## ---------- EM ----------
  iterNum <- 0
  epsilon <- Inf
  
  while (epsilon > tol && iterNum < max_iter) {
    
    ## ---- E-step (log-scale, stable) ----
    log_pdf0 <- dnorm(beta, delta, sqrt(nu0), log = TRUE)
    log_pdf1 <- dnorm(beta, delta + l1, sqrt(nu0 + kappa1), log = TRUE)
    log_pdf2 <- dnorm(beta, delta + l2, sqrt(nu0 + kappa2), log = TRUE)
    
    log_mat <- cbind(
      log(pi0) + log_pdf0,
      log(pi1) + log_pdf1,
      log(pi2) + log_pdf2
    )
    
    log_denom <- matrixStats::rowLogSumExps(log_mat)
    
    r0i <- exp(log_mat[,1] - log_denom)
    r1i <- exp(log_mat[,2] - log_denom)
    r2i <- exp(log_mat[,3] - log_denom)
    
    ## ---- M-step ----
    pi0_new <- mean(r0i, na.rm = TRUE)
    pi1_new <- mean(r1i, na.rm = TRUE)
    pi2_new <- mean(r2i, na.rm = TRUE)
    
    ## ---- weight clipping（关键稳定）----
    w0 <- pmin(1/nu0, quantile(1/nu0, 0.95, na.rm = TRUE))
    w1 <- pmin(1/(nu0 + kappa1), quantile(1/(nu0 + kappa1), 0.95, na.rm = TRUE), na.rm = TRUE)
    w2 <- pmin(1/(nu0 + kappa2), quantile(1/(nu0 + kappa2), 0.95, na.rm = TRUE), na.rm = TRUE)
    
    delta_new <- sum(r0i*w0*beta +
                       r1i*w1*(beta - l1) +
                       r2i*w2*(beta - l2), na.rm = TRUE) /
      sum(r0i*w0 + r1i*w1 + r2i*w2, na.rm = TRUE)
    
    l1_new <- min(
      sum(r1i*w1*(beta - delta_new), na.rm = TRUE) / sum(r1i*w1, na.rm = TRUE),
      0
    )
    if (is.na(l1_new)) l1_new <- 0
    
    l2_new <- max(
      sum(r2i*w2*(beta - delta_new), na.rm = TRUE) / sum(r2i*w2, na.rm = TRUE),
      0
    )
    if (is.na(l2_new)) l2_new <- 0
    
    ## ---- Update kappa----
    obj_kappa1 <- function(x){
      x <- max(x, 1e-6, na.rm = TRUE)
      -sum(r1i * dnorm(beta, delta_new + l1_new,
                       sqrt(nu0 + x), log=TRUE), na.rm = TRUE)
    }
    
    obj_kappa2 <- function(x){
      x <- max(x, 1e-6, na.rm = TRUE)
      -sum(r2i * dnorm(beta, delta_new + l2_new,
                       sqrt(nu0 + x), log=TRUE), na.rm = TRUE)
    }
    
    kappa1_new <- nloptr::neldermead(x0 = kappa1, fn = obj_kappa1)$par
    kappa2_new <- nloptr::neldermead(x0 = kappa2, fn = obj_kappa2)$par
    
    ## constrain kappa（关键）
    kappa1_new <- min(max(kappa1_new, 1e-6, na.rm = TRUE), quantile(nu0, 0.9, na.rm = TRUE))
    kappa2_new <- min(max(kappa2_new, 1e-6, na.rm = TRUE), quantile(nu0, 0.9, na.rm = TRUE))
    
    ## ---- convergence ----
    epsilon <- sqrt((pi0_new-pi0)^2 + (pi1_new-pi1)^2 + (pi2_new-pi2)^2 +
                      (delta_new-delta)^2 + (l1_new-l1)^2 + (l2_new-l2)^2 +
                      (kappa1_new-kappa1)^2 + (kappa2_new-kappa2)^2)
    
    ## update
    pi0 <- pi0_new; pi1 <- pi1_new; pi2 <- pi2_new
    delta <- delta_new
    l1 <- l1_new; l2 <- l2_new
    kappa1 <- kappa1_new; kappa2 <- kappa2_new
    
    iterNum <- iterNum + 1
  }
  
  ## ---------- WLS estimator（稳定版） ----------
  nu <- nu0
  nu <- nu + r1i*kappa1 + r2i*kappa2
  
  w <- 1 / nu
  w <- pmin(w, quantile(w, 0.95, na.rm = TRUE))  # clip
  
  delta_wls <- sum(w * beta, na.rm = TRUE) / sum(w, na.rm = TRUE)
  var_delta <- 1 / sum(w, na.rm = TRUE)
  
  ## ---------- Output ----------
  list(
    delta_em = delta,
    delta_wls = delta_wls,
    var_delta = var_delta,
    w=w 
  )
}