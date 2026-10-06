iter_nb <- function(y, x, metadata, fix_formula, rand_formula=NULL, 
                    theta0, max_iter = 50, tol = 0.01, verbose=FALSE){
  n <- nrow(y)
  d <- ncol(y)
  taxa <- colnames(y) 
  samp_id <- rownames(y)
  fix_eff <- colnames(x)
  p <- length(fix_eff)  
  
  if(!is.null(rand_formula)) {
    tformula <- paste0("y_crt ~ ", fix_formula, "+", rand_formula) 
  } else{
    tformula <- paste0("y_crt ~ ", fix_formula)
  }
  
  data_list <- lapply(seq_len(d), function(j){ 
    data.frame(metadata, y_crt=y[, j])
  })
  fit_list <- list(fit_taxon_each) 
  ## ---------- Initialization ---------- 
  theta <- orthogonalize_to_X(theta0, x)
  gamma <- matrix(0, d, p) 
  epsilon <- Inf
  iterNum <- 0
  
  ## ---------- Iterative update ---------- 
  while (epsilon > tol && iterNum < max_iter) {  
    #=== 1. Refit feature-wise models using current theta====
    fits <- foreach(j = seq_len(d)) %dopar% {
      fit <- fit_list[[1]](tformula=tformula, 
                           df=data_list[[j]], 
                           theta=theta) }
    
    ## Extract gamma  
    gamma_new <- matrix(NA, d, p) 
    colnames(gamma_new) <- fix_eff
    rownames(gamma_new) <- taxa
    for(j in seq_len(d)){ 
      fit <- fits[[j]] 
      if(inherits(fit, "glmmTMB")){ 
        coef_i <- fixef(fit)$cond 
        gamma_new[j, match(names(coef_i), fix_eff) ] <- coef_i  
      }
    } 
    
    ## Compute mu 
    mu_mat <- sapply(fits, function(fit) {
      y_crt_hat_i <- rep(NA, n)
      if(inherits(fit, "glmmTMB")){
        fitted_i <- predict(fit, re.form = NULL, type="response") 
        names(fitted_i) <- rownames(fit$frame)
      } else {
        fitted_i <- rep(NA, n)
        names(fitted_i) <- samp_id
      } 
      y_crt_hat_i[match(names(fitted_i), samp_id)] <- fitted_i
      return(y_crt_hat_i)
    })
    
    ## Extract observation-specific phi for each taxon
    phi_mat <- sapply(fits, function(fit) {
      phi_i <- rep(NA_real_, n)
      if (inherits(fit, "glmmTMB")) {
        disp_i <- predict(fit, type = "disp")
        names(disp_i) <- rownames(fit$frame)
        phi_i[match(names(disp_i), samp_id)] <- as.numeric(disp_i)
      }
      phi_i
    })
    
    #==2. update theta ===== 
    mu_inner <- mu_mat 
    theta_new <- theta
    
    eps_inner <- Inf
    max_inner <- 100
    inner_iter <- 0
    while (eps_inner > 1e-5 && inner_iter < max_inner) { 
      inner_iter <- inner_iter+1 
      denom <- 1 + mu_inner / phi_mat
      score <- rowSums((y - mu_inner) / denom, na.rm = TRUE)
      score[!is.finite(score)] <- 0
      
      fisher <- rowSums(mu_inner/denom, na.rm = TRUE)
      fisher[!is.finite(fisher)] <- NA_real_
      fisher <- pmax(fisher, 1e-8)
      
      # contrained Newton step, contraint X'theta=0
      # Newton step : delta=H^{-1} [U-X lambda]
      # lambda=(X'H^{-1} X)^{-1} X'H^{-1} U 
      fisher_inv <- 1 / fisher
      HinvC <- x * fisher_inv
      CtHinvC <- crossprod(x,HinvC)
      CtHinvU <- crossprod(x,fisher_inv * score)
      lambda <- tryCatch(
        {
          qr.solve(CtHinvC,CtHinvU)
        },
        error = function(e) {
          MASS::ginv(CtHinvC) %*% CtHinvU
        }
      ) 
      raw_step <- fisher_inv * (score -drop(x %*% lambda))
      
      # scalar damping/step control
      q95_step <- as.numeric(quantile(abs(raw_step),probs = 0.95,na.rm = TRUE))
      step_cap <- min(2, max(0.5, q95_step))
      max_step <- max(abs(raw_step),na.rm = TRUE)
      if (is.finite(max_step) && max_step > step_cap) {
        scale_factor <- step_cap / max_step
      } else {
        scale_factor <- 1
      }
      step_use <- scale_factor * raw_step
      
      # update theta
      theta_candidate <- theta_new + step_use
      theta_candidate <- orthogonalize_to_X(theta_candidate,x)
      
      actual_step <- theta_candidate - theta_new
      theta_new <- theta_candidate
      mu_inner <- sweep(mu_inner,1,exp(actual_step),`*`)
      eps_inner <- sqrt(mean(actual_step^2,na.rm = TRUE))
    } 
    #==3. outer convergence diagnosistic======   
    delta_gamma <- gamma_new - gamma 
    delta_theta <- theta_new-theta
    eps_gamma <- sqrt(mean(delta_gamma^2, na.rm=TRUE))
    eps_theta <- sqrt(mean(delta_theta^2, na.rm=TRUE))
    epsilon <- max(eps_gamma, eps_theta)
    
    #==update parameters===
    iterNum <- iterNum + 1
    gamma <- gamma_new
    theta <- theta_new   
    if (verbose) {
      message(
        "Iteration: ", iterNum, 
        "  epsilon = ", signif(epsilon,5) 
      )
    }  
  }  
  
  vcov_hat <- lapply(fits, function(fit) { 
    Sigma_hat_i <- matrix(NA_real_, nrow = p, ncol = p)
    dimnames(Sigma_hat_i) <- list(fix_eff, fix_eff) 
    if (inherits(fit, "glmmTMB")) { 
      vcov_hat_i <- as.matrix(vcov(fit)$cond) 
      r_idx <- match(rownames(vcov_hat_i), fix_eff)
      c_idx <- match(colnames(vcov_hat_i), fix_eff) 
      Sigma_hat_i[r_idx, c_idx] <- vcov_hat_i
    }
  })
  var_hat <- t(sapply(vcov_hat, function(Sigma_hat_i){
    if(!is.null(Sigma_hat_i)){
      diag(Sigma_hat_i)
    } else {
      rep(NA, p)
    }
  })) 
  
  ## ---------- Output ---------- 
  rownames(gamma) <- taxa
  colnames(gamma) <- fix_eff
  rownames(var_hat) <- taxa
  colnames(var_hat) <- fix_eff
  names(theta) <- rownames(y)
  names(vcov_hat) <- taxa 
  
  para1 <- list(
    gamma = gamma,
    theta = theta,
    vcov_hat = vcov_hat,
    var_hat = var_hat
  ) 
  
  return(para1)
}
