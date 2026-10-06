 estimate_alpha_beta_theta <- function(R, G, X, fix_formula, rand_formula=NULL, 
                                       alpha=0.05, p_adj_method="BH", emp_null_calibration=TRUE,
                                       verbose=TRUE, iter_control=list()) { 
   X <- as.data.frame(X, stringsAsFactors = FALSE) 
   if(!is.null(R) & !is.null(G)){ 
     #====Initialization for theta========   
     theta_G <- log(rowSums(G, na.rm=TRUE)+1); theta_G <- theta_G-mean(theta_G, na.rm=TRUE)
     theta_R <- log(rowSums(R, na.rm=TRUE)+1); theta_R <- theta_R-mean(theta_R, na.rm=TRUE) 
     theta0 <- c(theta_G, theta_R)  
     #====reparameter for estimate======
     rownames(G) <- paste0(rownames(G), "_G") 
     rownames(R) <- paste0(rownames(R), "_R") 
     y <- rbind(G, R) 
     Z <- rbind(X, X); rownames(Z) <- rownames(y)
     isR <- factor(c(rep(0, nrow(G)), rep(1, nrow(R))))
     Z <- data.frame(Z, isR=isR)
     formula <- paste0("(", fix_formula, ")*(isR)")
     
     if(is.null(rand_formula)){
       rand_formula <- "(1|paired_sample)"
     }else{  
       x <- gsub("^\\s*\\(|\\)\\s*$", "", rand_formula)
       parts <- strsplit(x, "\\|")[[1]] 
       rand_var   <- trimws(parts[1])
       group_var  <- trimws(parts[2])
       rand_formula <- paste0("(isR*", "(", rand_var,")",  "|", group_var, ")",  "+", "(1|paired_sample)") 
     }  
     metadata <- Z
     paired_sample <- rep(paste0("paired_sample", 1:nrow(X)), times = 2)
     metadata$paired_sample <- factor(paired_sample)
     options(na.action = "na.pass") # Keep NA's in rows of x
     Xmodel <- model.matrix(formula(paste0("~", formula)), data = Z)
     options(na.action = "na.omit") # Switch it back   
     
   } else {
     data <- if(!is.null(G)) G else R 
     y <- data
     rownames(X) <- rownames(y)
     formula <- fix_formula
     rand_formula <- rand_formula
     metadata <- X
     options(na.action = "na.pass")  
     Xmodel <- model.matrix(formula(paste0("~", formula)), data = X)
     options(na.action = "na.omit") 
   }
   
   if (ncol(y) < 50) {
     warn_txt <- sprintf(paste0("The number of taxa used for estimating ",
                                "sample-specific biases is: ",
                                ncol(y),
                                "\nA large number of taxa (>50) is required ",
                                "for the consistent estimation of biases"))
     warning(warn_txt, call. = FALSE)
   }
   
   ## ======== First MLE for gamma======
   para1 <- iter_nb(y=y, 
                    x=Xmodel, 
                    metadata=metadata, 
                    fix_formula=formula, 
                    rand_formula=rand_formula,
                    theta0=theta0,
                    max_iter = iter_control$max_MLE_iter,
                    tol = iter_control$tol_MLE,
                    verbose=verbose) 
   gamma1 <- para1$gamma
   var_hat1 <- para1$var_hat
   var_hat1[which(var_hat1<=0)] <- NA 
   
   ## ===== EM correction =========
   fun_list <- list(bias_correction_em) 
   bias1 <- foreach(i = seq_len(ncol(gamma1))) %dorng% {
     output <- fun_list[[1]](beta = gamma1[, i],
                             var_hat = var_hat1[, i],
                             tol = iter_control$tol_EM,
                             max_iter = iter_control$max_EM_iter) 
   }
   delta_em <- sapply(bias1, function(bias){bias$delta_em})
   delta_wls <- sapply(bias1, function(bias){bias$delta_wls})
   var_delta <- sapply(bias1, function(bias){bias$var_delta}) 
   names(delta_em) <- names(delta_wls) <- names(var_delta) <- colnames(Xmodel)
    
   # Obtain the final estimates and sample-specific biases
   gamma_hat <- t(t(gamma1)-delta_em)
   theta <- para1$theta
   samp_frac <- theta+as.numeric(Xmodel %*% delta_em)  
     
   var_hat <- t(t(var_hat1)+var_delta)+2*t(t(sqrt(var_hat1))*sqrt(var_delta)) 
   var_hat[is.na(gamma_hat)] <- NA 
   vcov_hat <- para1$vcov_hat 
   if(emp_null_calibration){
     #==robust empirical-null scale calibration====
     W0 <- gamma_hat/sqrt(var_hat) 
     p0 <- 2 * pnorm(abs(W0),lower.tail = FALSE)
     p0[is.na(p0)] <- 1 
     q0 <- apply(p0, 2, function(x) p.adjust(x, method = p_adj_method))  
     sigma0 <- numeric(ncol(W0))
     cut_grid <- seq(0, 0.30, by = 0.05)
     for(i in 1:ncol(W0)){
       sigmas <- sapply(cut_grid, function(cut) {
         id <- which(q0[, i] > cut & is.finite(W0[, i]))
         if (length(id) < 10) return(NA_real_) 
         s <- mad(W0[id, i],center = median(W0[id, i], na.rm = TRUE),constant = 1.4826,na.rm = TRUE)  
         max(s, 1)
       })
       n_keep <- sapply(cut_grid, function(cut) {
         sum(q0[, i] > cut & is.finite(W0[, i]))
       })
       cutoff <- choose_cutoff(cut_grid, sigmas, n_keep, eps = 0.02, K = 2)
       null_id <- which(q0[, i] > cutoff & is.finite(W0[, i]))
       if (length(null_id) < 5) null_id <- which(is.finite(W0[, i])) 
       x <- W0[null_id, i]
       s <- mad(x, center = median(x, na.rm = TRUE), constant = 1.4826, na.rm = TRUE)
       sigma0[i] <- max(s, 1, na.rm = TRUE)
     }
     names(sigma0) <- colnames(Xmodel)
     sigma1 <- sigma0
     var_hat <- sweep(var_hat, 2, sigma0^2,FUN = "*") 
   } 
   se_hat <- sqrt(var_hat) 
   #=======primary results======= 
   W <- gamma_hat/se_hat 
   p_hat <- 2 * pnorm(abs(W),lower.tail = FALSE)
   colnames(p_hat) <- colnames(W)
   p_hat[is.na(p_hat)] <- 1 
   q_hat <- apply(p_hat, 2, function(x) p.adjust(x, method = p_adj_method))
   diff_abn <- (q_hat <= alpha & !is.na(q_hat))
   gamma_prim <- data.frame(gamma_hat, check.names = FALSE)
   se_prim <- data.frame(se_hat, check.names = FALSE)
   W_prim <- data.frame(W, check.names = FALSE)
   p_prim <- data.frame(p_hat, check.names = FALSE) 
   q_prim <- data.frame(q_hat, check.names = FALSE)
   diff_prim <- data.frame(diff_abn, check.names = FALSE)
   colnames(gamma_prim) <- paste0("lfc_", colnames(gamma_prim))
   colnames(se_prim) <- paste0("se_", colnames(se_hat))
   colnames(W_prim) <- paste0("W_", colnames(W))
   colnames(p_prim) <- paste0("p_", colnames(p_hat)) 
   colnames(q_prim) <- paste0("q_", colnames(q_hat))
   colnames(diff_prim) <- paste0("diff_", colnames(diff_abn))
   
   #========MTX differential expression effect======
   if(!is.null(R) & !is.null(G)){ 
     #-------alpha+beta test----
     var_est <- colnames(gamma_hat)  
     main_vars <- var_est[!grepl(":", var_est)& !grepl("Intercept", var_est)  & !grepl("^isR", var_est)]
     pairs <- lapply(main_vars, function(v) {
       beta_name <- paste0(v, ":isR1")
       if (beta_name %in% var_est) {
         c(v, beta_name)
       } else {
         NULL
       }
     })
     pairs <- pairs[!sapply(pairs, is.null)]
     names(pairs) <- main_vars
     
     global_est <- matrix(NA, nrow(gamma_hat), length(pairs))
     global_var <- matrix(NA, nrow(gamma_hat), length(pairs))
     for (i in seq_along(pairs)) {
       v <- names(pairs)[i]
       alpha <- pairs[[v]][1]
       beta  <- pairs[[v]][2] 
       for(j in 1:nrow(gamma_hat)){
         theta_hat <- gamma_hat[j, c(alpha, beta)] 
         vcov_sub <- vcov_hat[[j]][c(alpha, beta), c(alpha, beta)]
         if (any(is.na(theta_hat))) {
           global_est[j, i] <- NA
         } else {
           global_est[j, i] <- sum(theta_hat)
         }
         
         if(!is.null(vcov_sub)){
           v1 <- vcov_sub[1,1]
           v2 <- vcov_sub[2,2]
           c12 <- vcov_sub[1,2]
           V <- v1 + v2 + 2*c12
           if (is.na(V) || V <= 1e-12) {
             global_var[j, i] <- NA 
           } else {
             global_var[j, i] <- V
           } 
         }else{
           global_var[j,i] <- NA
         }  
       }  
     }   
     
     if(emp_null_calibration){
       #==robust empirical-null scale calibration====
       W0 <- global_est/sqrt(global_var)
       p0 <- 2 * pnorm(abs(W0),lower.tail = FALSE)
       p0[is.na(p0)] <- 1 
       q0 <- apply(p0, 2, function(x) p.adjust(x, method = p_adj_method))  
       sigma0 <- numeric(ncol(W0))
       cut_grid <- seq(0, 0.30, by = 0.05)
       for(i in 1:ncol(W0)){
         sigmas <- sapply(cut_grid, function(cut) {
           id <- which(q0[, i] > cut & is.finite(W0[, i]))
           if (length(id) < 10) return(NA_real_) 
           s <- mad(W0[id, i],center = median(W0[id, i], na.rm = TRUE),constant = 1.4826,na.rm = TRUE) 
           max(s, 1)
         })
         n_keep <- sapply(cut_grid, function(cut) {
           sum(q0[, i] > cut & is.finite(W0[, i]))
         })
         cutoff <- choose_cutoff(cut_grid, sigmas,n_keep, eps = 0.02, K = 2)
         null_id <- which(q0[, i] > cutoff & is.finite(W0[, i]))
         if (length(null_id) < 5) null_id <- which(is.finite(W0[, i])) 
         x <- W0[null_id, i]
         s <- mad(x, center = median(x, na.rm = TRUE), constant = 1.4826, na.rm = TRUE)
         sigma0[i] <- max(s, 1, na.rm = TRUE)
       }
       global_var <- sweep(global_var, 2, sigma0^2, FUN = "*")
     } 
     global_se <- sqrt(global_var)
     global_w <- global_est/global_se 
     global_p_hat <- 2 * pnorm(abs(global_w), lower.tail = FALSE)
     global_q_hat <- apply(global_p_hat, 2, function(x) p.adjust(x, method = p_adj_method))
     colnames(global_est) <- paste0("lfc_sum_", names(pairs))
     colnames(global_se) <- paste0("se_sum_", names(pairs))
     colnames(global_w) <- paste0("W_sum_", names(pairs))
     colnames(global_p_hat) <- paste0("p_sum_", names(pairs))
     colnames(global_q_hat) <- paste0("q_sum_", names(pairs))
     
     res <- do.call("cbind", list(data.frame(taxon = colnames(y)),
                                  gamma_prim, global_est, se_prim, global_se, W_prim, global_w,
                                  p_prim, global_p_hat, q_prim, global_q_hat, diff_prim)) 
     rownames(res) <- NULL  
   }else{
     res <- do.call("cbind", list(data.frame(taxon = colnames(y)),
                                  gamma_prim, se_prim, W_prim, 
                                  p_prim, q_prim, diff_prim)) 
     rownames(res) <- NULL  
   } 
   
   if(emp_null_calibration){
     out <- list(res=res,  
                 samp_frac = samp_frac,
                 delta_em = delta_em, 
                 var_delta=var_delta,
                 sigma=sigma1)
   }else{
     out <- list(res=res,  
                 samp_frac = samp_frac,
                 delta_em = delta_em, 
                 var_delta=var_delta)
   }
   
   return(out)  
 }
 