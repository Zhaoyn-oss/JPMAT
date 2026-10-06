# data pre-processing
data_check <- function(G, R, X, formula){ 
  #---- Check sample id-----
  X <- as.data.frame(X)
  sample.X <- rownames(X) 
  if (is.null(sample.X)) {
    stop("X must have rownames as sample IDs.")
  }
  if (!is.null(G)) {
    G <- as.data.frame(G)
    sample.G <- rownames(G)
    if (is.null(sample.G)) {
      stop("G must have rownames as sample IDs.")
    }
  }
  if (!is.null(R)) {
    R <- as.data.frame(R)
    sample.R <- rownames(R)
    if (is.null(sample.R)) {
      stop("R must have rownames as sample IDs.")
    }
  }
  sample.list <- list(sample.X)
  if (!is.null(G)) sample.list <- c(sample.list, list(sample.G))
  if (!is.null(R)) sample.list <- c(sample.list, list(sample.R)) 
  sample.common <- Reduce(intersect, sample.list)
  if (length(sample.common) == 0) {
    stop("No common samples found among input objects.")
  }
  # Keep stable order
  X <- X[sample.common, , drop = FALSE]
  if (!is.null(G)) G <- G[sample.common, , drop = FALSE]
  if (!is.null(R)) R <- R[sample.common, , drop = FALSE]
  
  # Drop unused factor levels
  X[] <- lapply(X, function(x) {
    if (is.factor(x)) droplevels(x) else x
  })
   
  # Check formula variables
  if (!is.null(formula)) {  
    vars <- all.vars(formula)
    missing_vars <- vars[!vars %in% colnames(X)] 
    if (length(missing_vars) > 0) {
      stop("The following variables specified are not in X: ",
           paste(missing_vars, collapse = ", "))
    }
  }
   
  # Align features (G & R)
  if (!is.null(G) && !is.null(R)) { 
    feature.G <- colnames(G)
    feature.R <- colnames(R) 
    if (is.null(feature.G) || is.null(feature.R)) {
      stop("Both G and R must have colnames as feature IDs.")
    } 
    feature.common <- intersect(feature.G, feature.R) 
    if (length(feature.common) == 0) {
      stop("No common features found between G and R.")
    } 
    G <- G[, feature.common, drop = FALSE]
    R <- R[, feature.common, drop = FALSE]
  } 
  return(list(G = G, R = R, X = X))
}

# Structural zero helper
handle_struc_zero <- function(mat, X, struc_zero, group=NULL) {
  if (!struc_zero) {
    return(list(zero_ind = NULL,
                tax_keep = seq_len(ncol(mat))))
  }else{
    zero_ind <- get_struc_zero(
      data = mat,
      meta_data = X,
      group = group 
    )
    
    tax_idx <- apply(zero_ind[, -1], 1,
                     function(x) all(x == FALSE))
    
    return(list(zero_ind = zero_ind, 
                tax_keep = which(tax_idx)))
  } 
}

# Prevalence filtering helper
filter_core <- function(mat, X, prv_cut, lib_cut, tax_keep = NULL) { 
  core1 <- data_core(
    data = mat,
    meta_data = X,
    prv_cut = prv_cut,
    lib_cut = lib_cut,
    tax_keep = NULL,
    samp_keep = NULL
  )
  
  mat1 <- t(core1$feature_table)
  samp_keep <- rownames(mat1)
  
  core2 <- data_core(
    data = mat,
    meta_data = X,
    prv_cut = prv_cut,
    lib_cut = lib_cut,
    tax_keep = tax_keep,
    samp_keep = samp_keep
  )
  
  list(
    mat1 = mat1,
    mat2 = t(core2$feature_table),
    meta = core2$meta_data
  )
}

# Identify structural zeros
get_struc_zero <- function(data, meta_data, group) {
  feature_table <- data
  tax_name <- colnames(data)
  group_data <- factor(meta_data[, group])
  present_table <- as.matrix(feature_table)
  present_table[is.na(present_table)] <- 0
  present_table[present_table != 0] <- 1
  n_tax <- ncol(feature_table)
  n_group <- nlevels(group_data)
  
  p_hat <- matrix(NA, nrow = n_tax, ncol = n_group)
  rownames(p_hat) <- colnames(feature_table)
  colnames(p_hat) <- levels(group_data)
  for (i in seq_len(n_tax)) {
    p_hat[i, ] <- tapply(present_table[, i], group_data,
                        function(x) mean(x, na.rm = TRUE))
  } 
  output <- (p_hat == 0)  
  output <- cbind(tax_name, output)
  colnames(output) <- c("taxon",
                       paste0("structural_zero (", group,
                              " = ", colnames(output)[-1], ")"))
  output <- data.frame(output, check.names = FALSE, row.names = NULL)
  output[, -1] <- apply(output[, -1], 2, as.logical)
  return(output)
}

# Filter data by prevalence and library size
data_core <- function(data, meta_data, prv_cut, lib_cut,
                      tax_keep = NULL, samp_keep = NULL) {
  feature_table <- data
  
  # Discard taxa with prevalences < prv_cut
  if (is.null(tax_keep)) {
    prevalence <- apply(feature_table, 1, function(x) sum(x != 0, na.rm = TRUE)/length(x[!is.na(x)]))
    tax_keep <- which(prevalence >= prv_cut)
  }else if (length(tax_keep) == 0) {
    stop("All taxa contain structural zeros", call. = FALSE)
  } else {
    # Discard taxa with structural zeros
    feature_table <- feature_table[tax_keep, , drop = FALSE]
    prevalence <- apply(feature_table, 1, function(x) sum(x != 0, na.rm = TRUE)/length(x[!is.na(x)]))
    tax_keep <- which(prevalence >= prv_cut)
  }
  
  if (length(tax_keep) > 0) {
    feature_table <- feature_table[tax_keep, , drop = FALSE]
  } else {
    stop("No taxa remain under the current cutoff", call. = FALSE)
  }
  
  # Discard samples with library sizes < lib_cut
  if (is.null(samp_keep)) {
    lib_size <- colSums(feature_table, na.rm = TRUE)
    samp_keep <- which(lib_size >= lib_cut)
  }
  if (length(samp_keep) > 0){
    feature_table <- feature_table[, samp_keep, drop = FALSE]
    meta_data <- meta_data[samp_keep, , drop = FALSE]
  } else {
    stop("No samples remain under the current cutoff", call. = FALSE)
  }
  
  output <- list(feature_table = feature_table,
                meta_data = meta_data,
                tax_keep = tax_keep,
                samp_keep = samp_keep)
  return(output)
}
 

# NB regression for each taxon
fit_taxon_each <- function(tformula, df, theta, start=NULL){ 
  tformula <- as.formula(tformula)
  fit <- tryCatch({
    glmmTMB::glmmTMB(tformula, dispformula = ~ isR,
                     family = glmmTMB::nbinom2(), 
                     data = df, offset =theta, start = start)
  }, 
  error=function(e) {
    message(conditionMessage(e)) 
    NULL
  })
  fit
}


orthogonalize_to_X <- function(v, X) {
  # remove projection of v onto col(X)
  fit <- lm.fit(x = X, y = v)
  res <- fit$residuals
  as.numeric(res)
}

choose_cutoff <- function(cut_grid, sigmas, n_keep, eps = 0.02, K = 2) { 
  ok <- is.finite(cut_grid) & is.finite(sigmas) 
  cut_grid <- cut_grid[ok]
  sigmas <- sigmas[ok] 
  n_keep <- n_keep[ok]
  if (length(sigmas) < 2) return(cut_grid[1])  
  # Relative change between adjacent sigma estimates
  rel_change <- abs(diff(sigmas) / sigmas[-length(sigmas)]) 
  change_n <- abs(diff(n_keep))
  # Only count comparisons where the retained taxon set actually changes
  informative <- change_n >= 1
  for (i in seq_along(cut_grid)) { 
    ids <- i:length(rel_change) 
    # Only informative transitions
    ids <- ids[informative[ids]] 
    if (length(ids) >= K) { 
      ids <- ids[1:K] 
      if (all(rel_change[ids] <= eps, na.rm = TRUE)){
        return(cut_grid[i]) 
      } 
    }
  } 
  # conservative fallback
  cut_grid[1]
}
