#' Joint Analysis of Paired MGX and MTX Count Data
#'
#' Fits a joint model for paired metagenomic (MGX) and
#' metatranscriptomic (MTX) counts to assess differential genomic
#' abundance, transcriptional activity, and transcript expression.
#'
#' @param R Numeric matrix of raw MTX counts, with samples in rows
#'   and features in columns.
#' @param G Numeric matrix of raw MGX counts, with samples in rows
#'   and features in columns. Feature and sample identifiers should
#'   correspond to those in \code{R}.
#' @param metadata Data frame containing sample metadata and variables
#'   used in the model formulas.
#' @param fix_formula Specification of the fixed effects.
#' @param rand_formula Specification of the random effects.
#'   Default is \code{NULL}.
#' @param struc_zero Logical; whether to identify structural zeros.
#'   Default is \code{TRUE}.
#' @param group Grouping variable used for structural-zero identification.
#'   Default is \code{NULL}.
#' @param lib_cut Library-size filtering threshold. Default is \code{0}.
#' @param prv_cut Minimum feature prevalence threshold.
#'   Default is \code{0.05}.
#' @param alpha Significance level. Default is \code{0.05}.
#' @param p_adj_method Method for adjusting p values for multiple testing.
#'   Default is \code{"BH"}.
#' @param n_cl Number of workers used for parallel computation.
#'   Default is \code{4}.
#' @param emp_null_calibration Logical; whether to perform empirical
#'   null calibration. Default is \code{TRUE}.
#' @param verbose Logical; whether to display progress messages.
#'   Default is \code{FALSE}.
#' @param iter_control List of iteration controls:
#'   \describe{
#'     \item{max_MLE_iter}{Maximum number of MLE iterations.}
#'     \item{tol_MLE}{Convergence tolerance for MLE iterations.}
#'     \item{max_EM_iter}{Maximum number of EM iterations.}
#'     \item{tol_EM}{Convergence tolerance for EM iterations.}
#'   }
#'
#' @return The fitted analysis results. See the returned object's
#'   components for model estimates and statistical tests.
#' @importFrom foreach %dopar%
#'
#' @export
JPMAT <- function(R, G, metadata, fix_formula, rand_formula=NULL,
                  struc_zero=TRUE, group=NULL,
                  lib_cut=0, prv_cut = 0.05, alpha=0.05,
                  p_adj_method = "BH", n_cl=4,
                  emp_null_calibration=TRUE, verbose=FALSE,
                  iter_control=list(max_MLE_iter=100, tol_MLE=0.01,
                                    max_EM_iter=100, tol_EM=1e-5)){

  #----0. Parallel setup------
  if (n_cl > 1) {
    cl <- parallel::makeCluster(n_cl)
    doParallel::registerDoParallel(cl)
  } else {
    foreach::registerDoSEQ()
  }
  #---- 1 data check -----
  formula <- if(is.null(rand_formula)) as.formula(paste0("~", fix_formula)) else as.formula(paste0("~", fix_formula, "+", rand_formula))
  check <- data_check(G=G, R=R, X=metadata, formula=formula)
  G <- check$G
  R <- check$R
  metadata <- check$X

  if (!is.null(G) && !is.null(R)) {
    # structure zero for MGX and MTX
    sz_G <- handle_struc_zero(mat=G, X=metadata, struc_zero=struc_zero, group=group)
    sz_R <- handle_struc_zero(mat=R, X=metadata, struc_zero=struc_zero, group=group)
    tax_keep <- intersect(sz_G$tax_keep, sz_R$tax_keep)

    # prevalence and library size check for MGX and MTX
    f_G <- filter_core(t(G), X=metadata, prv_cut=prv_cut, lib_cut=lib_cut, tax_keep=tax_keep)
    f_R <- filter_core(t(R), X=metadata, prv_cut=prv_cut, lib_cut=lib_cut, tax_keep=tax_keep)
    taxon_mgx <- colnames(f_G$mat2)
    taxon_mtx <- colnames(f_R$mat2)
    common_taxa <- intersect(taxon_mgx, taxon_mtx)

    sample_mgx <- rownames(f_G$meta)
    sample_mtx <- rownames(f_R$meta)
    common_sample <- union(sample_mgx, sample_mtx)

    G <- G[common_sample, common_taxa, drop = FALSE]
    R <- R[common_sample, common_taxa, drop = FALSE]
    Xf <- metadata[common_sample, , drop = FALSE]

    zero_ind <- list(zero_ind_G = sz_G$zero_ind,
                     zero_ind_R = sz_R$zero_ind)

    output <- estimate_alpha_beta_theta(
      R = R,
      G = G,
      X = Xf,
      fix_formula = fix_formula,
      rand_formula = rand_formula,
      alpha = alpha,
      p_adj_method = p_adj_method,
      emp_null_calibration=emp_null_calibration,
      verbose=verbose,
      iter_control = iter_control
    )

  } else {
    mat <- if(!is.null(R)) R else G
    # prevalence and library size check for MGX or MTX
    f_mat <- data_core(t(mat), meta_data=metadata, prv_cut=prv_cut, lib_cut=lib_cut)
    common_taxon <- rownames(f_mat$feature_table)
    common_sample <- rownames(f_mat$meta_data)

    mat <- mat[common_sample, common_taxon, drop = FALSE]
    Xf <- metadata[common_sample, , drop = FALSE]

    zero_ind <- sz_mat$zero_ind

    output <- estimate_alpha_beta_theta(
      R = if (!is.null(R)) mat else NULL,
      G = if (!is.null(G)) mat else NULL,
      X = Xf,
      fix_formula = fix_formula,
      rand_formula = rand_formula,
      alpha = alpha,
      p_adj_method = p_adj_method,
      iter_control = iter_control
    )
  }

  out <- list(
    samp_frac = output$samp_frac,
    delta_em = output$delta_em,
    var_delta = output$var_delta,
    res = output$res,
    vcov_hat=output$vcov_hat,
    sigma0 =output$sigma
  )

  if (n_cl > 1) {
    parallel::stopCluster(cl)
  }
  return(out)
}







