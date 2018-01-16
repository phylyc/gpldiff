#' Compute the Bayesian false discovery rate from posterior probability.
#'
#' @param v  sorted (by descending order) posterior probablity that a gene 
#'           is in the positive class
#' @return false discovery rate incurred by declaring each gene as significant 
#'         in order
bayesian_fdr <- function(v) {
	cumsum(1 - v) / 1:length(v)
}
