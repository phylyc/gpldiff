#' Root mean squared error
#'
#' Calcuate the root mean squared error
#'
#' @param x  first set of values
#' @param y  second set of values
#' @return \code{numeric} value
#' @export
#'
rmse <- function(x, y) {
	sqrt(mean((x - y)^2))
}

coverage <- function (model, data, ...) {
	UseMethod("coverage", model)
}

#' Coverage probability
#'
#' Calculate the coverage probability by a normally distributed confidnece
#' interval on observed data.
#'
#' @param model  numeric vector of mean and standard error
#' @param data  observed data
#' @param level  confidence level
#' @return empirical coverage probability
#' @export
#'
coverage.default <- function(model, data, level=0.95) {
	if (length(model) != 2 || !is.list(model) || is.null(model$mean) || is.null(model$sd)) {
		stop("basic model must be a normal confidence interval specified by list(mean, sd)");
	}
	m <- model$mean;
	s <- model$sd;
	alpha <- 1 - level;
	lower <- qnorm(alpha/2, m, s);
	upper <- qnorm(1 - alpha/2, m, s);

	mean(lower <= data & data <= upper)
}

# lag-one autocorrelation
autocorl1 <- function(x) {
	cor(x[-length(x)],x[-1])
}

#' Compute the Bayesian false discovery rate from posterior probability.
#'
#' @param v  sorted (by descending order) posterior probablity that a gene 
#'           is in the positive class
#' @return false discovery rate incurred by declaring each gene as significant 
#'         in order
bayesian_fdr <- function(v) {
	cumsum(1 - v) / 1:length(v)
}

#' Signal to noise ratio
#'
#' Estimate the signal to noise ratio.
#'
#' @param object  \code{gpldiff} object or \code{ldiff_data} object
#' @return \code{numeric} value
#' @export
#'
snr <- function(object) {
	UseMethod("snr", object)
}
