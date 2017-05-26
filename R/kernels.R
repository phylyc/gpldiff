#' Squared exponential kernel function
#'
#' This function evaluates the squared exponential kernel function with
#' parameters \code{nu^2} and \code{lambda^2}.
#'
#' k(xi, xj) = nu^2 exp( -1/(2 lambda^2) ||x_i − x_j||^2 )
#'
#' @param xi       first data point
#' @param xj       second data point
#' @param nu2      scale parameter (squared)
#' @param lambda2  characteristic length scale parameter (squared)
#' @export
squared_exponential_kernel <- function(xi, xj, nu2, lambda2) {
	d <- xi - xj;
	nu2 * exp(-0.5 * sum(d * d) / lambda2)
}

# Squared kernel function
# Part of the partial derivative of the squared exponential kernel function
# w.r.t. lambda^2
# k(xi, xj) = -1/(2 lambda^2) ||x_i − x_j||^2
squared_kernel <- function(xi, xj, lambda2) {
	d <- xi - xj;
	0.5 * sum(d * d) / lambda2
}

#' Apply the covariance function and construct the kernel matrix
#'
#' Kernel matrix is symmetric positive definite.
#'
#' @param x        \code{vector} or \code{matrix} (rows are observations)
#' @param covar_f  covariance function \code{k(xi, xj)}
#' @param ...      other parameters for the covariance function
#' 
#' @export
kernel_matrix <- function(x, covar_f, ...) {
	x <- as.matrix(x);
	n <- nrow(x);
	K <- matrix(0, nrow=n, ncol=n);
	K[lower.tri(K, diag=TRUE)] <-
		unlist(lapply(1:n,
			function(i) {
				unlist(lapply(i:n,
					function(j) {
						covar_f(x[i,], x[j,], ...)
					}
				))
			}
		))
	K[upper.tri(K)] <- t(K)[upper.tri(K)];
	K
}
