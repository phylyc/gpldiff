#' Random regeneration of a spline curve
#'
#' Generate a set of random anchor data points and fit a spline curve through them.
#'
#' To generate curves with more variation, increase \code{m}, decrease \code{alpha}
#' and \code{beta} parameters towards 0.
#'
#' @param x  values of the independent variable
#' @param m  number of anchor data points
#' @param alpha  positive shape parameter 1 of beta distribution
#' @param beta   positive shape parameter 2 of beta distribution
#' @param closed  whether to anchor start and end points at 0
#' @param plot   whether to plot the curve
#' @return list of curve values \code{y} and curve function \code{f}.
#' @export
#'
rcurve <- function(x, m=10, alpha=0.5, beta=0.5, closed=FALSE, plot=FALSE) {
	N <- length(x);
	steps <- rpois((m-2), N/(m-1));
	idx <- c(1, cumsum(steps), N);
	valid <- idx <= N;
	idx <- idx[valid];
	m <- sum(valid);
	if (closed) {
		z <- c(0, rbeta(m-2, alpha, beta), 0);
	} else {
		z <- rbeta(m, alpha, beta);
	}
	sf <- splinefun(x[idx], z, method="natural");

	if (plot) {
		base::plot(x, sf(x), type="l", col="red")
		points(x[idx], z)
	}

	list(
		y = sf(x),
		f = sf
	)
}

#' Random regeneration of a polynomial curve
#'
#' Generate a random polynomial curve by drawing coefficients from the beta
#' distribution, scaled to range (-1, 1).
#'
#' @param d  degree of the polynomial
#' @param alpha  positive shape parameter 1 of beta distribution
#' @param beta   positive shape parameter 2 of beta distribution
#' @param plot   whether to plot the curve
#' @return list of curve values \code{y}, curve function \code{f}, and
#'         polynomial coefficients \code{beta}
#' @export
#'
rpolynomial <- function(x, d=10, alpha=1, beta=1, plot=FALSE) {
	# coefficients
	beta <- rbeta(d, alpha, beta) * 2 - 1;;

	f <- function(x) {
		X <- matrix(unlist(lapply(0:(d-1), function(n) x^n)), nrow=length(x));
		X %*% beta
	}

	if (plot) {
		plot(x, f(x), type="l", col="red")
	}

	list(
		y = f(x),
		f = f,
		beta = beta
	)
}

