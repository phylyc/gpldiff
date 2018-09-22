#' Random generation of latent difference data
#'
#' Generate two non-linear data series with noise and missingness.
#' Response values y are observed with Gaussian noise, and never
#' simultaneously from both groups (missingness).
#' In other words, At each x value, only one of two groups is observed.
#' 
#' @param n       number of data points
#' @param xlim    domain limits
#' @param xsigma  standard deviation of noise to add to x values
#' @param pg      probabilities of each group
#' @param mu      global mean of data points
#' @param sigma   standard deviation of observation noise in y values
#' @return a \code{list} object with both observed and latent variables
#' @export
#'
rldiff <- function(n, xlim=c(0, 10), xsigma=0.1, pg=c(0.5, 0.5), mu=0.5, sigma=0.05) {
	x <- sort(seq(xlim[1], xlim[2], length.out=n) + rnorm(n, sd = xsigma));

	m_a <- rcurve(x)$y;
	m_a <- m_a - mean(m_a) + mu;
	m_b <- rcurve(x)$y;
	m_b <- m_b - mean(m_b) + mu;

	y_a <- rnorm(n, mean = m_a, sd = sigma);
	y_b <- rnorm(n, mean = m_b, sd = sigma);

	g <- sample(c(-0.5, 0.5), n, prob=pg, replace=TRUE);

	# either one group or the other is observed
	y <- ifelse(g > 0, y_b, y_a);

	# latent variable to be inferred
	f <- m_b - m_a;

	structure(
		list(
			# observed
			J = n,
			x = x,
			y = y,
			g = g,
			# latent
			f = f,
			m_a = m_a,
			m_b = m_b,
			y_a = y_a,
			y_b = y_b,
			mu = mu,
			sigma = sigma
		),
		class = "ldiff_data"
	)	
}

plot.ldiff_data <- function(object, fit) {
	par(mfrow=c(4, 1));
	
	with(object,
		{
			plot(c(x, x), c(m_a, m_b), xlab="x", ylab="y", type="n", main = "Truth data");
			lines(x, m_a, col=3);
			lines(x, m_b, col=4);

			plot(c(x, x), c(y_a, y_b), xlab="x", ylab="y", type="n", main = "data with noise");
			points(x, y_a, col=3);
			points(x, y_b, col=4);

			plot(x, y, col=as.numeric(g > 0)+3, main = "Observed data with noise and missingness");

			plot(x, f, main = "difference in mean", type="l", ylim=c(-1.5, 1.5));

			if (!is.null(fit)) {
				fsd <- sqrt(fit$predict$fvar);
				points(x, fit$params$f, main = "E[f]", col="red")
				lines(x, fit$params$f - 2 * fsd, main = "E[f]", col="red")
				lines(x, fit$params$f + 2 * fsd, main = "E[f]", col="red")
			}
		}
	)
}

#' Signal to noise ratio
#'
#' Estimate the signal to noise ratio.
#'
#' @param object \code{ldiff_data} object
#' @return \code{numeric} value
#' @export
#'
snr.ldiff_data <- function(object) {
	object$mu / object$sigma
}

