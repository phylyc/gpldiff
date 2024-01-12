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
#' @param sigma2  variance of observation noise in y values
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
			sigma2 = sigma*sigma
		),
		class = "ldiff_data"
	)	
}

#' @export
plot.ldiff_data <- function(object, fit) {
	par(mfrow=c(4, 1));
	
	with(object,
		{
			ylim <- range(c(y_a, y_b));
			ylim[1] <- ylim[1] - diff(ylim)*0.5;

			plot(c(x, x), c(m_a, m_b), ylim=ylim, xlab="", ylab="y", type="n", main = "Truth data", las=1);
			lines(x, m_a, col="#0073C2FF", lwd=2);
			lines(x, m_b, col="#EFC000FF", lwd=2);
			legend("bottomright", inset=0.01, col=c("#EFC000FF", "#0073C2FF"), legend=c("case", "control"), lwd=2, bty="n");

			plot(c(x, x), c(y_a, y_b), ylim=ylim, xlab="", ylab="y", type="n", main = "Data observed with noise", las=1);
			points(x, y_a, pch=19, col="#0073C2FF");
			points(x, y_b, pch=19, col="#EFC000FF");
			segments(x, y_a, x, y_b, col="#86868633", lwd=2);
			legend("bottomright", inset=0.01, col=c("#EFC000FF", "#0073C2FF"), legend=c("case", "control"), pch=19, bty="n");

			cols <- c("#0073C2FF", "#EFC000FF")[as.numeric(g >= 0)+1];
			plot(x, y, ylim=ylim, xlab="", col=cols, pch=20, main = "Data observed with noise and missingness", las=1);
			points(x, y_a, pch=21, col="#0073C2FF");
			points(x, y_b, pch=21, col="#EFC000FF");
			segments(x, y_a, x, y_b, col="#86868633", lwd=2);
			legend("bottomright", inset=0.01, col=c("#EFC000FF", "#0073C2FF"), legend=c("case", "control"), pch=19, bty="n");
			legend("bottomleft", inset=0.01, legend=c("observed", "missing"), pch=c(19, 21), bty="n");

			flim <- range(f)*1.8;
			flim[1] <- flim[1] - diff(flim)*0.4;
			plot(x, f, ylim=flim, xlab="x", main = "Latent difference", type="p", col="#868686FF", pch=19, las=1);

			if (!is.null(fit)) {
				fsd <- sqrt(fit$predict$fvar);
				points(x, fit$params$f, main = "E[f]", col="#CD534CFF")
				lines(x, fit$params$f - 2 * fsd, main = "E[f]", col="#CD534C55", lwd=2)
				lines(x, fit$params$f + 2 * fsd, main = "E[f]", col="#CD534C55", lwd=2)
				legend("bottomleft", inset=0.01, col=c("#868686FF", "#CD534CFF"), legend=c("truth", "estimated"), pch=c(19, 21), bty="n");
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
	object$mu / sqrt(object$sigma2)
}

