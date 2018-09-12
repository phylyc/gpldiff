# NB K refers to the kernel matrix, unless specified otherwise
#    In the mathematical derivations, the \Sigma symbol is
#    used in place of K.

# Find mu^2 that maximizes the posterior of mu^2
# conditioned on all other parameters and hyperparameters
# Compute \hat{mu}
#   = \frac{\tau^2}{\tau^2 + \sigma^2} \sum_{j=1}^{J} y_j - g_j f_j
argmax_mu_lp <- function(y, g, f, sigma2, tau2) {
	tau2 / (length(y) * tau2 + sigma2) * sum(y - g * f)
}

# Find sigma^2 that maximizes the posterior of sigma^2
# conditioned on all other parameters and hyperparameters
# Compute \hat{sigma^2} = (\beta + (1/2) ||y - \mu - G f||^2) / (J/2 + \alpha + 1)
argmax_sigma2_lp <- function(y, g, f, mu, alpha, beta) {
	d <- y - mu - (g * f);
	( beta + 0.5 * sum(d*d) ) / ( length(y) / 2 + alpha + 1 )
}


# Find f that maximizes the posterior of f
# conditioned on all other parameters and hyperparameters
# Compute \hat{f} = (S^{-1} + D)^{-1} b
# by the Rasmussen-Williams method,
# where
#   S = K / sigma2
#   D = G^2
#   G = diag(g)
#   b = G (y - \mu \vec{1})
# For stability, compute (S^{-1} + D)^{−1} as S − S D^{1/2} B^{-1} D^{1/2} S
# where
#   B = I + D^{1/2} S D^{1/2}
# Therefore,
# (S^{-1} + D)^{-1} b = (S − S D^{1/2} B^{-1} D^{1/2} S) b
#                     = S (I − D^{1/2} B^{-1} D^{1/2} S) b
#                     = S (b − D^{1/2} B^{-1} D^{1/2} S b)
#                     = S a
# where
#   a = b − D^{1/2} B^{-1} D^{1/2} S b
# Rasmussen & Williams 2006
argmax_f_lp_rw <- function(y, g, mu, sigma2, K) {
	S <- K / sigma2;
	b <- g * (y - mu);
	a <- compute_a(S, b, g);
	S %*% a
}

# Find f that maximizes the posterior of f
# conditioned on all other parameters and hyperparameters
# Compute \hat{f} = (S^{-1} + D)^{-1} b directly
# NB Solve(K) costs O(n^3) time
# NB Without regularization of K, solve(K) may not be invertible
argmax_f_lp_direct <- function(y, g, mu, sigma2, K) {
	A <- solve(K) * sigma2 + diag(g*g);
	b <- g * (y - mu);
	# \hat{f} = A \ b
	solve(A, b)
}

# Compute a = b − D^{1/2} B^{-1} D^{1/2} S b
# where
#   B = I + D^{1/2} S D^{1/2}
# Rasmussen & Williams 2006
compute_a <- function(S, b, dsqrt) {
	n <- nrow(S);
	Dsqrt <- diag(dsqrt);

	B <- diag(1, nrow=n, ncol=n) + (Dsqrt %*% S %*% Dsqrt);

	# decompose B into L U
	U <- chol(B);
	L <- t(U);

	# thus, B^{-1} = (L U)^{-1} = U^{-1} L^{-1}
	# A^{-1} b == solve(A, b)

	b - Dsqrt %*% solve(U, solve(L, Dsqrt %*% (S %*% b)))
}

# Compute \hat{f} = (S^{-1} + D)^{-1} b
# by first inverting A = S^{-1} + D
# where
#   S = K / sigma2
#   D = G^2
#   G = diag(g)
#   b = G (y - \mu)
argmax_f_lp_invert <- function(y, g, mu, sigma2, K) {
	Ainv <- invert_add_invert_rw(K / sigma2, g);
	b <- g * (y - mu);
	Ainv %*% b
}

argmax_f_lp <- argmax_f_lp_rw;

# Compute the covariance matrix of the Laplace approximation
# of the posterior of f conditioned on the hyperparameters
# Compute \hat{V} = ( K^{-1} + G^2 / sigma2 )^{-1}
# NB define D = G^2 / sigma2
# NB invert_add_invert_rw takes Dsqrt
f_laplace_covariance <- function(g, sigma2, K) {
	invert_add_invert_rw(K, g / sqrt(sigma2))
}

# Compute variance vector of the Lapalce approximation
# of the posterior of f conditioned on the hyperparameters
# This vector consists of the diagonal elements of the covariance matrix
f_laplace_variance <- function(g, sigma2, K) {
	# TODO compute only diagonal elements
	diag(f_laplace_covariance(g, sigma2, K))
}

# Compute the laplace approximation of the log marginal likelihood (model evidence)
#
# Ideally, all parameters would be marginalized out
#   p(y, \phi | x) = \int p(y | x, \theta, \phi) d\theta
# where \theta encompasses all the parameters (f, \mu, \sigma^2)
#       \phi are the hyperparameters (\tau^2, \alpha, \beta, \nu^2, \lambda^2)
#       (\phi is often omitted in standard notation)
#       x is the fixed independent variable
#
# However, in order to derive a closed form approximation,
# we treat parameters \mu and \sigma^2 as fixed; thus, their hyperparameters
# \tau^2, \alpha, and \beta must also be treated as fixed.
#
# Then, following the approach in Rasmussen & Williams 2006 (RW06), we use second order
# Taylor expansion of the joint density \psi(f) at \hat{f} to approximate \psi(f),
# resulting an approximation of the log marginal likelihood similar to
# equation 3.32 (RW06):
#
# log q(y | x) \approx -(J + 2\alpha + 2) log \sigma - \frac{\beta}{\sigma^2}
# - \frac{1}{2} \frac{1}{\sigma^2} (y - \mu \vec{1} - G \hat{f})^\top(y - \mu \vec{1} - G \hat{f})
# - \frac{1}{2} \hat{f}^\top K^{-1} \hat{f} - log \tau
# - \frac{1}{2}\frac{\mu^2}{\tau^2} - \frac{1}{2} log | I + K D |
#
# Kernel matrix K is pre-calculated based on hyperparameters \nu^2 and \lambda^2
#
# @param params   maximum a posteriori estimates of parameters
# @param hparams  fixed hyperparameter values
# TODO  Only the explicit gradient is calculated
#       To be useful, the implicit gradient is needed too!
#       (Rasmussen & Williams 2006, p.125)
ll_g_hparams <- function(data, params, hparams, K=NULL, gradient=FALSE) {

	if (is.null(K)) {
		K <- kernel_matrix(data$x, squared_exponential_kernel, nu2=hparams$nu2, lambda2=hparams$lambda2);
	}

	mu <- params$mu;
	sigma2 <- params$sigma2;
	f <- params$f;

	J <- length(data$y);
	y <- data$y;
	g <- data$g;

	alpha = hparams$alpha;
	beta = hparams$beta;
	tau2 = hparams$tau2;

	d <- y - mu - (g * f);
	# NB  this D is different from D in other places, hence D'
	Dp <- diag(g * g / sigma2);
	KDp <- K %*% Dp;

	# compute f^\top K^{-1} f
	# since
	#   f = S a = (K / sigma2) a
	# thus
	#   f^\top K^{-1} f = ((K / sigma2) a)^\top K^{-1} f
	#                   = 1/sigma2 a^\top K K^{-1} f
	#                   = a^\top f / sigma2
	S <- K / sigma2;
	b <- g * (y - mu);
	a <- compute_a(S, b, g);
	ft.Kinv.f <- sum(a * f) / sigma2;

	eye <- diag(1, nrow=J, ncol=J);

	# laplace approximation of the log marginal likelihood
	ll <- - 0.5 * (
		(J + 2*alpha + 2) * log(sigma2) +
		2 * beta / sigma2 +
		sum(d * d) / sigma2 +
		log(tau2) +
		mu * mu / tau2 +
		ft.Kinv.f +
		det(KDp + eye, log=TRUE)
	);

	if (gradient) {
		# TODO verify
		# compute explicit gradient of marginal likelihood w.r.t. to nu^2 and lambda^2

		nu2 <- hparams$nu2;
		lambda2 <- hparams$lambda2;

		#E <- solve(eye + solve(KDp));
		E <- solve(eye + KDp, KDp);
		
		gll_nu2 <- 0.5 * (1/nu2) * (ft.Kinv.f - sum(diag(E)));

		V <- kernel_matrix(data$x, squared_kernel, lambda2=lambda2*lambda2);
		gll_lambda2 <- 0.5 * ((t(a) %*% V) %*% f - sum((diag(E) * diag(V))));

		c(ll = ll, gll_nu2 = gll_nu2, gll_lambda2 = gll_lambda2)
	} else {
		ll
	}
}

# Compute (S^{-1} + D)^{-1}
# where S is symmetric positive definite and D is diagonal
invert_add_invert_direct <- function(S, D) {
	# solve directly without approximation
	solve(solve(S) + D)
}

# Compute (S^{-1} + D)^{−1} as S − S D^{1/2} B^{-1} D^{1/2} S
# where S is symmetric positive definite and D is diagonal
# redefined expression is more numerically stable (and does not require K to
# be regularized), because B is guaranteed to be well-conditioned
# where B = I + D^{1/2} S D^{1/2}
# Rasmussen & Williams 2006
invert_add_invert_rw <- function(S, dsqrt) {
	n <- nrow(S);
	Dsqrt <- diag(dsqrt);
	B <- diag(1, nrow=n, ncol=n) + (Dsqrt %*% S %*% Dsqrt);
	# compute S - S %*% Dsqrt %*% solve(B) %*% Dsqrt %*% S
	U <- chol(B);
	L <- t(U);
	S - S %*% Dsqrt %*% solve(U, solve(L, Dsqrt %*% S))
}

# Default hyperparameter values
default_hparams <- function() {
	list(
		nu2 = 1,
		lambda2 = 1,
		alpha = 2,
		beta = 1,
		tau2 = 1
	);
}

# Fit model parameters with hyperparameters fixed
fit_params <- function(data, params, hparams, tol=1e-5, max.iter=10, predict=TRUE, plot=FALSE) {

	if (is.null(hparams)) {
		hparams <- default_hparams();
	}

	# initial guess
	if (is.null(params)) {
		params <- list(
			mu = 0,
			sigma2 = 1,
			f = rep(0, data$J)
		);
	}

	K <- kernel_matrix(data$x, squared_exponential_kernel, nu2=hparams$nu2, lambda2=hparams$lambda2);

	delta <- Inf;
	niters <- 0;
	while (delta > tol) {
		old <- unlist(params);
		params$f <- argmax_f_lp(data$y, data$g, params$mu, params$sigma2, K);
		params$mu <- argmax_mu_lp(data$y, data$g, params$f, params$sigma2, hparams$tau2);
		params$sigma2 <- argmax_sigma2_lp(data$y, data$g, params$f, params$mu, hparams$alpha, hparams$beta);

		if (plot) {
			graphics::plot(data$x, params$f)
		}

		delta <- norm(as.matrix(old - unlist(params)), "F");
		niters <- niters + 1;
		if (niters >= max.iter) break;
	}

	model <- list(params = params, hparams = hparams, niters=niters);
	if (predict) {
		model$predict <- list(
			fvar = f_laplace_variance(data$g, params$sigma2, K),
			K = K
		);
	}
	model$evidence <- ll_g_hparams(data, params, hparams, K);
	class(model) <- "gpldiff";

	model
}

#' Fit Gaussian process latent difference model
#'
#' Fit a Gaussian process latent difference (GPLDIFF) model using maximum a posteriori
#' estimates of parameters. Hyperparameters are estimated by maximizing the
#' marginal likelihood conditional on the hyperparameters. The primary
#' parameter of interest is \code{f}, representing the latent group difference.
#'
#' The model is
#'
#' \code{
#'    K_{ij} = k(x_i, x_j);
#'    f ~ multi_normal(0, K);
#'	  mu ~ normal(0, tau);
#'	  sigma^2 ~ inv_gamma(alpha, beta);
#'	  y ~ normal(mu + (f .* g), sigma);
#' }
#'
#' where \code{k(x_i, x_j)} is the squared exponential covariance function and
#' \code{.*} denotes elementwise multiplication.
#' 
#' Maximum a posteriori estimates of parameters are found by coordinate ascent 
#' using analytical gradients. The posterior of \code{f} (conditioned on the
#' hyperparameters) is approximated by Laplace approximation. The posterior of
#' parameters \code{mu} and \code{sigma} are approximated by delta functions.
#' Hyperparameters are optimized by numerical optimization.
#'
#' @param data      observed data; list of x, g, y
#' @param params    initial parameter values; list of f, mu, sigma
#' @param hparams   hyperparameters; list of nu2, lambda2, alpha, beta, tau2
#' @param adapt     whether to learn hyperparameters
#' @param tol       tolerance at parameter level optimization
#' @param tol2      tolerance at hyperparameter level optimization
#' @param max.iter  maximum number of iterations at parameter level
#' @param max.iter2 maximum number of iterations at hyperparameter level
#' @param predict   whether to save variables required for prediction
#' @return \code{gpldiff} object
#' @export
#' @examples
#' \dontrun{
#' # `data` has been read in
#' hparams <- list(
#' 	nu2 = 1.13^2,
#' 	lambda2 = 1.41^2,
#' 	alpha = 2,
#' 	beta = 1,
#' 	tau2 = 1
#' )
#' params <- NULL
#' fit <- gpldiff(data, params, hparams)
#' plot(fit, data);
#' }
#'
gpldiff <- function(data, params=NULL, hparams=NULL, adapt=FALSE, tol=1e-1, tol2=1e-1, max.iter=10, max.iter2=10, predict=TRUE, ...) {
	if (is.null(hparams)) {
		hparams <- default_hparams();
	}

	niters <- 0;
	if (adapt) {
		# find hyperparameter values lambada and nu that maximize the log marginal likelihood
		# NB other hyperparameters are fixed
		delta <- Inf;
		while (delta > tol) {
			old <- unlist(hparams);

			opt.lambda <- optimize(
				function(lambda) {
					hparams$lambda2 <- lambda^2;
					fit_params(data, params, hparams, tol=tol, max.iter=max.iter, predict=FALSE)$evidence
				},
				interval = c(0, max(2*IQR(data$x), 1)),
				maximum=TRUE,
				tol=tol2
			);
			hparams$lambda2 <- opt.lambda$maximum^2;

			opt.nu <- optimize(
				function(nu) {
					hparams$nu2 <- nu^2;
					fit_params(data, params, hparams, tol=tol, max.iter=max.iter, predict=FALSE)$evidence
				},
				interval = c(0, max(2*IQR(data$y), 1)),
				maximum=TRUE,
				tol=tol2
			);
			hparams$nu2 <- opt.nu$maximum^2;

			delta <- norm(as.matrix(old - unlist(hparams)), "F");

			niters <- niters + 1;
			if (niters >= max.iter2) break;
		}
	}

	model <- fit_params(data, params, hparams, tol=tol, max.iter=max.iter, predict=predict, ...);
	model$niters <- c(model$niters, niters);

	model
}

#' Calculate confidence interval for fitted GPLDIFF model
#'
#' @param object \code{gpldiff} object
#' @param parm   not used (only parameter \code{f} is assessed)
#' @param level  confidence level
#'
#' @export
#' @examples
#' \dontrun{
#' cint <- confint(fit);
#' with(data, mean(f >= cint$lower & f <= cint$upper));
#' }
#'
confint.gpldiff <- function(object, parm=NULL, level=0.95, ...) {
	alpha <- 1 - level;
	z <- qnorm(1 - alpha/2);
	if (is.null(object$predict)) {
		stop("gpldiff object must be created with `predict=TRUE`");
	}
	r <- z * sqrt(object$predict$fvar);
	cint <- with(object$param,
		data.frame(lower = f - r, upper = f + r)
	);
	cint
}

mean_center <- function(x) {
	x - mean(x)
}

#' Plot fitted GPLDIFF model
#'
#' @param model \code{gpldiff} object
#' @param data  data used to fit model
#' @export
#'
plot.gpldiff <- function(model, data, center=FALSE) {
	if (!is.null(model$predict)) {
		cint <- confint(model);
		ylim <- c(min(cint[,1]), max(cint[,2]));
	} else {
		cint <- NULL;
	}

	idx <- order(data$x);

	if (center) {
		data$y[data$g < 0] <- mean_center(data$y[data$g < 0]);
		data$y[data$g > 0] <- mean_center(data$y[data$g > 0]);
	}
		
	par(mfrow=c(3,1), mai=c(0.6, 0.7, 0.1, 0.5));	

	# plot observed data
	g <- data$g[idx] <= 0;
	plot(NA, xlim=range(data$x), ylim=range(data$y), xlab="x", ylab="observed responses");
	lines(data$x[idx][g], data$y[idx][g], col="grey", pch=20, type="b", lwd=2);
	lines(data$x[idx][!g], data$y[idx][!g], col="orange", pch=20, type="b", lwd=2);

	# plot latent difference f
	plot(NA, xlim=range(data$x), ylim=ylim, xlab="x", ylab="latent difference f", las=1);
	if (!is.null(data$f)) {
		points(data$x, data$f);
	}
	abline(h = 0, col="grey30", lty=2);
	lines(data$x[idx], model$params$f[idx], col="grey30", lwd=2, pch=20, type="b");
	if (!is.null(cint)) {
		lines(data$x[idx], cint[idx,1], col="grey");
		lines(data$x[idx], cint[idx,2], col="grey");
	}

	# plot log posterior odds
	lodds <- summary(model, log.odds=TRUE);
	plot(data$x[idx], lodds[idx],
		xlab="x", ylab="log posterior odds", col="red", type="b", pch=20, lwd=2, ylim=c(0, max(lodds[is.finite(lodds)])*1.1));
	abline(h = 0, col="grey30", lty=2);
}

#' Summarize GPLDIFF model
#'
#' Calculate the posterior probability that \code{f > 0}.
#' 
#' @param object \code{gpldiff} object
#' @export
#'
summary.gpldiff <- function(object, log.odds = FALSE, ...) {
	if (is.null(object$predict)) {
		stop("gpldiff object must be created with `predict=TRUE`");
	}

	f <- object$params$f;
	fsd <- sqrt(object$predict$fvar);

	if (log.odds) {
		pnorm(0, mean=f, sd=fsd, lower.tail=FALSE, log=TRUE) -
			pnorm(0, mean=f, sd=fsd, lower.tail=TRUE, log=TRUE)
	} else {
		# Pr(f > 0) = 1 - Pr(f <= 0)
		pnorm(0, mean=f, sd=fsd, lower.tail=FALSE)
	}
}

subset_gpldiff <- function(x, start, end) {
	idx <- which(x$data$x >= start & x$data$x <= end);
	d <- x$data;
	m <- x$model;

	x.sub <- x;
	x.sub$data <- list(
		J = length(idx),
		x = d$x[idx],
		g = d$g[idx],
		y = d$y[idx]
	);
	x.sub$model$params$f <- m$params$f[idx];
	x.sub$model$predict$fvar <-m$predict$fvar[idx];
	x.sub$model$predict$K <-m$predict$K[idx, idx];

	x.sub
}

