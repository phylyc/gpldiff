# NB K refers to the kernel matrix, unless specified otherwise
#    In the mathematical derivations, the \Sigma symbol is
#    used in place of K.

# Find mu^2 that maximizes the posterior of mu^2
# conditioned on all other parameters and hyperparameters
# Compute \hat{mu}
#   = \frac{\tau^2}{J \tau^2 + \sigma^2} \sum_{j=1}^{J} y_j - g_j f_j
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
#   S = K / sigma^2
#   D = G^2
#   G = I g
#   b = G (y - \mu \vec{1})
# For stability, compute (S^{-1} + D)^{-1} as S - S D^{1/2} B^{-1} D^{1/2} S
# where
#   B = I + D^{1/2} S D^{1/2}
# Therefore,
# (S^{-1} + D)^{-1} b = (S - S D^{1/2} B^{-1} D^{1/2} S) b
#                     = S (I - D^{1/2} B^{-1} D^{1/2} S) b
#                     = S (b - D^{1/2} B^{-1} D^{1/2} S b)
#                     = S a
# where
#   a = b - D^{1/2} B^{-1} D^{1/2} S b
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
	A <- solve(K %*% diag(g)) * sigma2 + diag(g);
	b <- y - mu;
	# \hat{f} = A \ b
	solve(A, b)
}

# Compute a = b - D^{1/2} B^{-1} D^{1/2} S b
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

	# NB Numerical difficulty can occur here!
	#    this may return a vector of NaNs
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
# - \frac{1}{2}\frac{\mu^2}{\tau^2} - \frac{1}{2} log | I + K W |
#
# where W = G^2 / sigma^2
#
# Kernel matrix K is pre-calculated based on hyperparameters \nu^2 and \lambda^2
#
# @param params   maximum a posteriori estimates of parameters
# @param hparams  fixed hyperparameter values
# @param hgradient  whether to calculate the gradient of the marginal
#                   likelihood with respect to the hyperparameters
#                   (Rasmussen & Williams 2006, p.125)
ll_g_hparams <- function(data, params, hparams, K=NULL, hgradient=FALSE) {

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
	W <- diag(g * g / sigma2);
	KW <- K %*% W;

	# Ideally,
	# ft.Kinv.f <- ft_kinv_f(y, g, f, mu, sigma2, K);
	# but vector a is needed later
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
		det(KW + eye, log=TRUE)
	);

	if (hgradient) {
		# compute explicit gradient of marginal likelihood w.r.t. to nu^2 and lambda^2

		nu2 <- hparams$nu2;
		lambda2 <- hparams$lambda2;

		E <- solve(eye + KW, KW);
		
		gll_nu2 <- 0.5 * (1/nu2) * (ft.Kinv.f - sum(diag(E)));

		V <- kernel_matrix(data$x, squared_kernel, lambda2=lambda2*lambda2);
		K.V <- K * V;
		gll_lambda2 <- 0.5 * 1/(sigma2*sigma2) * t(a) %*% K.V %*% a -
			0.5 * sum(diag( solve(eye + KW, K.V %*% W) ));

		list(evidence = ll, gradient = list(nu2 = gll_nu2, lambda2 = gll_lambda2))
	} else {
		ll
	}
}

# Compute f^\top K^{-1} f
# since
#   f = S a = (K / sigma2) a
# thus
#   f^\top K^{-1} f = ((K / sigma2) a)^\top K^{-1} f
#                   = 1/sigma2 a^\top K K^{-1} f
#                   = a^\top f / sigma2
ft_kinv_f <- function(y, g, f, mu, sigma2, K) {
	S <- K / sigma2;
	b <- g * (y - mu);
	a <- compute_a(S, b, g);
	sum(a * f) / sigma2
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
		alpha = 0.1,
		beta = 0.1,
		tau2 = 10
	);
}

# Fit model parameters with hyperparameters fixed
fit_params <- function(data, params, hparams, tol=1e-5, max.iter=10, predict=TRUE, plot=FALSE, fixed=NULL, hgradient=TRUE) {

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

	if (is.null(fixed)) {
		fixed <- list(
			mu = FALSE,
			sigma2 = FALSE,
			f = FALSE
		);
	}

	K <- kernel_matrix(data$x, squared_exponential_kernel, nu2=hparams$nu2, lambda2=hparams$lambda2);

	delta <- Inf;
	niters <- 0;
	while (delta > tol) {
		old <- unlist(params);

		if (!fixed$f) {
			params$f <- argmax_f_lp(data$y, data$g, params$mu, params$sigma2, K);
			if (any(is.nan(params$f))) {
				print(str(params$f));
				stop("Numerical difficulties encountered; `f` contains NaN values.")
			}
		}
		if (!fixed$mu) {
			params$mu <- argmax_mu_lp(data$y, data$g, params$f, params$sigma2, hparams$tau2);
			if (any(is.nan(params$mu))) {
				print(str(params$mu));
				stop("Numerical difficulties encountered; `mu` is NaN.")
			}
		}
		if (!fixed$sigma2) {
			params$sigma2 <- argmax_sigma2_lp(data$y, data$g, params$f, params$mu, hparams$alpha, hparams$beta);
			if (any(is.nan(params$sigma2))) {
				print(str(params$sigma2));
				stop("Numerical difficulties encountered: `sigma2` is NaN.")
			}
		}

		if (plot) {
			graphics::plot(data$x, params$f)
		}

		delta <- norm(as.matrix(old - unlist(params)), "F");
		if (is.na(delta)) {
			print(str(params));
			stop("Numerical difficulties encountered; `params` contain NaN values.")
		}

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

	res <- ll_g_hparams(data, params, hparams, K, hgradient=hgradient);
	model$evidence <- res$evidence;
	if (hgradient) {
		model$gradient <- res$gradient;
	}
	class(model) <- "gpldiff";

	model
}

# @return updated momentum
momentum <- function(v, g, learn.rate = 0.01, beta=0.9) {
	if (is.na(v)) {
		learn.rate * g
	} else {
		beta * v + (1 - beta) * learn.rate * g
	}
}

# Adaptive Moment Estimation
# Computes the first and second order moments of the momentum
# @param momentum   previous value of the moments
# @param g  gradient
# @param t  time step
# @return updated momentum
adam <- function(momentum, g, t, beta1=0.9, beta2=0.999) {
	c(
		(beta1 * momentum[1] + (1 - beta1) * g) / (1 - beta1^t),
		(beta2 * momentum[2] + (1 - beta2) * g*g) / (1 - beta2^t)
	)
}

# Calculate update value using ADAM learning rate
adam_step <- function(momentum, learn.rate = 0.01, eps = 1e-8) {
	learn.rate * momentum[1] / (sqrt(momentum[2]) + eps)
}

#' Fit Gaussian process latent difference model
#'
#' Fit a Gaussian process latent difference (GPLDIFF) model using maximum a posteriori
#' estimates of parameters. Hyperparameters are estimated by maximizing the
#' marginal likelihood conditional on the hyperparameters. The primary
#' parameter of interest is \code{f}, representing the latent group difference.
#' The worst-case time complexity is O(J^3) and space complexity is O(J^2).
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
#' @param data      observed data; a \code{list} of J, x, g, y;
#'                  or a \code{gpldiff_data} object
#' @param params    initial parameter values; a \code{list} of f, mu, sigma
#' @param hparams   hyperparameters; a \code{list} of nu2, lambda2, alpha, beta, tau2
#' @param adapt     whether to learn hyperparameters
#' @param learn.rate  learning rate for hyperparameters
#'                    (set to 0 to optimize hyperparameters using the slow
#'                    algorithm)
#' @param tol       minimum change in parameter values for convergence
#' @param tol2      minimum change in log likelihood for convergence
#' @param max.iter  maximum number of iterations at parameter level
#' @param max.iter2 maximum number of iterations at hyperparameter level
#' @param predict   whether to save variables required for prediction
#' @param verbose   verbosity level; none: 0, info: 1, debug: 2
#' @return \code{gpldiff} object
#' @export
#' @examples
#' \dontrun{
#' # `data` has been read in
#' hparams <- list(
#' 	nu2 = 1.13^2,
#' 	lambda2 = 1.41^2,
#' 	alpha = 0.1,
#' 	beta = 0.1,
#' 	tau2 = 1
#' )
#' params <- NULL
#' fit <- gpldiff(data, params, hparams)
#' plot(fit, data);
#' }
#'
gpldiff <- function(data, params=NULL, hparams=NULL, adapt=c("none", "GD", "GDM", "ADAM", "L-BFGS-B", "Brent"), learn.rate=0.2, tol=1e-1, tol2=1e-1, max.iter=10, max.iter2=10, predict=TRUE, verbose=1, ...) {
	if (is.null(hparams)) {
		hparams <- default_hparams();
	}

	if (verbose >= 1) {
		time.began <- proc.time();
		message("J: ", data$J);
	}

	adapt <- match.arg(adapt);

	check_data(data);

	niters <- 0;
	if (adapt != "none") {
		# find hyperparameter values lambda and nu that maximize the log marginal likelihood
		# NB other hyperparameters are fixed
		delta <- Inf;
		old.ll <- -Inf;

		if (adapt == "GDM") {
			m.lambda2 <- NA;
			m.nu2 <- NA;
		} else if (adapt == "ADAM") {
			m.lambda2 <- c(0.0, 0.0);
			m.nu2 <- c(0.0, 0.0);
		}

		while (delta > tol2) {
			niters <- niters + 1;
			if (verbose >= 2) {
				message("iteration ", niters)
			}

			if (adapt == "GD") {

				# Basic gradient descent
				# learning rate needs tuning!
				# Quits early if log likelihood decreases
				# TODO figure out which parameter is overshooting and decrease
				#      learning rate or step size for that parameter

				res <- fit_params(data, params, hparams, tol=tol, max.iter=max.iter, predict=FALSE,
													hgradient=TRUE);
				ll <- res$evidence;

				if (verbose >= 2) {
					message("hparams: ", hparams$nu2, " ", hparams$lambda2)
					message("hgrads: ", res$gradient$nu2, " ", res$gradient$lambda2)
				}

				if (ll < old.ll) {
					# ll has gotten worse: backtrack and reduce the learning rate
					hparams <- hparams.old;
					learn.rate <- learn.rate * 0.5;
					old.ll <- -Inf;
					if (verbose >= 2) {
						message("log evidence: ", ll)
						message("learn rate: ", learn.rate)
						message("backtracking ...")
					}

					next;
				} else {
					# continue updating hyperparameters
					hparams.old <- hparams;
					hparams$lambda2 <- hparams$lambda2 + learn.rate * res$gradient$lambda2;
					hparams$nu2 <- hparams$nu2 + learn.rate * res$gradient$nu2;
				}

			} else if (adapt == "GDM") {

				# Gradient descent with momentum update

				m.lambda2 <- momentum(m.lambda2, res$gradient$lambda2);
				m.nu2 <- momentum(m.nu2, res$gradient$nu2);

				hparams$lambda2 <- hparams$lambda2 + m.lambda2;
				hparams$nu2 <- hparams$nu2 + m.nu2;

			} else if (adapt == "ADAM") {

				# ADAM did not seem to work well: it keeps getting stuck due to
				# second moment becoming very large

				m.lambda2 <- adam(m.lambda2, -res$gradient$lambda2, niters);
				m.nu2 <- adam(m.nu2, -res$gradient$nu2, niters);
				step.lambda2 <- adam_step(m.lambda2, learn.rate=learn.rate);
				step.nu2 <- adam_step(m.nu2, learn.rate=learn.rate);
				hparams$lambda2 <- hparams$lambda2 - step.lambda2;
				hparams$nu2 <- hparams$nu2 - step.nu2;

				if (verbose >= 2) {
					message("ADAM m1: ", m.nu2[1], " ", m.lambda2[1]);
					message("ADAM m2: ", m.nu2[2], " ", m.lambda2[2]);
					message("ADAM step: ", step.nu2, " ", step.lambda2)
				}

			} else if (adapt == "L-BFGS-B") {

				# `optim` requires the objective function and gradient
				# function separately, this requires two expensive calls per
				# iteration within `optim`!

				# Due to lack of parameter constraint, CG and BFGS can often
				# run into numeric problems near -Inf and Inf,
				# causing the B matrix to be singular

				# With contraints added, optimize will force the
				# use of L-BFGS-B

				opt <- optim(
					c(
						lnu = log(sqrt(hparams$nu2)),
						llambda = log(sqrt(hparams$lambda2))
					),
					fn = function(par) {
						hparams$nu2 <- exp(par[1])^2;
						hparams$lambda2 <- exp(par[2])^2;
						if (verbose >= 2) {
							message("hparams: ", hparams$nu2, " ", hparams$lambda2)
						}
						fit_params(data, params, hparams, tol=tol, max.iter=max.iter,
											 predict=FALSE)$evidence
					},
					gr = function(par) {
						hparams$nu2 <- exp(par[1])^2;
						hparams$lambda2 <- exp(par[2])^2;
						unlist(fit_params(data, params, hparams, tol=tol, max.iter=max.iter,
											 predict=FALSE, hgradient=TRUE)$gradient)
					},
					lower = -10, upper = 10,
					method = "L-BFGS-B",
					control = list(abtol=tol, retol=tol, pgtol=tol, fnscale=-1, maxit=max.iter)
				);

				hparams$nu2 <- exp(opt$par[1])^2;
				hparams$lambda2 <- exp(opt$par[2])^2;
				ll <- opt$value;

			} else {

				opt.lambda <- optimize(
					function(lambda) {
						hparams$lambda2 <- lambda^2;
						fit_params(data, params, hparams, tol=tol, max.iter=max.iter,
											 predict=FALSE)$evidence
					},
					interval = c(0, max(2*IQR(data$x), 1)),
					maximum=TRUE,
					tol=tol
				);
				hparams$lambda2 <- opt.lambda$maximum^2;

				opt.nu <- optimize(
					function(nu) {
						hparams$nu2 <- nu^2;
						fit_params(data, params, hparams, tol=tol, max.iter=max.iter,
											 predict=FALSE)$evidence
					},
					interval = c(0, max(2*IQR(data$y), 1)),
					maximum=TRUE,
					tol=tol
				);
				hparams$nu2 <- opt.nu$maximum^2;

				ll <- opt.nu$objective;
			}

			if (verbose >= 2) {
				message("log evidence: ", ll);
			}

			# enforce bounds
			if (hparams$lambda2 <= 0) {
				hparams$lambda2 <- 1e-3;
			}
			if (hparams$nu2 <= 0) {
				hparams$nu2 <- 1e-3;
			}

			delta <- ll - old.ll;
			old.ll <- ll;

			if (niters >= max.iter2) break;
		}
	}

	model <- fit_params(data, params, hparams, tol=tol, max.iter=max.iter, predict=predict, ...);
	model$niters <- c(model$niters, niters);

	if (verbose >= 1) {
		message("final log evidence: ", model$evidence);
		elapsed <- proc.time() - time.began;
		message("elapsed wall time: ", elapsed[3]);
	}

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
		stop("`gpldiff` object must have been created by calling `gpldiff()` with `predict=TRUE`");
	}
	r <- z * sqrt(object$predict$fvar);
	cint <- with(object$param,
		data.frame(lower = f - r, upper = f + r)
	);
	cint
}

#' Empirical coverage probability of latent differences of a fitted GPLDIFF model
#'
#' @param model  \code{gpldiff} object
#' @param data   known latent differences
#' @param level  confidence level
#' @return empirical coverage probability
#'
coverage.gpldiff <- function(model, data, level=0.95) {
	if (is.null(model$predict)) {
		stop("`gpldiff` object must have been created by calling `gpldiff()` with `predict=TRUE`");
	}

	if (is.list(data)) {
		f <- data$f;
	} else {
		f <- data;
	}

	msd <- list(mean=model$params$f, sd=sqrt(model$predict$fvar));
	coverage(msd, f, level=level)
}

mean_center <- function(x) {
	x - mean(x)
}

#' Plot fitted GPLDIFF model
#'
#' @param model \code{gpldiff} object
#' @param data  data used to fit model
#' @param center  whether to center the response variable
#' @param estimated  whether to plot estimated points
#'
#' @method plot gpldiff
#' @export
#'
plot.gpldiff <- function(model, data, which=NULL, center=FALSE, estimated=FALSE, xlab="x") {
	if (!is.null(model$predict)) {
		cint <- confint(model);
		flim <- c(min(cint[,1]), max(cint[,2]));
		flim[1] <- flim[1] - diff(flim)*0.2;
	} else {
		cint <- NULL;
	}

	if (is.null(which)) {
		which <- c("response", "residual", "latent", "odds");
	}

	idx <- order(data$x);

	if (center) {
		data$y[data$g <= 0] <- mean_center(data$y[data$g <= 0]);
		data$y[data$g > 0] <- mean_center(data$y[data$g > 0]);
	}
		
	par(mfrow=c(length(which),1), mai=c(0.8, 0.9, 0.1, 0.5));

	g <- data$g[idx] > 0;
	yhat <- (model$params$mu + model$params$f * data$g);

	if ("response" %in% which) {
		# plot observed data
		ylim <- range(data$y);
		ylim[1] <- ylim[1] - diff(ylim)*0.2;
		plot(NA, xlim=range(data$x), ylim=ylim, xlab="", ylab="observed responses", las=1);
		lines(data$x[idx][!g], data$y[idx][!g], col="#0073C2FF", pch=20, type="b", lwd=2);
		lines(data$x[idx][g], data$y[idx][g], col="#EFC000FF", pch=20, type="b", lwd=2);
		legend("bottomright", inset=0.01, col=c("#EFC000FF", "#0073C2FF"), lwd=2, legend=c("case", "control"), bty="n");
		# plot estimated data
		if (estimated) {
			points(data$x[idx][!g], yhat[idx][!g], col="#0073C2FF", pch=21, type="b")
			points(data$x[idx][g], yhat[idx][g], col="#EFC000FF", pch=21, type="b")
			legend("bottomleft", inset=0.01, pch=c(20, 21), legend=c("observed", "estimated"), bty="n");
		}
	}

	if ("residual" %in% which) {
		# plot residual of y_hat
		r <- data$y - yhat;
		rlim <- range(r);
		rlim[1] <- rlim[1] - diff(rlim)*0.2;
		cols <- c("#0073C2FF", "#EFC000FF")[as.integer(g)+1];
		plot(data$x[idx], r[idx], col=cols, pch=20, xlab="", ylab="residual", ylim=rlim, las=1);
		abline(h = 0, col="grey30", lty=2);
		legend("bottomleft", inset=0.01, pch=20, col=c("#EFC000FF", "#0073C2FF"), legend=c("case", "control"), bty="n");
	}

	# plot latent difference f_hat
	if ("latent" %in% which) {
		plot(NA, xlim=range(data$x), ylim=flim, xlab="", ylab="latent difference f", las=1);
		if (!is.null(data$f)) {
			points(data$x, data$f, pch=20, col="#868686FF");
			legend("bottomleft", inset=0.01, col=c("#868686FF", "#CD534CFF"), pch=c(20, 21), legend=c("truth", "estimated"), bty="n");
		}
		abline(h = 0, col="grey30", lty=2);
		points(data$x[idx], model$params$f[idx], col="#CD534CFF", pch=21);
		if (!is.null(cint)) {
			lines(data$x[idx], cint[idx,1], col="#CD534C55", lwd=2);
			lines(data$x[idx], cint[idx,2], col="#CD534C55", lwd=2);
		}
	}

	# plot log posterior odds
	if ("odds" %in% which) {
		lodds.cut <- 5;
		lodds <- summary(model, log.odds=TRUE);
		plot(NA, xlim=range(data$x[idx]), ylim=range(lodds[idx]),
			xlab=xlab, ylab="log posterior odds", las=1);
		abline(h = 0, col="grey30", lty=2);
		abline(h = -lodds.cut, col="grey30");
		abline(h = lodds.cut, col="grey30");
		lines(data$x[idx], lodds[idx], col="#CD534CFF", pch=20, lwd=2);
	}
}

#' Summarize GPLDIFF model
#'
#' Calculate the posterior probability that \code{f > 0}.
#' 
#' @param object \code{gpldiff} object
#' @method summary gpldiff
#' @export
#'
summary.gpldiff <- function(object, log.odds = FALSE, ...) {
	if (is.null(object$predict)) {
		stop("`gpldiff` object must have been created by calling `gpldiff()` with `predict=TRUE`");
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

#' Signal to noise ratio
#'
#' Estimate the signal to noise ratio.
#'
#' @param object \code{gpldiff} object
#' @return \code{numeric} value
#'
#' @method snr gpldiff
#' @export
#'
snr.gpldiff <- function(object) {
	object$params$mu / sqrt(object$params$sigma2)
}

check_data <- function(data) {
	stopifnot(data$J > 0);
	stopifnot(!any(is.na(data$x)));
	stopifnot(!any(is.na(data$y)));
	stopifnot(!any(is.na(data$g)));
}

