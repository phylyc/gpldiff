library(devtools);
load_all("../R");

rcurve <- function(x, nknots=12) {
	z <- rep(0, N);
	nonzeros <- round(runif(rpois(1, 10), 1, N));
	z[nonzeros] <- runif(length(nonzeros));
	w <- rbeta(N, 1, 10);
	w[nonzeros] <- 1;
	s <- smooth.spline(x, z, w, lambda=1e-4, nknots=nknots, keep.data=FALSE);
	s$y
}

rpolynomial <- function(x, dmin=c(1, 10)) {
	# degree of polynomial
	d <- round(runif(1, dmin[1], dmin[2]));
	# coefficients
	beta <- rbeta(d, 0.5, 0.5) - 0.5;
	X <- matrix(unlist(lapply(0:(d-1), function(n) x^n)), nrow=length(x));

	list(
		y = X %*% beta,
		beta = beta
	)
}

# simulate non-linear data of two groups with missingness

set.seed(6);

# generate data

N <- 200;
offset <- 0;
sigma <- 0.05;

x <- sort(seq(0, 10, length.out=N) + rnorm(N, sd = 0.1));

m_a <- rcurve(x);
m_b <- rcurve(x) + offset;

#m_a <- sin(x);
#m_b <- sin(x/(pi/2) - 7*pi/2) + offset;
mu <- mean((m_b + m_a) / 2);

y_a <- rnorm(N, mean = m_a, sd = sigma);
y_b <- rnorm(N, mean = m_b, sd = sigma);

g <- sample(c(-0.5, 0.5), N, replace=TRUE);
#g <- sample(c(-1, 1), N, replace=TRUE);
# either one group or the other is observed
# group probabilities are equal; each group has ~50% missingness
y <- ifelse(g > 0, y_b, y_a);

# latent variable to be inferred
f <- m_b - m_a;

data <- list(
	# inputs
	J = N,
	x = x,
	y = y,
	g = g,
	# answers
	f = f,
	y_a = y_a,
	y_b = y_b,
	sigma = sigma
);

hparams <- list(
	nu2 = 1,
	lambda2 = 2,
	alpha = 0.1,
	beta = 0.1,
	tau2 = 10
);

# initial guess
params <- list(
	mu = 0,
	#mu = mu,
	sigma2 = 2,
	#sigma2 = sigma^2,
	f = rep(0, data$J)
	#f = f
);

fixed <- list(
	mu = FALSE,
	sigma2 = FALSE,
	f = FALSE
);

# sigma2 needs to be overestimated for better fit... why?
# set to 0.1^2 = 0.01 but estimate is ~0.3
# if sigma2 is constrained to known value, then 
# the fit of f is much poorer

fit <- gpldiff(data, params=params, hparams=hparams, fixed=fixed);
str(fit)

# Optimizing hyperparameters now works!
load_all("../R");
fit.adapt <- gpldiff(data, adapt="GD", tol2=1e-2, max.iter2=50, learn.rate=0.2, verbose=TRUE);
fit.adapt2 <- gpldiff(data, adapt="Brent", tol2=1e-2, max.iter2=20, learn.rate=0, verbose=TRUE);
fit.adapt3 <- gpldiff(data, adapt="L-BFGS-B", tol2=1e-2, max.iter2=20, learn.rate=0, verbose=TRUE);

str(fit.adapt)
str(fit.adapt2)

fit$evidence
fit.adapt$evidence
fit.adapt2$evidence
fit.adapt3$evidence

fit <- fit.adapt;
fit <- fit.adapt2;

# plot data

par(mfrow=c(4, 1));

plot(c(x, x), c(m_a, m_b), xlab="x", ylab="y", type="n", main = "Truth data");
lines(x, m_a, col=3);
lines(x, m_b, col=4);

plot(c(x, x), c(y_a, y_b), xlab="x", ylab="y", type="n", main = "data with noise");
points(x, y_a, col=3);
points(x, y_b, col=4);

plot(x, y, col=as.numeric(g > 0)+3, main = "Observed data with noise and missingness");

plot(x, f, main = "difference in mean", type="l", ylim=c(-0.6, 0.6));
fsd <- sqrt(fit$predict$fvar);
points(x, fit$params$f, main = "E[f]", col="red")
lines(x, fit$params$f - 2 * fsd, main = "E[f]", col="red")
lines(x, fit$params$f + 2 * fsd, main = "E[f]", col="red")


library(ggplot2);
library(io);

gcols <- c(control="royalblue", case="orangered");

truth <- data.frame(
	x = c(x, x), y = c(m_a, m_b),
	group = factor(
		c(rep("control", length(x)), rep("case", length(x))),
		levels=names(gcols)
	)
);

qdraw(
	ggplot(truth, aes(x, y, colour=group)) + theme_bw() +
		scale_colour_manual(values=gcols) +
		geom_line(lwd=1) +
		xlab("position") + ylab("response") +
		ylim(-2, 2)
	,
	height = 2, width = 6,
	file = "truth-data.pdf"
);

qdraw(
	ggplot(data.frame(x=x, y=y, group = factor(g, labels=names(gcols))), aes(x, y, colour=group)) + theme_bw() +
		scale_colour_manual(values=gcols) +
		geom_point(alpha=0.5) +
		geom_spoke(aes(y=-1.8, colour=group), angle=pi/2, radius=0.2) +
		xlab("position") + ylab("observed response") +
		ylim(-2, 2)
	,
	height = 2, width = 6,
	file = "observed-data.pdf"
);

z <- qnorm(1 - 0.05/2);
d <- data.frame(
	x = x,
	f = f,
	f_hat = fit$params$f,
	f_hat_min = fit$params$f - z * sqrt(fit$predict$fvar),
	f_hat_max = fit$params$f + z * sqrt(fit$predict$fvar)
);

cols <- c(truth="grey20", estimate="firebrick2");

qdraw(
	ggplot(d, aes(x = x)) + theme_bw() +
		geom_hline(aes(yintercept=0), linetype="dashed", colour="grey30") +
		scale_colour_manual(name="difference",values=cols) +
		geom_line(aes(y=f, colour="truth"), lwd=1) +
		xlab("position") + ylab("group difference") + guides(colour=FALSE)
		#ylim(-2, 3.5)
	,
	height = 2, width = 6,
	file = "truth-diff.pdf"
);

# identify regions where f > 0
idx <- d$f_hat_min > 0;
# mark boundaries
b <- diff(c(0, idx));
starts <- which(b > 0);
ends <- which(b < 0);
if (length(ends) < length(starts)) {
	ends <- c(ends, length(idx));
}
regions <- data.frame(
	start = d$x[starts], end = d$x[ends]
);

qdraw(
	ggplot(d, aes(x = x)) + theme_bw() +
		geom_hline(aes(yintercept=0), linetype="dashed", colour="grey30") +
		scale_colour_manual(name="difference",values=cols) +
		geom_line(aes(y=f, colour="truth"), lwd=1) +
		geom_line(aes(y=f_hat, colour="estimate"), lwd=1) +
		geom_ribbon(aes(ymin = f_hat_min, ymax = f_hat_max), fill=cols["estimate"], alpha=0.3) +
		geom_segment(aes(x=start, xend=end, y=3, yend=3), data=regions, lwd=2) +
		xlab("position") + ylab("group difference") + guides(colour=FALSE)
		#ylim(-2, 3.5)
	,
	height = 2, width = 6,
	file = "estimated-diff.pdf"
);

# write data

saveRDS(data, "gp-compare_sim-data.rds");

