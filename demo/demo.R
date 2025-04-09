library(io)
library(ggplot2)

lapply(list.files("../R/", "\\.R$", full.names=TRUE), source);

set.seed(5);
# simulate two non-linear data series with missingness
data <- rldiff(200);
f <- data$f;

# write data
saveRDS(data, "gp-compare_sim-data.rds");

# hyperparameters
hparams <- list(
	nu2 = 1,
	lambda2 = 2,
	alpha = 0.1,
	beta = 0.1,
	tau2 = 10
);

# initial guess values of parameters
params <- list(
	mu = 0,
	sigma2 = 1,
	f = rep(0, data$J)
);

fixed <- list(
	mu = FALSE,
	sigma2 = FALSE,
	f = FALSE
);

# ---

# sigma2 needs to be overestimated for better fit... why?
# set to 0.1^2 = 0.01 but estimate is ~0.3
# if sigma2 is constrained to known value, then 
# the fit of f is much poorer
#
#params <- list(
#	mu = 0,
#	sigma2 = sigma^2,
#	f = rep(0, data$J)
#);
#
#fixed <- list(
#	mu = FALSE,
#	sigma2 = TRUE,
#	f = FALSE
#);


fit.nadapt <- gpldiff(data, params=params, hparams=hparams, fixed=fixed);

# Optimize hyperparameters
fit.adapt <- gpldiff(data, adapt="GD", tol2=1e-2, max.iter2=50, learn.rate=0.2, verbose=TRUE);
#fit.adapt2 <- gpldiff(data, adapt="Brent", tol2=1e-2, max.iter2=20, learn.rate=0, verbose=TRUE);
#fit.adapt3 <- gpldiff(data, adapt="L-BFGS-B", tol2=1e-2, max.iter2=20, learn.rate=0, verbose=TRUE);

fit.nadapt$evidence
fit.adapt$evidence

str(fit.nadapt)
str(fit.adapt)

cor(f, fit.nadapt$params$f)
cor(f, fit.adapt$params$f)

rmse(f, fit.nadapt$params$f)
rmse(f, fit.adapt$params$f)

# plot data
qdraw(plot(data, fit.nadapt), height=10, "data_no-adapt.pdf");
qdraw(plot(data, fit.adapt), height=10, "data_adapt.pdf");

# plot fit
qdraw(plot(fit.nadapt, data), height=10, "fit_no-adapt.pdf");
qdraw(plot(fit.adapt, data), height=10, "fit_adapt.pdf");

cov.levels <- c(0.99, 0.95, 0.9, 0.8, 0.7, 0.6, 0.5);
cov.nadapt <- unlist(lapply(cov.levels, function(level) coverage(fit.nadapt, f, level)));
cov.adapt <- unlist(lapply(cov.levels, function(level) coverage(fit.adapt, f, level)));

plot(cov.levels, cov.nadapt, xlim=c(0.5, 1), ylim=c(0.5, 1))
abline(a=0, b=1)

plot(cov.levels, cov.adapt, xlim=c(0.5, 1), ylim=c(0.5, 1))
abline(a=0, b=1)

# prepare data for ggplot

fit <- fit.adapt;

gcols <- c(case="#EFC000FF", control="#0073C2FF");

truth <- with(data, data.frame(
	x = c(x, x), y = c(m_a, m_b),
	group = factor(
		c(rep("control", length(x)), rep("case", length(x))),
		levels=names(gcols)
	)
));

# calculate confidence interval of f estimate
cf <- confint(fit);

d <- with(data, data.frame(
	x = x,
	f = f,
	f_hat = fit$params$f,
	f_hat_min = cf$lower,
	f_hat_max = cf$upper,
	y = y,
	y_a = y_a,
	y_b = y_b,
	y_missing = ifelse(g <= 0, y_b, y_a),
	group = factor(g, levels=c(0.5, -0.5), labels=names(gcols)),
	group_missing = factor(g, levels=c(-0.5, 0.5), labels=names(gcols))
));

# identify regions where f > 0
regions <- find_sig_regions(fit, data);


# make prettier plots

ymin <- -0.30;
ymax <- 1.25;

qdraw(
	ggplot(truth, aes(x, y, colour=group)) + theme_bw() +
		scale_colour_manual(values=gcols) +
		geom_line(lwd=1) +
		xlab("position") + ylab("response") +
		ylim(ymin, ymax)
	,
	height = 2, width = 6,
	file = "truth-data.pdf"
);

qdraw(
	ggplot(d, aes(x, y, colour=group)) + theme_bw() +
		scale_colour_manual(values=gcols) +
		geom_point() +
		geom_spoke(aes(y=-0.30, colour=group), angle=pi/2, radius=0.05) +
		xlab("position") + ylab("response") +
		ylim(ymin, ymax)
	,
	height = 2, width = 6,
	file = "observed-data.pdf"
);

qdraw(
	ggplot(d, aes(x, y, colour=group)) + theme_bw() +
		scale_colour_manual(values=gcols) +
		geom_segment(aes(x=x, xend=x, y=y_a, yend=y_b), colour="grey", alpha=0.5) +
		geom_point(aes(shape="yes")) +
		geom_point(aes(y=y_missing, shape="no", colour=group_missing)) +
		scale_shape_manual(name="observed", values=c(yes=19, no=21)) +
		xlab("position") + ylab("response") +
		ylim(ymin, ymax)
	,
	height = 2, width = 6,
	file = "hidden-data.pdf"
);

cols <- c(truth="black", estimate="#CD534CFF");

ymin <- -1.5;
ymax <- 1.5;

qdraw(
	ggplot(d, aes(x = x)) + theme_bw() +
		geom_hline(aes(yintercept=0), linetype="dashed", colour="grey30") +
		scale_colour_manual(name="difference",values=cols) +
		geom_line(aes(y=f, colour="truth"), lwd=1) +
		xlab("position") + ylab("group difference") + guides(colour=FALSE) +
		ylim(ymin, ymax)
	,
	height = 2, width = 6,
	file = "truth-diff.pdf"
);

qdraw(
	ggplot(d, aes(x = x)) + theme_bw() +
		geom_hline(aes(yintercept=0), linetype="dashed", colour="grey30") +
		scale_colour_manual(name="difference",values=cols) +
		geom_line(aes(y=f, colour="truth"), lwd=1) +
		geom_line(aes(y=f_hat, colour="estimate"), lwd=1) +
		geom_ribbon(aes(ymin = f_hat_min, ymax = f_hat_max), fill=cols["estimate"], alpha=0.3) +
		geom_segment(aes(x=start, xend=end, y=ymax*0.9, yend=ymax*0.9, alpha=posterior), data=regions, lwd=2) +
		guides(alpha=FALSE) +
		xlab("position") + ylab("group difference") + guides(colour=FALSE) +
		ylim(ymin, ymax)
	,
	height = 2, width = 6,
	file = "estimated-diff.pdf"
);

