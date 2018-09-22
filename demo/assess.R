library(io)
library(ggplot2)

lapply(list.files("../R/", "\\.R$", full.names=TRUE), source);

signal_to_noise2 <- function(data) {
	sd_signal <- sd(c(data$m_b, data$m_a));
	sd_noise <- data$sigma;

	sd_signal^2 / sd_noise^2
}

assess_coverage <- function(conf.levels, B, N, sigma=0.05, adapt="none", coverage.clevel=0.95, ...) {
	coverages <- matrix(unlist(lapply(1:B,
		function(b) {
			data <- rldiff(N, sigma=sigma, ...);
			fit <- gpldiff(data, adapt=adapt, ...);
			unlist(lapply(conf.levels, function(cl) coverage(fit, data, cl)))
		}
	)), nrow=length(conf.levels));
	
	coverages_to_df(coverages, coverage.clevel=coverage.clevel)
}

coverages_to_df <- function(coverages, coverage.clevel=0.95) {
	cov.means <- rowMeans(coverages);
	cov.lowers <- apply(coverages, 1, quantile, prob=1-coverage.clevel);
	cov.uppers <- apply(coverages, 1, quantile, prob=coverage.clevel);

	data.frame(
		confidence_level = conf.levels,
		coverage = cov.means,
		coverage_min = cov.lowers,
		coverage_max = cov.uppers
	)
}

plot_coverage_profile <- function(d) {
	ggplot(d, aes(x=confidence_level, y=coverage, ymin=coverage_min, ymax=coverage_max)) + 
		theme_bw() +
		geom_abline(slope=1, size=2, colour="grey90") + geom_point() +
		geom_errorbar(width=0.01, colour="grey60") +
		xlim(0, 1.01) + ylim(0, 1.01)
}

####

set.seed(1);

conf.levels <- c(0.99, 0.95, 0.90, 0.80, 0.70, 0.60, 0.50, 0.40, 0.30, 0.20, 0.10);

d.nadapt <- assess_coverage(conf.levels, B=100, N=100);
qdraw(plot_coverage_profile(d.nadapt), file ="coverage_nadapt_snr-10.pdf");

d.adapt <- assess_coverage(conf.levels, B=100, N=100, adapt="GD");
qdraw(plot_coverage_profile(d.adapt), file ="coverage_adapt_snr-10.pdf");

####

set.seed(1);

d.nadapt.snr5 <- assess_coverage(conf.levels, B=100, N=100, sigma=0.1);
qdraw(plot_coverage_profile(d.nadapt.snr5), file ="coverage_nadapt_snr-5.pdf");

d.adapt.snr5 <- assess_coverage(conf.levels, B=100, N=100, sigma=0.1, adapt="GD");
qdraw(plot_coverage_profile(d.adapt.snr5), file ="coverage_adapt_snr-5.pdf");

####

set.seed(1);

d.nadapt.snr2 <- assess_coverage(conf.levels, B=100, N=100, sigma=0.25);
qdraw(plot_coverage_profile(d.nadapt.snr2), file ="coverage_nadapt_snr-2.pdf");

d.adapt.snr2 <- assess_coverage(conf.levels, B=100, N=100, sigma=0.25, adapt="GD");
qdraw(plot_coverage_profile(d.adapt.snr2), file ="coverage_adapt_snr-2.pdf");

####

set.seed(1);

d.nadapt.snr1 <- assess_coverage(conf.levels, B=100, N=100, sigma=0.5);
qdraw(plot_coverage_profile(d.nadapt.snr1), file ="coverage_nadapt_snr-1.pdf");

d.adapt.snr1 <- assess_coverage(conf.levels, B=100, N=100, adapt="GD", sigma=0.5);
qdraw(plot_coverage_profile(d.adapt.snr1), file ="coverage_adapt_snr-1.pdf");

####

set.seed(1);

N <- 100;
data <- rldiff(N, sigma=0.05);
snr(data)

fit <- gpldiff(data);
coverage(fit, data)
fit$params$sigma
snr(fit)
qdraw({plot(data, fit)}, height=10, file="fit_nadapt_snr-10.pdf")

fit.adapt <- gpldiff(data, adapt="GD");
coverage(fit.adapt, data)
fit.adapt$params$sigma
snr(fit.adapt)
qdraw({plot(data, fit.adapt)}, height=10, file="fit_adapt_snr-10.pdf");

####

set.seed(1);

N <- 100;
data <- rldiff(N, sigma=0.10);
snr(data)

fit <- gpldiff(data);
coverage(fit, data)
fit$params$sigma
snr(fit)
qdraw({plot(data, fit)}, height=10, file="fit_nadapt_snr-5.pdf")

fit.adapt <- gpldiff(data, adapt="GD");
coverage(fit.adapt, data)
fit.adapt$params$sigma
snr(fit.adapt)
qdraw({plot(data, fit.adapt)}, height=10, file="fit_adapt-gd_snr-5.pdf");

fit.adapt <- gpldiff(data, adapt="Brent");
coverage(fit.adapt, data)
fit.adapt$params$sigma
snr(fit.adapt)
qdraw({plot(data, fit.adapt)}, height=10, file="fit_adapt-brent_snr-5.pdf");

####

set.seed(1);

N <- 100;
data <- rldiff(N, sigma=0.25);
snr(data)

fit <- gpldiff(data);
coverage(fit, data)
fit$params$sigma
snr(fit)
qdraw({plot(data, fit)}, height=10, file="fit_nadapt_snr-2.pdf")

fit.adapt <- gpldiff(data, adapt="GD");
coverage(fit.adapt, data)
fit.adapt$params$sigma
snr(fit.adapt)
qdraw({plot(data, fit.adapt)}, height=10, file="fit_adapt_snr-2.pdf");

####

set.seed(1);

N <- 100;
data <- rldiff(N, sigma=0.5);
snr(data)

fit <- gpldiff(data);
coverage(fit, data)
fit$params$sigma
snr(fit)
qdraw({plot(data, fit)}, height=10, file="fit_nadapt_snr-1.pdf")

fit.adapt <- gpldiff(data, adapt="GD");
coverage(fit.adapt, data)
fit.adapt$params$sigma
snr(fit.adapt)
qdraw({plot(data, fit.adapt)}, height=10, file="fit_adapt_snr-1.pdf");

####

