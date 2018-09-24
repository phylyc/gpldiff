library(io)
library(ggplot2)
library(reshape2)

lapply(list.files("../R/", "\\.R$", full.names=TRUE), source);

signal_to_noise2 <- function(data) {
	sd_signal <- sd(c(data$m_b, data$m_a));
	sd_noise <- data$sigma;

	sd_signal^2 / sd_noise^2
}

# NB only works for scalar parameters
assess_bias <- function(B, N, sigma=0.05, pars=c("mu", "sigma2"), ...) {
	names(pars) <- pars;
	biases <- lapply(1:B,
		function(b) {
			data <- rldiff(N, sigma=sigma);
			fit <- gpldiff(data, ...);
			bias(fit$params, data, pars)
		}
	);
	# collect bias estimates across rounds from the same parameter together
	lapply(pars,
		function(p) {
			unlist(lapply(biases,
				function(b) {
					b[[p]]
				}
			))
		}
	)
}

# assess coverage probability
# since f is the only parameter for which we approximate the posterior,
# it is the parameter for which we can derive a coverage probability
assess_coverage <- function(conf.levels, B, N, sigma=0.05, adapt="none", coverage.clevel=0.95, 
	fixed=NULL, ...) {
	coverages <- matrix(unlist(lapply(1:B,
		function(b) {
			data <- rldiff(N, sigma=sigma, ...);
			params <- NULL;
			if (!is.null(fixed)) {
				params <- list(
					mu = 0,
					sigma2 = 1,
					f = rep(0, data$J)
				);
				if (fixed$mu) {
					params$mu <- data$mu;
				}
				if (fixed$sigma2) {
					params$sigma2 <- data$sigma2;
				}
				if (fixed$f) {
					params$f <- data$f;
				}
			}
			fit <- gpldiff(data, params=params, adapt=adapt, fixed=fixed, ...);
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

plot_bias <- function(d) {
	ggplot(d, aes(x=value)) + facet_grid(. ~ L1, scale="free_x") +
		geom_histogram() + theme_bw() + xlab("bias")
}

####


message("SNR = ", snr)
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

d.adapt.snr5 <- assess_coverage(conf.levels, B=100, N=100, sigma=0.1, adapt="Brent");
qdraw(plot_coverage_profile(d.adapt.snr5), file ="coverage_adapt-brent_snr-5.pdf");

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

snrs <- c(10, 5, 2, 1);

lapply(snrs,
	function(snr) {
		set.seed(3);

		N <- 100;
		mu <- 0.5;
		data <- rldiff(N, mu = mu, sigma=mu / snr);
		snr(data)

		fit <- gpldiff(data);
		qdraw({plot(data, fit)}, height=10, file=sprintf("data_nadapt_snr-%s.pdf", snr));
		qdraw({plot(fit, data)}, height=10, file=sprintf("fit_nadapt_snr-%s.pdf", snr));

		fit.adapt <- gpldiff(data, adapt="GD");
		qdraw({plot(data, fit.adapt)}, height=10, file=sprintf("data_adapt_snr-%s.pdf", snr));
		qdraw({plot(fit.adapt, data)}, height=10, file=sprintf("fit_adapt_snr-%s.pdf", snr));

		message("SNR = ", snr)

		message("No adapt")
		print(snr(fit))
		print(coverage(fit, data))
		print(fit$params$sigma)

		message("Adapt")
		print(snr(fit.adapt))
		print(coverage(fit.adapt, data))
		print(fit.adapt$params$sigma)
	}
);

####

set.seed(1);

fixed <- list(
	mu = FALSE,
	sigma2 = TRUE,
	f = FALSE
);

d.nadapt <- assess_coverage(conf.levels, B=100, N=100, fixed=fixed);
qdraw(plot_coverage_profile(d.nadapt), file ="coverage_nadapt_snr-10_fixed-sigma.pdf");

d.adapt <- assess_coverage(conf.levels, B=100, N=100, adapt="GD", fixed=fixed);
qdraw(plot_coverage_profile(d.adapt), file ="coverage_adapt_snr-10_fixed-sigma.pdf");

####

set.seed(1);

fixed <- list(
	mu = FALSE,
	sigma2 = TRUE,
	f = FALSE
);

d.nadapt <- assess_coverage(conf.levels, B=100, N=100, sigma=0.25, fixed=fixed);
qdraw(plot_coverage_profile(d.nadapt), file ="coverage_nadapt_snr-2_fixed-sigma.pdf");

d.adapt <- assess_coverage(conf.levels, B=100, N=100, adapt="GD", sigma=0.25, fixed=fixed);
qdraw(plot_coverage_profile(d.adapt), file ="coverage_adapt_snr-2_fixed-sigma.pdf");

####

set.seed(1);

fixed <- list(
	mu = TRUE,
	sigma2 = FALSE,
	f = FALSE
);

d.nadapt <- assess_coverage(conf.levels, B=100, N=100, fixed=fixed);
qdraw(plot_coverage_profile(d.nadapt), file ="coverage_nadapt_snr-10_fixed-mu.pdf");

d.adapt <- assess_coverage(conf.levels, B=100, N=100, adapt="GD", fixed=fixed);
qdraw(plot_coverage_profile(d.adapt), file ="coverage_adapt_snr-10_fixed-mu.pdf");

####

set.seed(1);

fixed <- list(
	mu = FALSE,
	sigma2 = FALSE,
	f = TRUE
);

d.nadapt <- assess_coverage(conf.levels, B=100, N=100, fixed=fixed);
qdraw(plot_coverage_profile(d.nadapt), file ="coverage_nadapt_snr-10_fixed.pdf");

d.adapt <- assess_coverage(conf.levels, B=100, N=100, adapt="GD", fixed=fixed);
qdraw(plot_coverage_profile(d.adapt), file ="coverage_adapt_snr-10_fixed.pdf");

####

set.seed(1);

b.nadapt <- assess_bias(B=100, N=100);
qdraw(plot_bias(melt(b.nadapt)), width=10, file = "bias_nadapt_snr-10.pdf");

b.adapt <- assess_bias(B=100, N=100, adapt="GD");
qdraw(plot_bias(melt(b.adapt)), width=10, file = "bias_adapt_snr-10.pdf");

b.adapt.brent <- assess_bias(B=100, N=100, adapt="Brent");
qdraw(plot_bias(melt(b.adapt.brent)), width=10, file = "bias_adapt-brent_snr-10.pdf");

####

set.seed(1);

b.nadapt <- assess_bias(B=100, N=100, sigma=0.1);
qdraw(plot_bias(melt(b.nadapt)), width=10, file = "bias_nadapt_snr-5.pdf");

b.adapt <- assess_bias(B=100, N=100, adapt="GD", sigma=0.1);
qdraw(plot_bias(melt(b.adapt)), width=10, file = "bias_adapt_snr-5.pdf");

b.adapt.brent <- assess_bias(B=100, N=100, adapt="Brent", sigma=0.1);
qdraw(plot_bias(melt(b.adapt.brent)), width=10, file = "bias_adapt-brent_snr-5.pdf");

