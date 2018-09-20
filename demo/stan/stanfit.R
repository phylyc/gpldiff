library(rstan);

# prerequisite:
#   generate `../gp-compare_sim-data.rds`
#   by running `../simulate.R`

# approximate time: each chain can take 1-2 hours

# Notes
# Divergent transitions were detected after warm-up.
# Increasing `adapt_delta` to 0.99 (maximum is 1) drastically
# increased the run time but did not reduce divergent transitions
# by much.

# use multiple cores to sample each chain in parallel
options(mc.cores=4);


# root mean square error
rmse <- function(x, y) {
	sqrt(mean((x - y)^2))
}

# Calculate observed coverage probability
coverage <- function(x, sample, alpha=0.05) {
	bounds <- t(apply(sample, 2, quantile, probs=c(alpha/2, 1 - alpha/2)));
	mean(bounds[,1] <= x & x <= bounds[,2])
}


data <- readRDS("../gpldiff_sim-data.rds");

# do not use Nystr\"om approximation
data$M <- data$J;

# set fixed hyperparameters
data$alpha <- 0.1;
data$beta <- 0.1;
data$tau2 <- 10;

# hyperparameters nu2 and lambda2 will be learned

fit <- stan(
	file="gpldiff.stan",
	data=data, iter=2000, chains=4,
	pars=c("f", "sigma2", "nu", "lambda"),
	control = list(adapt_delta = 0.99)
);
saveRDS(fit, "stanfit.rds")

print(fit)

# Rhat > 1 for many parameters:
# We will need much longer chains to represent the posterior well

f.mcmc <- extract(fit, "f")[[1]];
f.mcmc.mean <- apply(f.mcmc, 2, mean);
f.mcmc.median <- apply(f.mcmc, 2, median);

cor(data$f, f.mcmc.mean)
rmse(data$f, f.mcmc.mean)

cor(data$f, f.mcmc.median)
rmse(data$f, f.mcmc.median)

coverage(data$f, f.mcmc, alpha=0.05)
coverage(data$f, f.mcmc, alpha=0.1)
coverage(data$f, f.mcmc, alpha=0.25)
coverage(data$f, f.mcmc, alpha=0.5)

