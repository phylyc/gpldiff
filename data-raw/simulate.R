# simulate non-linear data of two groups with missingness

set.seed(1);

# generate data

N <- 100;
N <- 500;
sigma <- 0.1;

x <- sort(seq(0, 4*pi, length.out=N) + rnorm(N, sd = 0.1));

m_a <- sin(x);
m_b <- sin(x/(pi/2) - 7*pi/2) + 0.5;

y_a <- rnorm(N, mean = m_a, sd = sigma);
y_b <- rnorm(N, mean = m_b, sd = sigma);

g <- sample(c(-0.5, 0.5), N, replace=TRUE);
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


# plot data

par(mfrow=c(4, 1));

plot(c(x, x), c(m_a, m_b), xlab="x", ylab="y", type="n", main = "original data");
lines(x, m_a, col=3);
lines(x, m_b, col=4);

plot(c(x, x), c(y_a, y_b), xlab="x", ylab="y", type="n", main = "data with noise");
points(x, y_a, col=3);
points(x, y_b, col=4);

plot(x, y, col=as.numeric(g > 0)+3, main = "data with noise and missingness");

plot(x, f, main = "difference in mean");


# write data

saveRDS(data, "gp-compare_sim-data.rds");
saveRDS(data, "gp-compare_sim-data_med.rds");

