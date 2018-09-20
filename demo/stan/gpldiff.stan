// Stan implementation of the GPLDIFF model

// NB  Nystrom approximation as implemented below is not
//     well supported by Stan: it can be very slow.

functions {
	
	// time complexity is O(J^3) due to cholesky decomposition of L
	vector calculate_f(matrix K, vector f_eta, int J) {
		// lower triangular matrix
		matrix[J, J] L;

		L = cholesky_decompose(K);
		return L * f_eta;
	}

	// time complexity is O(J M^2) due to matrix multiplication for B
	vector approximate_f(matrix K, vector f_eta, vector u, int J, int M) {
		// subset index
		int idx[M];

		// inverse of subset kernel matrix
		matrix[M, M] W;

		idx = sort_indices_asc(u)[1:M];

		// approximate K as \tilde{K} = C W C^T 
		// where
		//   C \in R^{J \times M} is subset of columns
		//   G \in R^{M \times M} is subset of columns and rows
		//   W = G^{-1}
		W = inverse_spd(K[idx, idx]);

		// define B \in R^{J \times M} s.t. \tilde{K} = B B^T
		//   B = C * cholesky_decompose(W)
		// then
		//   f = B * f_eta
		return K[, idx] * (cholesky_decompose(W) * f_eta);
	}

}

data {
	// number of loci
	int<lower=1> J;

	// number of data points to subset for Nystr\"om approximation
	int<lower=1,upper=J> M;

	// positions
	real x[J];

	// scores
	vector[J] y;

	// group membership
	vector[J] g;

	// hyperparameters
	
	// parameters on the prior of observation variance sigma2
	real<lower=0> alpha;
	real<lower=0> beta;

	// variance of the observation mean mu
	real<lower=0> tau2;
}

parameters {
	// covariance function parameters
	// scale parameter of the kernel
	real<lower=0> nu;
	// characteristic length scale
	real<lower=0> lambda;

	// mean of the observation model
	real mu;
	// variance of the observation model
	real<lower=0> sigma2;

	// dummy variable for constructing the Gaussian process
	vector[M] f_eta;

	// dummy variable for subsampling
	vector<lower=0,upper=1>[J] u;
}

transformed parameters {
	// latent function values sampled from the Gaussian process
	// main parameter of interest
	vector[J] f;

	// define auxiliary parameters inside private scope
	{
		// kernel matrix
		matrix[J, J] K;

		// evaluate covariance function
		K = cov_exp_quad(x, nu, lambda);
		// add slight regularization for stability
		// otherwise, we would encounter:
		//   Error evaluating the log probability at the initial value
		for (j in 1:J) {
			K[j, j] = K[j, j] + 1e-5;
		}
		
		if (J == M) {
			f = calculate_f(K, f_eta, J);
		} else {
			f = approximate_f(K, f_eta, u, J, M);
		}
	}
}

model {
  y ~ normal(mu + (f .* g), sqrt(sigma2));

	mu ~ normal(0, sqrt(tau2));
  sigma2 ~ inv_gamma(alpha, beta);

  f_eta ~ normal(0, 1);

	// for Nystr\"om approximation; only used if M < J
	u ~ uniform(0, 1);
}

