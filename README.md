# Gaussian process latent difference model #

Package for fitting a Gaussian process model to detect latent difference between two populations, under the assumption that this latent difference is distributed as a multivariate Gaussian.

### Installation ###

The required documentation files for this package needs to generated before
this package can be installed using `devtools`.

```
library(devtools)
document()
install()
```

Note: Do not simply run devtools::install_bitbucket() as it will not generated the documentation files.

### Example ###

The observed `data` should be inputted as a `list` consisting of

- `x`, a vector containing values for the continuous independent variable
- `g`, a vector containing `-0.5` or `0.5` for membership in control group or case group, respectively
- `y`, a vector containing values for the continuous response variable

The model hyperparameters should be fixed, while the model parameters may be inferred by the model from the data.

```
library(gpldiff)

hparams <- list(
	nu2 = 1,
	lambda2 = 1,
	alpha = 2,
	beta = 1,
	tau2 = 1
)

params <- NULL
```

We can now fit the model by

```
fit <- gpldiff(data, params, hparams)
```

We can view the summary statistics and generate plots:

```
summary(fit)
plot(fit, data)
```

See `?gpldiff` for more information.
