# Calculate the determinant of a matrix
det <- function(x, log=FALSE, ...) {
	z <- determinant(x, logarithm=TRUE, ...)
	if (log) {
		c(z$modulus)
	} else {
		c(z$sign * exp(z$modulus))
	}
}
