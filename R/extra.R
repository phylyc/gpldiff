coordinate_to_string <- function(x) {
	sprintf("%s:%d-%d",
		x$chromosome,
		x$start,
		x$end
	)
}
