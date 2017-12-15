region_center <- function(start, end) {
	start + floor((end - start)/2)
}

#' @param split  whether to split chromosome into chromosome arms
read_gistic <- function(fname, split=TRUE, genome="hg19") {
	x <- read.table(fname, sep="\t", header=TRUE);
	colnames(x) <- c("type", "chromosome", "start", "end", "nlq", "g_score", "mean_amplitude", "frequency");
	# calculate position as the center of each segment
	x$position <- with(x, region_center(start, end));

	if (split) {
		if (genome != "hg19") {
			stop("Only hg19 is currently supported");
		}

		# split chromosomes
		centromeres$hg19
	}

	x
}

prepare_gistics <- function(case, control) {
	# concatenate control data with case data
	data <- list(
		# total sample size
		J = nrow(control) + nrow(case),
		# convert position to Mbp
		x = c(control$position, case$position) / 1e6,
		# group membership
		g = c(rep(-0.5, nrow(control)), rep(0.5, nrow(case))),
		y = c(control$g_score, case$g_score)
	);

	# data are currently sorted by cohort but need to be sorted by position
	idx <- order(data$x);
	data$x <- data$x[idx];
	data$g <- data$g[idx];
	data$y <- data$y[idx];

	data
}

#' Compare GISTIC scores using GPLDIFF.
#'
#' @param case     file name of GISTIC scores table for case cohort
#' @param control  file name of GISTIC scores table for control cohort
#' @param param    initial parameter values to \code{gpldiff()}
#' @param hparams  hyperparameter values to \code{gpldiff()}
#' @param ...      other paramsters to \code{gpldiff()}
#' @return a list of \code{gpldiff} objects
compare_gistics <- function(case, control, params=NULL, hparams=NULL, adapt=FALSE, ...) {

	if (is.character(case)) {
		case <- read_gistic(case);
	}

	if (is.character(control)) {
		control <- read_gistic(control);
	}

	if (is.null(hparams)) {
		hparams <- list(
			nu2 = 5^2,
			lambda2 = 1^2,
			alpha = 10,
			beta = 1,
			tau2 = 1^2
		);
	}

	case.split <- split(case, list(type = case$type, chromosome = case$chromosome));
	control.split <- split(control, list(type = control$type, chromosome = control$chromosome));

	data.sets <- mcmapply(prepare_gistics, case.split, control.split, SIMPLIFY=FALSE);
	models <- mclapply(names(data.sets), function(i) { 
		message("Processing ", i)
		gpldiff(data.sets[[i]], params=params, hparams=hparams, adapt=adapt, ...)
	});

	fits <- mapply(function(d, m) list(data = d, model = m), data.sets, models, SIMPLIFY=FALSE)
	
	structure(fits, class="gistic_gpldiffs")
}

#' Summarize \code{gistic_gpldiffs} object.
#' @param fits  a list of \code{gpldiff} objects
#'              (returned from \code{compare_gistics})
#' @return list of data.frame
summary.gistic_gpldiffs <- function(fits) {
	regions.all <- lapply(
		fits,
		function(fit) {
			find_sig_regions(fit$model, fit$data)
		}
	);

	regions.amp <- regions.all[grep("Amp", names(regions.all))];
	names(regions.amp) <- sub("Amp.", "", names(regions.amp));
	regions.amp <- process_regions(combine_regions(regions.amp));

	regions.del <- regions.all[grep("Del", names(regions.all))];
	names(regions.del) <- sub("Del.", "", names(regions.del));
	regions.del <- process_regions(combine_regions(regions.del));

	rbind(
		data.frame(
			type = "Amp",
			regions.amp
		),
		data.frame(
			type = "Del",
			regions.del
		)
	)
}

