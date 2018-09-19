#' @param model  \code{gpldiff} model
#' @param data   data object to which \code{gpldiff} model was fitted
#' @param lodds.cut   threshold for log odds for determining candidate regions
#' @param max.gap     if the gap size between two adjacent candidate regions
#'                    is less than this threshold, then these regions are
#'                    merged together
#' @param min.obs     minimum number of observations required for any significant region
#' @param direction   if \code{direction > 0}, test for \code{f > 0};
#'                    otherwise, test for \code{f <= 0}
#' @return \code{data.frame} of significant regions described by fields:
#' \itemize{
#'    \item \code{start}, x value at the start of region
#'    \item \code{end}, x value at the end of region
#'    \item \code{start_idx}, index at start of region
#'    \item \code{end_idx}, index at end of region
#'    \item \code{n_obs}, number of observations within region
#'    \item \code{diff}, observed difference in mean of data points within
#'          region between case and control
#'    \item \code{ldiff}, estimated latent difference \code{f}
#'    \item \code{posterior}, posterior probability that latent difference
#'          \code{f > 0} (or \code{f <= 0} if \code{direction <= 0})
#' }
#' start, x value at start of region;
#' end, 
find_sig_regions <- function(model, data, lodds.cut=2, max.gap=5, min.obs=2, direction=1) {
	# find candidate regions

	prob <- summary(model);
	if (direction <= 0) {
		prob <- 1 - prob;
	}

	lodds <- log10(prob) - log10(1 - prob);
	idx <- which(lodds > lodds.cut);

	if (length(idx) <= 1) return(NULL);

	# mark contiguous start and ends with 1 and -1 respectively
	# gap size between two contiguous regions is j - i - 1
	# where i is the end index of the first region
	#       j is the start index of the second region
	# therefore, j - i - 1 <= max.gap
	#            j - i <= max.gap + 1
	# note the start (1) and end (-1) markers inserted at both ends
	boundaries <- c(1, diff(diff(idx) <= max.gap + 1), -1);

	start_idx <- idx[which(boundaries == 1)];
	end_idx <- idx[which(boundaries == -1)];

	if (length(start_idx) >= 2 && end_idx[1] > start_idx[2]) {
		# second region starts before the first segment ends
		# first region is a singleton and not marked by an end marker
		start_idx <- start_idx[-1];
	}

	if (length(end_idx) >= 2 && start_idx[length(start_idx)] < end_idx[length(end_idx)-1]) {
		# last region is a singleton and not marked by a start marker
		end_idx <- end_idx[-length(end_idx)];
	}

	# construct start and end indices of candidate regions
	regions <- data.frame(
		# convert position from Mbp back to bp
		# FIXME this only makes sense for GISTIC!
		start = data$x[start_idx] * 1e6,
		end = data$x[end_idx] * 1e6,
		start_idx = start_idx,
		end_idx = end_idx,
		n_obs = end_idx - start_idx + 1
	);

	# filter regions for number of observations
	regions <- regions[regions$n_obs >= min.obs, ];

	if (nrow(regions) > 0) {
		# TODO avoid duplication

		# calculate observed difference in means
		regions$diff <- unlist(lapply(
			1:nrow(regions),
			function(ri) {
				ridx <- regions$start_idx[ri]:regions$end_idx[ri];
				y <- data$y[ridx];
				g <- data$g[ridx];
				mean(y[g > 0]) - mean(y[g < 0])
			}
		));

		# calculate mean latent difference
		regions$ldiff <- unlist(lapply(
			1:nrow(regions),
			function(ri) {
				ridx <- regions$start_idx[ri]:regions$end_idx[ri];
				mean(model$params$f[ridx])
			}
		));

		# calculate posterior probabilies
		regions$posterior <- unlist(lapply(
			1:nrow(regions),
			function(ri) {
				ridx <- regions$start_idx[ri]:regions$end_idx[ri];
				p <- region_posterior(ridx, data, model)
				if (direction > 0) {
					p
				} else {
					1 - p
				}
			}
		));
	}

	regions
}

# calculate the posterior probability that mean f in region > 0
region_posterior <- function(ridx, data, model) {
	n <- length(ridx);
	mean_region <- mean(model$params$f[ridx]);
	# TODO compute full covariance matrix once and subset it as needed?
	covar_region <- f_laplace_covariance(data$g[ridx], model$params$sigma2, model$predict$K[ridx, ridx]);
	se_region <- sqrt(sum(covar_region) / n);

	# Pr(mean_region > 0)
	posterior<- pnorm(0, mean=mean_region, sd=se_region, lower.tail=FALSE);
	posterior
}

# check whether two regions overlap
overlap <- function(s1, e1, s2, e2) {
	# use bitwise OR to allow all arguments to be vectors
	!( e1 < s2 | e2 < s1 )
}

# filter out regions that overlap with padded centromere regions
filter_centromere_regions <- function(regions, padding=10e6, genome="hg19") {
	cen_chroms <- centromeres[[genome]]$chromosome;
	cen_starts <- centromeres[[genome]]$start - padding;
	cen_ends <- centromeres[[genome]]$end + padding;

	idx <- match(regions$chromosome, cen_chroms);
	regions[!overlap(regions$start, regions$end, cen_starts[idx], cen_ends[idx]), ]
}



# combine regions from different chromosomes together
combine_regions <- function(regions) {
	combined <- do.call(rbind,
		mapply(
			function(d, chrom) {
				if (!is.null(d)) {
					data.frame(
						chromosome = chrom,
						d
					)
				} else {
					NULL
				}
			},
			regions,
			names(regions),
			SIMPLIFY = FALSE
		)
	);
	rownames(combined) <- NULL;

	combined
}

# process regions of the same type
process_regions <- function(regions, direction=1, ...) {
	regions <- regions[order(regions$posterior, decreasing=TRUE), ];

	# apply non-secific filters
	regions.f <- filter_centromere_regions(regions, ...);
	# filter problematic regions
	# regions with NaN diff are usually spurious
	if (direction > 0) {
		# regions with negative diff are contradictory... bug or numeric instability in the code???
		regions.f <- regions.f[is.finite(regions.f$diff) & regions.f$diff > 0, ]
	} else {
		regions.f <- regions.f[is.finite(regions.f$diff) & regions.f$diff < 0, ]
	}

	# calculate Bayesian FDR
	regions.f$fdr <- bayesian_fdr(regions.f$posterior);

	regions.f
}

