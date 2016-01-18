#' getPrComps
#'
#' Get principal components that explain 90% of variance.
#' @param spectrums nxm matrix representing frequencies for n individuals over m integer m/z values.
#' @export getPrComps
#' @examples
#' getPrComps()

getPrComps <- function(spectrums) {
	require(stats)
	result <- prcomp(spectrums)
	variances <- sapply(result$sdev, function(x)x^2)
	sum <- 0
	totake <- 1
	for (i in 1:length(variances)) {
		if (sum > .9) {
			totake <- i-1
			break
		}
		sum <- sum + variances[i]
	}
	return(result$x[,1:totake])
}

