#' getPrComps
#'
#' Get principal components that explain 90% of variance.
#' @param spectrums nxd matrix with n vectors of d dimensions
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

#' kmeans CH
#'
#' Run k-means, using the k resulting in best CH value
#' @param nxd matrix with n individuals and d dimensions
#' @export kmeansCH
#' @examples
#' kmeansCH()
kmeansCH <- function(f) {
	require(vegan)
	fit <- cascadeKM(f, 2, 6)
	k <- which.max(fit$results[2,])+1
	kmeans(f, k)
}

#' forestSelect
#'
#' Get top r features using data matrix and labels for data
#' @param nxd matrix with n individuals and d dimensions
#' @param vector of length n with integer labels
#' @param number of features to select (r<d)
#' @export forestSelect
#' @examples
#' forestSelect()
forestSelect <- function(m, lbl, r) {
  require(randomForest)
  lbl <- as.factor(lbl)
  myforest <- randomForest(m, lbl, importance = T)
  head(sort(myforest$importance[,7], decreasing = T,)n=r)
}