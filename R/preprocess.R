#' Round data frame
#'
#' Round all m/z and intensity values to integers.
#' @param df Data frame
#' @export round_df

round_df <- function(df) {
    #Round each column to the nearest integer
    data.frame(apply(df, 2, function(x) sapply(x, round)))
}

#' Find binary peaks
#'
#' Find peaks in window of size 2*neighbors + 1 and label m/z integers within the error as peaks.  Returns vector of peak m/z integers.
#' @param df Data frame of m/z and intensities
#' @param neighbors Number of neighboring m/z values to compare on right and left
#' @param error m/z Decimal error value
#' @export binary_peaks

binary_peaks <- function(df, neighbors, error=0) {
    #For each m/z value, see if it is max of its neighbors
    #peaks are true/false values
    peaks <- sapply(1:nrow(df), function(x) {
                    window <- intersect((x - neighbors):(x + neighbors), 1:nrow(df))
                    all(df[[2]][x] >= df[[2]][window])
}
    )
    #All the m/z values within the error of a peak should also be labled as peaks
    peakmzs <- df[[1]][peaks]
    #newpeaks are all o f
    toReturn <- unique(unlist(lapply(peakmzs, function(x) {
                       round(x - error/2 * x):round(x + error/2 * x) 

})))
    #remove items that are not in original m/z range (assumes df[[1]] sorted)
    toReturn[toReturn >= df[[1]][1] & toReturn  <= df[[1]][nrow(df)]]
}

#' Combine peak vectors
#'
#' Create a binary matrix, each column represents an m/z value, and each row represents a mass spectra.  The value indicates whether or not the m/z of this spectra is a peak.
#' @param list_mz_peaks List of m/z peak vectors
#' @export combine_peaks

combine_peaks <- function(list_mz_peaks)  {
    mzs <- sort(unique(unlist(list_mz_peaks)))
    mymatrix <- do.call(rbind, lapply(list_mz_peaks, function(x) {
                          mzs %in% x
}))
    colnames(mymatrix) <- mzs
    mymatrix
}

#' Classifier Accuracies
#'
#' Find the best classifier using leave-one-out cross validation (svm) and out-of-bag error (random forests).  Returns a vector of classification accuracies
#' @param peaks Boolean matrix of mass spectra rows with m/z columns, indicating if an m/z value corresponds to a peak.
#' @param labels The correct classifications of the peaks.
#' @param minpeaks How many "true" values must show up for a given m/z value for it to be considered a feature.
#' @export classifier_accuracies

classifier_accuracies <- function(peaks, labels, minpeaks) {
  peaks <- peaks[,(apply(peaks, 2, sum)) > minpeaks]
  optsvm <- tune.svm(peaks, labels, kernel = "radial", cost = sapply(seq(-5, 15, by=2), function(x)2^x), gamma = sapply(seq(-15, 3, by=2), function(x)2^x), tunecontrol=tune.control(nrow(peaks)))
  result <- svm(peaks, labels, kernel="radial", gamma=optsvm$best.parameters$gamma, cost=optsvm$best.parameters$cost, cross=nrow(peaks))
  # resrf <- randomForest(peaks, labels, ntree = 1000, mtry = 9)
  optimalRF <- tuneRF(peaks, labels, doBest = T)
  toReturn <- c(sum(labels==optimalRF$predicted)/length(labels), sum(result$accuracies)/100/length(labels))
  names(toReturn) <- c("RF", "SVM")
  toReturn
}
