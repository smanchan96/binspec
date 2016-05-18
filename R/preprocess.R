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
#' Create a binary matrix, each column represents and m/z value, and each row represents a mass spectra.  The value indicates whether or not the m/z of this spectra is a peak.
#' @param list_mz_peaks List of m/z peak vectors
#' @export combine_peaks

combine_peaks <- function(list_mz_peaks)  {
    NULL
}

#' Random Forest
#'
#' Build a random forest classifier and test it.  No need for test set because out-of-bag error measurement.
#' @param training_set
#' @param training_labels
#' @export random_forest_classify

#' Artificial Neural Network
#'
#' Build a ANN classifier and test it.
#' @param training_set
#' @param training_labels
#' @param test_set
#' @param test_labels
#' @export ann_classify

#' Support Vector Machine
#'
#' Build a SVM classifier and test it.
#' @param training_set
#' @param training_labels
#' @param test_set
#' @param test_labels
#' @export svm_classify

