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
#' Find the best classifier using leave-one-out cross validation (svm) and out-of-bag error (random forests).  Returns a list of classifier results
#' @param peaks Boolean matrix of mass spectra rows with m/z columns, indicating if an m/z value corresponds to a peak.
#' @param labels The correct classifications of the peaks.
#' @param training The rows to actually use to train
#' @param minpeaks How many "true" values must show up for a given m/z value for it to be considered a feature.
#' @export classifier_accuracies

classifier_accuracies <- function(peaks, labels, training, min_peak_percentage) {
  require("e1071")
  require("randomForest")
  minpeaks <- floor(min_peak_percentage * length(labels))
  peaks <- peaks[,(apply(peaks, 2, sum)) > minpeaks]
  peaks <- peaks[training,]
  labels <- labels[training]
  optsvm <- tune.svm(peaks, labels, kernel = "radial", cost = sapply(seq(-5, 5, by=2), function(x)2^x), gamma = sapply(seq(-9, 3, by=2), function(x)2^x), tunecontrol=tune.control(nrow(peaks)))
  svm_result <- svm(peaks, labels, kernel="radial", gamma=optsvm$best.parameters$gamma, cost=optsvm$best.parameters$cost, cross=nrow(peaks))
  optimalRF <- tuneRF(peaks, labels, doBest = T)
  #toReturn <- c(sum(labels==optimalRF$predicted)/length(labels), sum(result$accuracies)/100/length(labels))
  list(RF=optimalRF, SVM=svm_result)
}

#' SVM and RF Accuracies
#' 
#' Given a vector of neighbor values and a vector of the minimum number of peaks to be considered, this function finds the peak mz values for a data set by running binary_peaks using each of the neighbor vector values, runs SVM and RF on the peaks for each of the min_peak_count values, and returns the accuracies of each test in a table. The table's rows are the number of neighbors, and the columns are the min_peak_count values.
#' @param list_of_dfs The first data frame of mz values and frequencies
#' @param labels The labels of the two states the first data frame's values could be classified as
#' @param neighbors A vector of the number of neighbors to be considered in the binary_peaks function
#' @param min_peaks_percentage A vector of the minimum percent of times an m/z must be a peak to be considered in the classifier_accuracies function
#' @param errow_window A vector of percentage of nearby peaks that should be also labeled as peaks when one is found
#' @param training Vector of the data frames to be used for training
#' @param multiple_cores Number of cores to use
#' @export svm_rf

svm_rf <- function(list_of_dfs, labels, neighbors, min_peaks_percentage, error_window=0.005, training=1:length(list_of_dfs), multiple_cores=1) {
  require(parallel)
  cl <- makeCluster(mc <- getOption("cl.cores", multiple_cores))
  toReturn <- vector(mode="list", length=(length(neighbors) * length(min_peaks_percentage) * length(error_window)))
  cartprod <- expand.grid(neighbors, min_peaks_percentage, error_window)
  for (i in 1:length(toReturn)) {
          toReturn[[i]] <- list(lst=list_of_dfs, lbl=labels, trn=training, indices=cartprod[i,])
  }
  newReturn <- parLapply(cl, toReturn,
  #newReturn <- lapply(head(toReturn),
         function(x) {
           neighbor_index <- x$indices[,1]
           min_peaks_index <- x$indices[,2]
           error_index <- x$indices[,3]
           peaks <- lapply(x$lst, function(y) binary_peaks(y, neighbor_index, error_index))
           #Creates table where columns are mass spectras and rows are peak mz values
           cpeaks <- combine_peaks(peaks)
           #trainset <- cpeaks[training,]
           #trainlabels <- x$lbl[training]
           results <- classifier_accuracies(cpeaks, labels, training, min_peaks_index)
           list(RF=results[[1]], SVM=results[[2]], error_window=error_index, neighbors=neighbor_index, min_peaks_freq=min_peaks_index)
         }
         )
  stopCluster(cl)
  newReturn
}

#' Rank importance of features
#' 
#' Given a matrix of binary peaks and each row's corresponding labels, this function takes returns the absolute difference between the proportion of times an m/z value was labeled as a peak within each of the two classes.
#' @param peaks A matrix of peaks
#' @param labels A factor vector of labels whose length is equal to the number of rows of peaks
#' @export naive_feature_importance

naive_feature_importance <- function(peaks, labels) {
  importance <- apply(peaks, 2, function(x) {
                      tmp <- (tapply(x, labels, mean))
                      abs(tmp[1]-tmp[2])
  })
  names(importance) <- colnames(peaks)
  importance 
}

#' Plot importance of randomForest 
#' 
#' Plot importance of each m/z value according to a randomForest object.  Returns a ggplot object
#' @param rf A randomForest object
#' @export plot_rf_importance

plot_rf_importance<- function(rf, count) {
  rf.importance <- rf$importance[,1]
  importancedf <- data.frame(mz = as.numeric(names(rf.importance)), importance = rf.importance)
  list(ggplot(importancedf, aes(mz, importance)) + geom_point() + ggtitle("Random Forest: m/z vs importance"), head(sort(rf.importance, decreasing=T), n=count))
}

#' Plot importance of naive importance vector
#' 
#' Plot importance of each m/z value according to the naive ranking vector.  Returns a ggplot object
#' @param naive_importance A vector returned from naive_feature_importance()
#' @export plot_naive_importance

plot_naive_importance <- function(naive_importance, count=10) {
  importancedf <- data.frame(mz = as.numeric(names(naive_importance)), importance = naive_importance)
  list(ggplot(importancedf, aes(mz, importance)) + geom_point() + ggtitle("Naive Rank: m/z vs importance") + ylim(0, 1), head(sort(naive_importance, decreasing=T), n=count))
}

#' Get frequency of peaks
#' 
#' Return frequency of peaks within each label
#' @param peaks
#' @param labels
#' @export peak_frequencies

peak_frequencies <- function(peaks, labels) {
  apply(peaks, 2, function(x) tapply(x, labels, mean))
}

#' Plot important peaks on each day
#' 
#' Returns a ggplot object of the most important peak frequencies
#' @param peaks_freqs The frequency of each peak within each label
#' @param rf A random forest object, to be used for feature importance
#' @param count Number of peaks to select
#' @export plot_rf_differences

plot_rf_differences <- function(peak_freqs, rf, count) {
  mzs <- names(head(sort(rf$importance[,1], decreasing=T), count))
  filteredData <- peak_freqs[,colnames(peak_freqs) %in% mzs]
  tmpdf <- melt(filteredData)
  names(tmpdf) <- c("Label", "mz", "Peak_Frequency")
  tmpdf$Label <- as.factor(tmpdf$Label)
  list(ggplot(tmpdf, aes(mz, Peak_Frequency, color=Label)) + geom_point() + ggtitle("Most important m/z's vs Peak Frequency"), filteredData)
}

#' Generates matrix of accuracy values
#'
#' Returns a matrix of accuracy values where the columns represent the classifier type, and the rows represent the min_peak_count.neighbor values used for the classifier
#' @param rf A random forest object
#' @param svm A support vector machine object
#' @export accuracy_matrix

accuracy_matrix <- function(rf, svm) {
  RFaccuracy <- unlist(lapply(rf, function(x) sum(diag(x$confusion))/sum(x$confusion[,-3])))
  SVMaccuracy <- unlist(lapply(svm, function(x) mean(x$accuracies/100)))
  rbind(RFaccuracy, SVMaccuracy)
}

