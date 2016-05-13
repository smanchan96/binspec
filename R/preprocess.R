#' Find binary peaks
#'
#' Find the intensity peaks
#' @param df Data frame of m/z and intensities
#' @param neighbors Number of neighboring m/z values to compare on right and left
#' @param smooth Smoothing parameter
#' @export binary_peaks
#' @examples
#' binary_peaks()

binary_peaks <- function(df, neighbors, smooth) {
  sapply(1:nrow(df), function(x) {
    window <- intersect((x - neighbors):(x + neighbors), 1:nrow(df))
    all(df[[2]][x] >= df[[2]][window])
  })
}