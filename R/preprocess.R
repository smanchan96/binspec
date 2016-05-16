#' Round data frame
#'
#' Find the intensity peaks
#' @param df Data frame
#' @export round_df
#' @examples
#' integer_round_df()

round_df <- function(df) {
    #Round each column to the nearest integer
    data.frame(apply(df, 2, function(x) sapply(x, round)))
}

#' Find binary peaks
#'
#' Find binary peaks
#' @param df Data frame of m/z and intensities
#' @param neighbors Number of neighboring m/z values to compare on right and left
#' @param error m/z measuring error parameter
#' @export binary_peaks
#' @examples
#' binary_peaks()

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
    unlist(lapply(peakmzs, function(x) {
                       floor(x - error * x):ceiling(x + error * x) 

    }))
}

