#' Person constructor
#'
#' Creates a data frame with ids, m/z values, frequencies, and days. 
#' @param id What is this person's name
#' @param days Vector of days since initial sample taken
#' @param ... All the data frames of m/z values and frequencies
#' @export persondf
#' @examples
#' persondf()

persondf <- function(id, dates, ...) {
  mydf <- do.call(rbind, list(...))
  lengths <- unlist(lapply(list(...), nrow))
  names(mydf) = c("mz", "freq")
  mydf$ids <- rep(id, length(mydf$mz))
  mydf$days <- unlist(lapply(1:length(lengths), function(i) rep(dates[i], lengths[i])))
  return(mydf)
}
