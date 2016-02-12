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
  p1  <- list(...)
  for (x in 1:5) {
    names(p1[[x]])<-c("m/z", "freq")
  }
  mydf <- do.call(rbind, p1)
  lengths <- unlist(lapply(p1, nrow))
  names(mydf) = c("mz", "freq")
  mydf$ids <- rep(id, length(mydf$mz))
  mydf$days <- unlist(lapply(1:length(lengths), function(i) rep(dates[i], lengths[i])))
  return(mydf)
}

#' Plot person
#'
#' Plots the person by day
#' @param persondf Person data frame
#' @export persondf.plot
#' @examples
#' persondf.plot()
persondf.plot <- function(persondf){
  require(lattice)
  xyplot(persondf$freq~persondf$mz|sapply(persondf$days, function(x)paste("Day ", x)), xlab="m/z", ylab="Frequency", main=persondf$id[1])
}

# Use fread to read the csv files