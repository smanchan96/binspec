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

#' Reduce to one variable
#'
#' Reduces each persons m/z, frequencies, and days to just frequencies
#' @param ... Each parameter should be a person df
#' @export reduceAndConcat
#' @examples
#' reduceAndConcat()
reduceAndConcat <- function(...) {
  mmdb <- minMaxByDay(...)
  rowlen <- sum(sapply(1:nrow(mmdb), function(x) mmdb$max[x]-mmdb$min[x]+1))
  mymatrix <<- matrix(nrow=length(list(...)), ncol=rowlen)
  mydf <- data.frame("mz"=as.vector(unlist(sapply(1:nrow(mmdb), function(x) mmdb[x,1]:mmdb[x,2]))), "freq"=rep(0, rowlen), "ids"=rep("anon", rowlen), days=as.vector(unlist(sapply(1:nrow(mmdb), function(x) rep(x, mmdb[x,2]-mmdb[x,1]+1)))))
  for (i in 1:nrow(mymatrix)) {
    mymatrix[i,] <<- dfmerger(mydf, list(...)[[i]])
  }
  mymatrix
}

#' Reduce to one variable helper
#'
#' Takes in a data frame whose freq column hasn't been modified and fills it with matching rows from srcdf
#' @param filldf 
#' @param srcdf
#' @export dfmerger
#' @examples
#' dfmerger()
dfmerger <- function(filldf, srcdf) {
  frq <- rep(0, nrow(filldf))
  sapply(1:nrow(filldf), function(x) {
    tmp <- srcdf$freq[srcdf$mz==filldf$mz[x] & srcdf$day==filldf$day[x]]
    if (length(tmp)==1){
      frq[x]<<-tmp
    }
})
  frq
}

#' Returns the min and max of each day
#'
#' Reduces each persons m/z, frequencies, and days to just frequencies
#' @param ... Each arameter should be a person df
#' @export minMaxByDay
#' @examples
#' minMaxByDay()
minMaxByDay <- function(...) {
  myDays <- unique(unlist(lapply(list(...), function(x) unique(x$days))))
  mins <- sapply(myDays, function(x)min(unlist(lapply(list(...), function(y)min(y$freq[y$days==x])))))
  maxs <- sapply(myDays, function(x)max(unlist(lapply(list(...), function(y)max(y$freq[y$days==x])))))
  data.frame("min"=mins, "max"=maxs)
}
