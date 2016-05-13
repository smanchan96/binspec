setwd("/scratch/mentors/nkong/smanchan/ibp100txt/Separated/BloodDraw_Day1/")
filenames0 <- list.files("./0d/", full.names = T)
filenames1 <- list.files("./2d/", full.names = T)
filenames2 <- list.files("./3d/", full.names = T)
filenames3 <- list.files("./6d/", full.names = T)
filenames4 <- list.files("./11d/", full.names = T)

day0 <- lapply(filenames0, read.csv)
day1 <- lapply(filenames1, read.csv)
day2 <- lapply(filenames2, read.csv)
day3 <- lapply(filenames3, read.csv)
day4 <- lapply(filenames4, read.csv)
rm(filenames1, filenames2, filenames3, filenames4)
setwd("/scratch/mentors/nkong/smanchan/Code/clustspec/")

library(clustspec)
ppl <- lapply(1:length(day0), function(x) list(day0[[x]], day1[[x]], day2[[x]], day3[[x]], day4[[x]]))
persons <- lapply(1:length(ppl), function(x) do.call(function(...) persondf(toString(x), c(0,2,3,6,11), ...), ppl[[x]]))

adjustdf <- function(df, min, max) {
  newdf <- data.frame(min:max, min:max)
  j <- 1
  for(i in 1:length(newdf[[1]])) {
    if (j <= length(df[[1]]) && newdf[[1]][i] == round(df[[1]][j])) {
      newdf[[2]][i] = df[[2]][j]
      j = j + 1
    }
    else {
      newdf[[2]][i] = 0
    }
  }
  return(newdf[[2]])
}
adjustlist <- function(daylists) {
  max <- round(max(unlist(lapply(daylists, function(x)max(x[[1]])))))
  min <- round(min(unlist(lapply(daylists, function(x)min(x[[1]])))))
  # Populate dfs with either the actual value or 0 if not present
  lapply(daylists, function(x)adjustdf(x, min, max))
}

day0<-adjustlist(day0)
day1<-adjustlist(day1)
day2<-adjustlist(day2)
day3<-adjustlist(day3)
day4<-adjustlist(day4)

vectors <- do.call(cbind, list(do.call(rbind,day0),do.call(rbind,day1),do.call(rbind,day2),do.call(rbind,day3),do.call(rbind,day4)))

################################################################
#                                                              #
# Use of package to find results begins here.  n=53, d=15803   #
#                                                              #
################################################################

pcs <- getPrComps(vectors)
km <- kmeansCH(pcs)
fs <- forestSelect(vectors, km$cluster, 10)

getmz <- function(index) {
  freqs <- vectors[,index]
  eachmz <- lapply(1:length(freqs), function(x) if (freqs[x]!= 0) {
    persons[[x]][persons[[x]]$freq==freqs[x],]
  } else return(NULL))
  return(eachmz)
  # newlst <- matrix(ncol=4)
  # colnames(newlst) <- c("mz", "freq", "ids", "days")
  # for (i in eachmz) {
  #   if (!is.null(i) & nrow(i)!=0) {
  #     newlst<-rbind(newlst, i)
  #   }
  # }
  # c(names(sort(table(round(newlst$mz)), decreasing = T)[1]), names(sort(table(newlst$days))[1]))
}

cl <- sapply(c(5, 10, 20, 50, 100, 200, 500), function(r) {
  fs <- forestSelect(vectors, km$cluster, r)
  indices <- sapply(names(fs), strtoi)
  newvec <- vectors[,indices]
  sort(table(kmeansCH(newvec)$cluster), decreasing = T)
})
res <- sapply(c(5, 10, 20, 50, 100, 200, 500), function(r) {
  fs <- forestSelect(vectors, km$cluster, r)
  indices <- sapply(names(fs), strtoi)
  newvec <- vectors[,indices]
  newkm <- kmeansCH(newvec)
  randIndex(km$cluster, newkm$cluster)
})
colnames(res) <- c(5,10,20,50,100,200,500)
colnames(cl) <- c(5,10,20,50,100,200,500)

################################################################
#                                                              #
# Visualization of results begins here                         #
#                                                              #
################################################################
plot(strtoi(colnames(res)), res[1,], xlab="r", ylab="Rand Index", main="Features vs Rand Index")


if (is.list(cl)) {
  n <<- max(unlist(lapply(cl, length)))
  for(i in 1:length(cl)) {
    length(cl[[i]]) <- n
  }
  newcl <<- do.call(cbind, cl)
} else {
  newcl <<- cl
  n <<- ncol(newcl)
}
m <- sort(table(km$cluster), decreasing = T)
length(m) <- n
barplot(t(cbind(m, newcl)), beside = T, xlab="Cluster", ylab="Member count", main="Cluster Size", names.arg = 1:nrow(cbind(m, newcl)), legend.text = c("C", "5", "10", "20", "50", "100", "200", "500"))

