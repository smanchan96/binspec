#load packages
library(binspec)
library(e1071)
library(randomForest)
#Load the data
filenamesRoom <- list.files("/scratch/mentors/nkong/smanchan/Newest_Data/room_temp_106/", full.names = T)
roomtemps <- lapply(filenamesRoom, function(x) {
  y <- round_df(read.csv(x))
  names(y) <- c("mz", "freq")
  y
})

extractname <- function(str) {
  hyphsplit <- strsplit(str, "-")[[1]]
  id <- strsplit(hyphsplit[2], "\\.")[[1]][1]
  mystr <- tail(strsplit(hyphsplit[1], "/")[[1]], n=1)
  m <- gregexpr('[0-9]+', mystr)
  day <- strtoi(regmatches(mystr, m))
  val <- c(id, day)
  names(val) <- c("id", "day")
  val
}
roomtempnames <- lapply(filenamesRoom, extractname)
roomdays <- unlist(lapply(roomtempnames, function(x) x[2]))
roomtemps <- roomtemps[roomdays == "1" | roomdays == "0"]
room_labels <- as.factor(roomdays[roomdays == "1" | roomdays == "0"])

filenames_four <- list.files("/scratch/mentors/nkong/smanchan/Newest_Data/4_degree_117/days0_1_2_5_10/", full.names = T)
fourtemps <- lapply(filenames_four, function(x) {
  y <- round_df(read.csv(x))
  names(y) <- c("mz", "freq")
  y
})
fourtempnames <- lapply(filenames_four, extractname)
fourdays <- unlist(lapply(fourtempnames, function(x) x[2]))
fourtemps <- fourtemps[fourdays == "1" | fourdays == "0"]
four_labels <- as.factor(fourdays[fourdays == "1" | fourdays == "0"])

###Use of library functions starts here
min_peaks_count <- c(30, 50, 70, 90, 110)
neighbors <- c(5, 10, 15, 20)
four_rf <- matrix(data=NA, nrow = length(neighbors), ncol=length(min_peaks_count))
four_svm <- matrix(data=NA, nrow = length(neighbors), ncol=length(min_peaks_count))
room_rf <- matrix(data=NA, nrow = length(neighbors), ncol=length(min_peaks_count))
room_svm <- matrix(data=NA, nrow = length(neighbors), ncol=length(min_peaks_count))
rownames(four_rf) <- as.character(neighbors)
rownames(four_svm) <- rownames(four_rf)
rownames(room_rf) <- rownames(four_rf)
rownames(room_svm) <- rownames(four_rf)
colnames(four_rf) <- as.character(min_peaks_count)
colnames(four_svm) <-colnames(four_rf)
colnames(room_rf) <- colnames(four_rf)
colnames(room_svm) <- colnames(four_rf)

for (i in 1:length(neighbors)) {
  roompeaks <- lapply(roomtemps, function(x) binary_peaks(x, neighbors[i], .005))
  fourpeaks <- lapply(fourtemps, function(x) binary_peaks(x, neighbors[i], .005))
  
  room_peaks <- combine_peaks(roompeaks)
  four_peaks <- combine_peaks(fourpeaks)
  
  for (j in 1:length(min_peaks_count)) {
    room_results <- classifier_accuracies(room_peaks, room_labels, min_peaks_count[j])
    four_results <- classifier_accuracies(four_peaks, four_labels, min_peaks_count[j])
    four_rf[i,j] <- four_results[1]
    four_svm[i,j] <- four_results[2]
    room_rf[i,j] <- room_results[1]
    room_svm[i,j] <- room_results[2]
  }
}
write.table(four_svm, "four_svm.txt")
write.table(four_rf, "four_rf.txt")
write.table(room_svm, "room_svm.txt")
write.table(room_rf, "room_rf.txt")
