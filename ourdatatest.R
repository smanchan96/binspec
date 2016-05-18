#load packages
library(binspec)
#Test our functions
filenamesRoom <- list.files("/scratch/mentors/nkong/smanchan/Newest_Data/room_temp_106/", full.names = T)
roomtemps <- lapply(filenamesRoom, function(x) {
                    y <- round_df(read.csv(x))
                    names(y) <- c("mz", "freq")
                    y
}
)

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
labels <- roomdays[roomdays == "1" | roomdays == "0"]

roompeaks <- lapply(roomtemps, function(x) binary_peaks(x, 10, .005))
#Start test
peaks <- combine_peaks(roompeaks)
peaks
labels
