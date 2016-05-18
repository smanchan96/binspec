#load packages
library(binspec)
##Test our functions
set.seed(42)
df <- data.frame(x=1:100, y=sample(1:100, 100, replace=T))
df <- round_df(df)
df2 <- data.frame(x=1:200, y=sample(1:100, 200, replace=T))
df2 <- round_df(df2)
peaklist <- list(binary_peaks(df, 10, .05), binary_peaks(df2, 10, .05))
peaklist
combine_peaks(peaklist)
