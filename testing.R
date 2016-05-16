#load packages
library(binspec)
##Test our functions
#set.seed(42)
df <- data.frame(x=1:100, y=runif(100))
df
round_df(df)
binary_peaks(df, 1, 0)
