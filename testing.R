#load packages
library(binspec)
##Test our functions
#set.seed(42)
df <- data.frame(1:100, runif(100))
round_df(df)

