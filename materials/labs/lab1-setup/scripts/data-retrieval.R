library(tidyverse)

# retrieve pollution data
url <- 'https://raw.githubusercontent.com/pstat197/pstat197a/main/materials/labs/lab1-setup/data/pollution.csv'
pollution <- read_csv(url)

# write as csv to file
write_csv(pollution, file = 'data/pollution.csv')

# clear environment
rm(list = ls())