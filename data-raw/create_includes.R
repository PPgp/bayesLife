library(data.table)
library(usethis)

# Run this from the data-raw directory

include_2022 <- as.data.frame(fread("include_2022.txt", sep = "\t")[, .(country_code, include_code)])
use_data(include_2022, overwrite = TRUE)

include_2019 <- as.data.frame(fread("include_2019.txt", sep = "\t")[, .(country_code, include_code)])
use_data(include_2019, overwrite = TRUE)

include_2017 <- as.data.frame(fread("include_2017.txt", sep = "\t")[, .(country_code, include_code)])
use_data(include_2017, overwrite = TRUE)

include_2015 <- as.data.frame(fread("include_2015.txt", sep = "\t")[, .(country_code, include_code)])
use_data(include_2015, overwrite = TRUE)

include_2012 <- as.data.frame(fread("include_2012.txt", sep = "\t")[, .(country_code, include_code)])
use_data(include_2012, overwrite = TRUE)

include_2010 <- as.data.frame(fread("include_2010.txt", sep = "\t")[, .(country_code, include_code)])
use_data(include_2010, overwrite = TRUE)


