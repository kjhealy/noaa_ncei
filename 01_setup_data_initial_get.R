### NOAA NCEI Global Sea Surface Temperature Anomaly

library(tidyverse)
library(here)

## Functions
source(here("R", "functions.R"))

### Initial get.
## This will take quite some time, may time out, and we try to be nice.

# Data collection starts in September 1981
first_yr <- paste0("1981", sprintf('%0.2d', 9:12))
yrs <- 1982:2024
months <- sprintf('%0.2d', 1:12)
subdirs <- c(first_yr, paste0(rep(yrs, each = 12), months))

politely_get_nc_files <- slowly(get_nc_files)

walk(subdirs,
     \(x) politely_get_nc_files(subdir = x))
