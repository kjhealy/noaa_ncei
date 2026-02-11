# NOAA NCEI Global Sea Surface Temperature Anomaly

#install.packages("ncdf4")
#install.packages("CFtime")

library(tidyverse)
library(here)

# Set up parallel processing (6 background processes)
library(mirai)


# Summon daemons
daemons(15)
status()

## For the raster data
library(ncdf4)
library(CFtime)
library(terra)


## Functions
source(here("R", "functions.R"))

## Seasons for plotting
season <- function(in_date) {
  br = yday(as.Date(c("2019-03-01", "2019-06-01", "2019-09-01", "2019-12-01")))
  x = yday(in_date)
  x = cut(x, breaks = c(0, br, 366))
  levels(x) = c("Winter", "Spring", "Summer", "Autumn", "Winter")
  x
}


### Make the dfs

# For the initial get, use 01_setup_data_initial_get.R

# Update for a specific month
#get_nc_files(subdir = "202512")
#clean_prelims(subdir = "202512")

# Get filenames
# All the daily .nc files we downloaded:
all_fnames <- fs::dir_ls(here("raw"), recurse = TRUE, glob = "*.nc")
length(all_fnames)
chunked_fnames <- chunk(all_fnames, 25)


## North Atlantic
## Atlantic box
crop_bb <- c(-80, 0, 0, 60)

# Try one only
tmp <- process_raster(all_fnames[1])

chk <- process_raster(chunked_fnames[[200]])
chk |> as_tibble()

## wheeee
tictoc::tic("Terra Method")
df <- chunked_fnames |>
  map(in_parallel(
    \(x) process_raster(x),
    process_raster = process_raster
  )) |>
  list_rbind() |>
  as_tibble() |>
  mutate(
    date = ymd(date),
    year = lubridate::year(date),
    month = lubridate::month(date),
    day = lubridate::day(date),
    yrday = lubridate::yday(date),
    season = season(date)
  )
tictoc::toc()

## Save out as a csv
write_csv(df, file = here("data", "natl_means.csv"))


## All the seas

# All the seas with rasterize() and zonal()
# Seas of the world polygons from https://www.marineregions.org/downloads.php,
# IHO Sea Areas V3 shapefile.
seas <- sf::read_sf(here("raw", "World_Seas_IHO_v3"))

terra::rast(all_fnames[100]) |>
  terra::rotate()


## Rasterize the seas polygons using one of the nc files
## as a reference grid for the rasterization process
one_raster <- all_fnames[1]
seas_vect <- terra::vect(seas)
tmp_tdf_seas <- terra::rast(one_raster)["sst"] |>
  rotate()
seas_zonal <- rasterize(seas_vect, tmp_tdf_seas, "NAME")

# Take a look
#plot(seas_zonal)

# Need to wrap the object (because C++ pointers)
# If we don't do this it can't be passed around
# across the processes that in_parallel() will spawn
seas_zonal_wrapped <- terra::wrap(seas_zonal)


## Parallelized, but even so, be patient.
tictoc::tic("big op")
seameans_df <- chunked_fnames |>
  map(in_parallel(
    \(x) process_raster_zonal(x),
    process_raster_zonal = process_raster_zonal,
    seas_zonal_wrapped = seas_zonal_wrapped
  )) |>
  list_rbind() |>
  mutate(
    date = ymd(date),
    year = lubridate::year(date),
    month = lubridate::month(date),
    day = lubridate::day(date),
    yrday = lubridate::yday(date),
    season = season(date)
  )
tictoc::toc()

write_csv(seameans_df, file = here("data", "oceans_means_zonal.csv"))
saveRDS(seameans_df, here("data", "seameans_df.rda"), compress = "xz")


## World average

# World box (60S 60N)
world_crop_bb <- c(-180, 180, -60, 60)

tictoc::tic("Terra Method")
world_df <- chunked_fnames |>
  map(in_parallel(
    \(x) process_raster(x, crop_area = world_crop_bb),
    process_raster = process_raster,
    world_crop_bb = world_crop_bb
  )) |>
  list_rbind() |>
  as_tibble() |>
  mutate(
    date = ymd(date),
    year = lubridate::year(date),
    month = lubridate::month(date),
    day = lubridate::day(date),
    yrday = lubridate::yday(date),
    season = season(date)
  )
tictoc::toc()

write_csv(world_df, file = here("data", "global_means_60S60N.csv"))

## Banish daemons
daemons(0)
