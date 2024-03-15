# NOAA NCEI Global Sea Surface Temperature Anomaly

#install.packages("ncdf4")
#install.packages("CFtime")


library(tidyverse)
library(here)
library(furrr)
plan(multisession)

## Plot setup
library(systemfonts)
clear_registry()

register_variant(
  name = "Myriad Pro SemiCondensed",
  family = "Myriad Pro",
  width = "semicondensed",
  weight = c("normal", "semibold"),
)

library(showtext)
showtext_opts(dpi = 300)
showtext_auto()

library(myriad)
import_myriad_semi()
import_myriad_condensed()

theme_set(theme_myriad_semi())


season <-  function(in_date){
  br = yday(as.Date(c("2019-03-01",
                      "2019-06-01",
                      "2019-09-01",
                      "2019-12-01")))
  x = yday(in_date)
  x = cut(x, breaks = c(0, br, 366))
  levels(x) = c("Winter", "Spring", "Summer", "Autumn", "Winter")
  x
}


season_lab <-  tibble(yrday = yday(as.Date(c("2019-03-01",
                                             "2019-06-01",
                                             "2019-09-01",
                                             "2019-12-01"))),
                      lab = c("Spring", "Summer", "Autumn", "Winter"))


## Data

library(ncdf4)
library(CFtime)
library(terra)



## Single slice. Defaults to North Atlantic rectangle
## I'll keep this fn here, but I won't use it: the Terra way is faster
get_anom <- function(fname, lat_1=0, lat_2=60, lon_1=0, lon_2=280) {
  dname_sst <- "sst"
  dname_anom <- "anom"

  ncin <- nc_open(fname)

  lon <- ncvar_get(ncin,"lon")
  lat <- ncvar_get(ncin,"lat")
  time <- ncvar_get(ncin,"time")
  tunits <- ncatt_get(ncin,"time","units")

  tmp_array_sst <- ncvar_get(ncin,dname_sst)
  fillvalue_dl_sst <- ncatt_get(ncin,dname_sst,"_FillValue")

  tmp_array_anom <- ncvar_get(ncin,dname_anom)
  fillvalue_dl_anom <- ncatt_get(ncin,dname_anom,"_FillValue")

  timestamp <- CFtimestamp(CFtime(tunits$value, calendar = "proleptic_gregorian", time))

  tmp_array_sst[tmp_array_sst==fillvalue_dl_sst$value] <- NA
  tmp_array_anom[tmp_array_anom==fillvalue_dl_anom$value] <- NA

  lon_sset <- lon[lon > lon_2]
  lat_sset <- lat[lat > lat_1 & lat < lat_2]
  lonlat_sset <- expand.grid(lon = lon_sset, lat = lat_sset)

  tmp_vec_anom <- as.vector(tmp_array_anom[lon > lon_2, lat > lat_1 & lat < lat_2])
  tmp_vec_sst <- as.vector(tmp_array_sst[lon > lon_2, lat > lat_1 & lat < lat_2])

  tibble(
    lon = lonlat_sset$lon,
    lat = lonlat_sset$lat,
    anom = tmp_vec_anom,
    sst = tmp_vec_sst,
    wt = cos(lat * (pi/180))) |>
    summarize(date = lubridate::as_datetime(timestamp),
              mean_anom = mean(anom, na.rm = TRUE),
              mean_sst = mean(sst, na.rm = TRUE),
              wt_mean_anom = weighted.mean(anom, wt, na.rm = TRUE),
              wt_mean_sst = weighted.mean(sst, wt, na.rm = TRUE))

}

## All the daily .nc files we downloaded:
all_fnames <- fs::dir_ls(here("raw"), recurse = TRUE, glob = "*.nc")
chunked_fnames <- split(all_fnames, ceiling(seq_along(all_fnames)/25))

## The Terra way. Should be considerably faster. (And it is)
## Even faster with the chunked names, to maximize layered raster processing
# Atlantic
crop_bb <- c(-80, 0, 0, 60)

process_raster <- function(fnames, crop_area = crop_bb, var = "sst") {
  tdf <- terra::rast(fnames)[var] |>
    terra::rotate() |>   # Fix 0-360 lon
    terra::crop(crop_bb) # Atlantic region

  wts <- terra::cellSize(tdf, unit = "km") # For scaling

  data.frame(
    date = terra::time(tdf),
    # global() calculates a quantity for the whole grid on a particular SpatRaster
    # so we get one weighted mean per file that comes in
    wt_mean_sst = global(tdf, "mean", weights = wts, na.rm=TRUE)$weighted_mean
  )

}

# Try one only
# fname <- all_fnames[5000]
# chk <- process_raster(fname)

chk <- process_raster(chunked_fnames[[1]])


## wheeee
## Process >15,000 files, but in chunks!
tictoc::tic("Terra Method")
df <- future_map(chunked_fnames, process_raster) |>
  list_rbind() |>
  as_tibble() |>
  mutate(date = ymd(date),
         year = lubridate::year(date),
         month = lubridate::month(date),
         day = lubridate::day(date),
         yrday = lubridate::yday(date),
         season = season(date))
tictoc::toc()


# tictoc::tic("NCDF4 Method")
# ## wheeee
# ## Process >15,000 files
# df_ncdf <- future_map(all_fnames, get_anom) |>
#   list_rbind()
#
# df_ncdf
# tictoc::toc()

## The plot

dfp <- df


dfp |>
  mutate(year_flag = case_when(
    year == 2023 ~ "2023",
    year == 2024 ~ "2024",
    .default = "All other years"
  )) |>
  ggplot(aes(x = yrday, y = wt_mean_sst, group = year, color = year_flag)) +
  geom_line() +
  scale_x_continuous(breaks = season_lab$yrday, labels = season_lab$lab) +
  scale_color_manual(values = c("orange", "firebrick", "gray70")) +
  labs(x = "Season", y = "Mean Temperature (Celsius)",
       color = "Year",
       title = "Mean Daily Sea Surface Temperature, North Atlantic Ocean, 1981-2024",
       subtitle = "Gridded and weighted NOAA OISST v2.1 estimates",
       caption = "Data processed with R; Figure made with ggplot")





## Checks
atl_array <- tmp_array_sst[lon > lon_2, lat > lat_1 & lat < lat_2]

atl_grid <- expand.grid(lon = lon[lon > lon_2],
                        lat = lat[lat > lat_1 & lat < lat_2])

lattice::levelplot(atl_array ~ lon * lat, data = atl_grid, pretty = T)

grid <- expand.grid(lon = lon,
                    lat = lat)
lattice::levelplot(tmp_array_sst ~ lon * lat, data = grid, pretty = T)
