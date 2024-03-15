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

## Functions
get_nc_files <- function(url = "https://www.ncei.noaa.gov/data/sea-surface-temperature-optimum-interpolation/v2.1/access/avhrr/",
                       subdir) {
  local <- here::here("raw/www.ncei.noaa.gov/data/sea-surface-temperature-optimum-interpolation/v2.1/access/avhrr/")

  localdir <- paste0(local, subdir)

  if(!fs::dir_exists(localdir)) fs::dir_create(localdir)

  files <- rvest::read_html(paste0(url,subdir)) |>
    rvest::html_elements("a") |>
    rvest::html_text2()
  files <- subset(files, str_detect(files, "nc"))

  full_urls <- paste0(url, subdir, "/", files)
  full_outpaths <- paste0(localdir, "/", files)

  walk2(full_urls, full_outpaths, \(x,y) httr::GET(x, httr::write_disk(y, overwrite = TRUE)))

}

# Remove prelim files if finalized version exists
clean_prelims <- function(subdir) {
  local <- here::here("raw/www.ncei.noaa.gov/data/sea-surface-temperature-optimum-interpolation/v2.1/access/avhrr/")
  path <- paste0(local, subdir)

  complete_nc_files <- basename(fs::dir_ls(path, regexp = "index|._preliminary.nc", invert = TRUE))
  prelim_nc_files <- basename(fs::dir_ls(path, glob = "*_preliminary.nc"))

  # Complete exists
  complete_dates <- str_extract(complete_nc_files, paste0(subdir, "\\d{2}"))
  complete_dates_regexp <- paste(complete_dates, collapse = "|")

  ## Deletion
  deletion_candidates <- str_detect(prelim_nc_files, complete_dates_regexp)
  delete_these <- prelim_nc_files[deletion_candidates]
  if(!rlang::is_empty(delete_these)) fs::file_delete(paste0(path, "/", prelim_nc_files[delete_these]))

}

# Update the files
# February
# get_nc_files(subdir = "202402")
# clean_prelims(subdir = "202402")

## Seasons for plotting
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

## For the filename processing
## This one gives you an unknown number of chunks each with approx n elements
chunk <- function(x, n) split(x, ceiling(seq_along(x)/n))

## This one gives you n chunks each with an approx equal but unknown number of elements
chunk2 <- function(x, n) split(x, cut(seq_along(x), n, labels = FALSE))




## For the raster data
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

## The Terra way. Should be considerably faster. (And it is)
## Even faster with the chunked names, to maximize layered raster processing,
## Chunks of 25 elements or so seem to work quickly enough.

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

## Filenames
## All the daily .nc files we downloaded:
all_fnames <- fs::dir_ls(here("raw"), recurse = TRUE, glob = "*.nc")
chunked_fnames <- chunk(all_fnames, 25)


## Atlantic box
crop_bb <- c(-80, 0, 0, 60)

# Try one only
chk <- process_raster(chunked_fnames[[500]])
chk

## wheeee
## Process >15,000 files
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

## Save out as a csv
write_csv(df, file = here("data", "sst_means.csv"))

## The plot
dfp <- df


dfp |>
  mutate(year_flag = case_when(
    year == 2023 ~ "2023",
    year == 2024 ~ "2024",
    .default = "All other years"
  )) |>
  ggplot(aes(x = yrday, y = wt_mean_sst, group = year, color = year_flag)) +
  geom_line(linewidth = rel(0.8)) +
  scale_x_continuous(breaks = season_lab$yrday, labels = season_lab$lab) +
  scale_color_manual(values = c("orange", "firebrick", "gray70")) +
  guides(
    x = guide_axis(minor.ticks = TRUE, cap = "both"),
    y = guide_axis(minor.ticks = TRUE, cap = "both"),
    color = guide_legend(override.aes = list(linewidth = 1.4))
  ) +
  labs(x = "Season", y = "Mean Temperature (Celsius)",
       color = "Year",
       title = "Mean Daily Sea Surface Temperature, North Atlantic Ocean, 1981-2024",
       subtitle = "Gridded and weighted NOAA OISST v2.1 estimates",
       caption = "Data processed with R; Figure made with ggplot") +
  theme(axis.line = element_line(color = "gray30", linewidth = rel(1)))

# ## Checks
# atl_array <- tmp_array_sst[lon > lon_2, lat > lat_1 & lat < lat_2]
# atl_grid <- expand.grid(lon = lon[lon > lon_2],
# lattice::levelplot(atl_array ~ lon * lat, data = atl_grid, pretty = T)
# grid <- expand.grid(lon = lon,
#                     lat = lat)
# lattice::levelplot(tmp_array_sst ~ lon * lat, data = grid, pretty = T)
