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

  complete_nc_files <- c(basename(fs::dir_ls(path, glob = "index", invert = TRUE)),
                         basename(fs::dir_ls(path, glob = "preliminary.nc")))
  prelim_nc_files <- basename(fs::dir_ls(path, glob = "*_preliminary.nc"))

  # Complete exists
  complete_dates <- str_extract(complete_nc_files, paste0(subdir, "\\d{2}"))
  complete_dates_regexp <- paste(complete_dates, collapse = "|")

  ## Deletion
  deletion_candidates <- str_detect(prelim_nc_files, complete_dates_regexp)
  delete_these <- prelim_nc_files[deletion_candidates]
  if(!rlang::is_empty(delete_these)) fs::file_delete(paste0(path, "/", delete_these))

}

# Update the files
# February
#get_nc_files(subdir = "202402")
#clean_prelims(subdir = "202402")

#get_nc_files(subdir = "202403")
#clean_prelims(subdir = "202403")


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

layerinfo <- tibble(
  num = c(1:4),
  raw_name = c("anom_zlev=0", "err_zlev=0", "ice_zlev=0", "sst_zlev=0"),
  name = c("anom", "err", "ice", "sst"))


## The Terra way. Should be considerably faster. (And it is)
## Even faster with the chunked names, to maximize layered raster processing,
## Chunks of 25 elements or so seem to work quickly enough.
process_raster <- function(fnames, crop_area = c(-80, 0, 0, 60), layerinfo = layerinfo) {

  tdf <- terra::rast(fnames) |>
    terra::rotate() |>   # Fix 0-360 lon
    terra::crop(crop_area) # Manually crop to a defined box. Default is Atlantic lat/lon box

  wts <- terra::cellSize(tdf, unit = "km") # For scaling

  # global() calculates a quantity for the whole grid on a particular SpatRaster
  # so we get one weighted mean per file that comes in
  out <- data.frame(date = terra::time(tdf),
             means = terra::global(tdf, "mean", weights = wts, na.rm=TRUE))
  out$var <- rownames(out)
  out$var <- gsub("_.*", "", out$var)
  out <- reshape(out, idvar = "date", timevar = "var",
          direction = "wide")

  colnames(out) <- gsub("weighted_mean\\.", "", colnames(out))
  out
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

## Save out as a csv
write_csv(df, file = here("data", "natl_means.csv"))

### Alternative: do ALL the seas FAST with rasterize() and zonal()
## Seas of the world polygons
seas <- sf::read_sf(here("raw", "World_Seas_IHO_v3")) #|>
  # filter(NAME %in% c("South Atlantic Ocean",
  #                    "South Pacific Ocean",
  #                    "Indian Ocean",
  #                    "North Pacific Ocean",
  #                    "North Atlantic Ocean"))

seas_ids <- tibble(
  ID = c(1:5),
  sea = c("South Atlantic Ocean",
          "South Pacific Ocean",
          "Indian Ocean",
          "North Pacific Ocean",
          "North Atlantic Ocean")
)

## Rasterize the seas polygons using one of the nc files
## as a reference grid for the rasterization process
one_raster <- all_fnames[1]
seas_vect <- terra::vect(seas)
tmp_tdf_seas <- terra::rast(one_raster)["sst"] |>
  rotate()
seas_zonal <- rasterize(seas_vect, tmp_tdf_seas, "NAME")

# Take a look
plot(seas_zonal)

# Need to wrap the object (because C++ pointers)
# If we don't do this it can't be passed around
# across the processes that future_map() will spawn
#seas_vect <- wrap(terra::vect(seas))
seas_zonal_wrapped <- wrap(seas_zonal)


# This is WAY Quicker than the extract() method
process_raster_zonal <- function(fnames) {

  d <- terra::rast(fnames)
  wts <- terra::cellSize(d, unit = "km") # For scaling

  layer_varnames <- terra::varnames(d)
  date_seq <- rep(terra::time(d))
  # These are the new colnames zonal post-calculation below
  new_colnames <- c("sea", paste(layer_varnames, date_seq, sep = "_"))

  tdf_seas <- d |>
    terra::rotate() |>   # Fix 0-360 lon
    terra::zonal(unwrap(seas_zonal_wrapped), mean, na.rm = TRUE)

  colnames(tdf_seas) <- new_colnames

  tdf_seas |>
    tidyr::pivot_longer(-sea,
                        names_to = c("measure", "date"),
                        values_to = "value",
                        names_pattern ="(.*)_(.*)") |>
    tidyr::pivot_wider(names_from = measure, values_from = value)

}

## Parallelized, but even so, be patient.
tictoc::tic("big op")
seameans_df <- future_map(chunked_fnames, process_raster_zonal) |>
  list_rbind() |>
  mutate(date = ymd(date),
         year = lubridate::year(date),
         month = lubridate::month(date),
         day = lubridate::day(date),
         yrday = lubridate::yday(date),
         season = season(date))
tictoc::toc()

write_csv(seameans_df, file = here("data", "oceans_sst_means_zonal.csv"))

month_labs <- seameans_df |>
  filter(sea == "North Atlantic Ocean",
         year == 2023,
         day == 15) |>
  select(date, year, yrday, month, day) |>
  mutate(month_lab = month(date, label = TRUE, abbr = TRUE))

main_oceans <- c(
  "North Atlantic Ocean",
  "South Atlantic Ocean",
  "North Pacific Ocean",
  "South Pacific Ocean",
  "Indian Ocean")

out <- seameans_df |>
  filter(sea %in% main_oceans) |>
  mutate(year_flag = case_when(
    year == 2023 ~ "2023",
    year == 2024 ~ "2024",
    .default = "All other years"
  ),
  sea_f = factor(sea, levels = main_oceans, ordered = TRUE)) |>
  filter(sea != "Indian Ocean") |>
  ggplot(aes(x = yrday, y = sst_wt_mean, group = year, color = year_flag)) +
  geom_line(linewidth = rel(0.5)) +
  scale_x_continuous(breaks = month_labs$yrday, labels = month_labs$month_lab) +
  scale_color_manual(values = c("orange", "firebrick", "skyblue")) +
  guides(
    x = guide_axis(cap = "both"),
    y = guide_axis(minor.ticks = TRUE, cap = "both"),
    color = guide_legend(override.aes = list(linewidth = 1.4))
  ) +
  facet_wrap(~ sea_f, axes = "all_x", axis.labels = "all_y") +
  labs(x = "Month of the Year", y = "Mean Temperature (Celsius)",
       color = "Year",
       title = "Mean Daily Sea Surface Temperatures, 1981-2024",
       subtitle = "Area-weighted 0.25° grid estimates; NOAA OISST v2.1; IHO Sea Boundaries",
       caption = "Data processed with R; Figure made with ggplot by Kieran Healy / @kjhealy") +
  theme(axis.line = element_line(color = "gray30", linewidth = rel(1)),
        strip.text = element_text(face = "bold", size = rel(1.4)),
        plot.title = element_text(size = rel(1.525)),
        plot.subtitle = element_text(size = rel(1.1)))

ggsave(here("figures", "four_oceans.pdf"), out, width = 10, height = 10)

ggsave(here("figures", "four_oceans.png"), out, width = 10, height = 10, dpi = 300)



## All the world's oceans and seas
out <- seameans_df |>
  mutate(year_flag = case_when(
    year == 2023 ~ "2023",
    year == 2024 ~ "2024",
    .default = "All other years")) |>
  ggplot(aes(x = yrday, y = sst_wt_mean, group = year, color = year_flag)) +
  geom_line(linewidth = rel(0.5)) +
  scale_x_continuous(breaks = month_labs$yrday, labels = month_labs$month_lab) +
  scale_color_manual(values = c("orange", "firebrick", "skyblue")) +
  guides(
    x = guide_axis(cap = "both"),
    y = guide_axis(minor.ticks = TRUE, cap = "both"),
    color = guide_legend(override.aes = list(linewidth = 1.4))
  ) +
  facet_wrap(~ reorder(sea, sst_wt_mean), axes = "all_x", axis.labels = "all_y") +
  labs(x = "Month of the Year", y = "Mean Temperature (Celsius)",
       color = "Year",
       title = "Mean Daily Sea Surface Temperatures, 1981-2024",
       subtitle = "Area-weighted 0.25° grid estimates; NOAA OISST v2.1; IHO Sea Boundaries",
       caption = "Data processed with R; Figure made with ggplot by Kieran Healy / @kjhealy") +
  theme(axis.line = element_line(color = "gray30", linewidth = rel(1)),
        strip.text = element_text(face = "bold", size = rel(1.4)),
        plot.title = element_text(size = rel(1.525)),
        plot.subtitle = element_text(size = rel(1.1)))

ggsave(here("figures", "all_seas.pdf"), out, width = 40, height = 40)

ggsave(here("figures", "all_seas.png"), out, width = 40, height = 40, dpi = 300)



## North Atlantic only
dfp <- df


out_atlantic <- dfp |>
  mutate(year_flag = case_when(
    year == 2023 ~ "2023",
    year == 2024 ~ "2024",
    .default = "All other years"
  )) |>
  ggplot(aes(x = yrday, y = wt_mean_sst, group = year, color = year_flag)) +
  geom_line(linewidth = rel(1.1)) +
  scale_x_continuous(breaks = month_labs$yrday, labels = month_labs$month_lab) +
  scale_color_manual(values = c("orange", "firebrick", "lightblue")) +
  guides(
    x = guide_axis(cap = "both"),
    y = guide_axis(minor.ticks = TRUE, cap = "both"),
    color = guide_legend(override.aes = list(linewidth = 2))
  ) +
  labs(x = "Month", y = "Mean Temperature (Celsius)",
       color = "Year",
       title = "Mean Daily Sea Surface Temperature, North Atlantic Ocean, 1981-2024",
       subtitle = "Gridded and weighted NOAA OISST v2.1 estimates",
       caption = "Kieran Healy / @kjhealy") +
  theme(axis.line = element_line(color = "gray30", linewidth = rel(1)))

ggsave(here("figures", "north_atlantic.png"), out_atlantic, height = 7, width = 10, dpi = 300)


