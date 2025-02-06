#' Get NCDF Files from NOAA
#'
#' @param url The endpoint URL of the AVHRR data, <https://www.ncei.noaa.gov/data/sea-surface-temperature-optimum-interpolation/v2.1/access/avhrr/>
#' @param local A local file path for the raw data folders, i.e. where all the year-month dirs go. Defaults to a local version under raw/ of the same path as the NOAA website.
#' @param subdir The subdirectory of monthly data to get. A character string of digits, of the form "YYYYMM". No default.
#'
#' @return  A directory of NCDF files.
#' @export
#'
#'
get_nc_files <- function(url = "https://www.ncei.noaa.gov/data/sea-surface-temperature-optimum-interpolation/v2.1/access/avhrr/",
                         local = here::here("raw/www.ncei.noaa.gov/data/sea-surface-temperature-optimum-interpolation/v2.1/access/avhrr/"),
                         subdir) {

  localdir <- here::here(local, subdir)

  if(!fs::dir_exists(localdir)) fs::dir_create(localdir)

  files <- rvest::read_html(paste0(url,subdir)) |>
    rvest::html_elements("a") |>
    rvest::html_text2()
  files <- subset(files, str_detect(files, "nc"))

  full_urls <- paste0(url, subdir, "/", files)
  full_outpaths <- paste0(localdir, "/", files)

  walk2(full_urls, full_outpaths, \(x,y) httr::GET(x, httr::write_disk(y, overwrite = TRUE)))

}

#
#' Remove -prelim nc files if final nc version exists
#'
#' @param subdir A character vector. The name of the subdirectory of .nc files you want to clean up.
#' @param local A file path. The name of the subdirectory where all the year-month dirs are. Defaults to a local version under raw/ of the same path as the NOAA website.
#'
#' @return Returns "No duplicates" if there are none; otherwise silently deletes dupes.
#' @export
#'
#' @examples
clean_prelims <- function(subdir,
                          local = here::here("raw/www.ncei.noaa.gov/data/sea-surface-temperature-optimum-interpolation/v2.1/access/avhrr/")) {

  path <- here::here(local, subdir)

  all_nc_files <- basename(fs::dir_ls(path, regexp = "[.]html", invert = TRUE))
  prelim_nc_files <- basename(fs::dir_ls(path, glob = "*_preliminary.nc"))
  final_nc_files <- all_nc_files[str_detect(all_nc_files, "_preliminary",
                                            negate = TRUE)]

  prelim_nc_dates <- str_extract(prelim_nc_files,
                                 paste0(subdir, "\\d{2}"))
  final_nc_dates <- str_extract(final_nc_files,
                                paste0(subdir, "\\d{2}"))

  dupes <- intersect(prelim_nc_dates, final_nc_dates)

  if(rlang::is_empty(dupes)) return(message("No duplicates"))


  ## Deletion
  dupes <- paste0(dupes, collapse = "|")
  deletion_candidates <- str_detect(prelim_nc_files, dupes)
  delete_these <- prelim_nc_files[deletion_candidates]
  if(!rlang::is_empty(delete_these)) fs::file_delete(paste0(path, "/", delete_these))

}


## For the filename processing
## This one gives you an unknown number of chunks each with approx n elements
chunk <- function(x, n) split(x, ceiling(seq_along(x)/n))

## This one gives you n chunks each with an approx equal but unknown number of elements
chunk2 <- function(x, n) split(x, cut(seq_along(x), n, labels = FALSE))

## The Terra way. Should be considerably faster. (And it is)
## Even faster with the chunked names, to maximize layered raster processing,
## Chunks of 25 elements or so seem to work quickly enough.
layerinfo <- tibble(
  num = c(1:4),
  raw_name = c("anom_zlev=0", "err_zlev=0",
               "ice_zlev=0", "sst_zlev=0"),
  name = c("anom", "err",
           "ice", "sst"))

process_raster <- function(fnames, crop_area = c(-80, 0, 0, 60), layerinfo = layerinfo) {

  tdf <- terra::rast(fnames) |>
    terra::rotate() |>   # Convert 0 to 360 lon to -180 to +180 lon
    terra::crop(crop_area) # Manually crop to a defined box.  Default is roughly N. Atlantic lat/lon box

  wts <- terra::cellSize(tdf, unit = "km") # For scaling

  # global() calculates a quantity for the whole grid on a particular SpatRaster
  # so we get one weighted mean per file that comes in
  out <- data.frame(date = terra::time(tdf),
                    means = terra::global(tdf, "mean", weights = wts, na.rm=TRUE))
  out$var <- rownames(out)
  out$var <- gsub("_.*", "", out$var)
  out <- reshape(out, idvar = "date",
                 timevar = "var",
                 direction = "wide")

  colnames(out) <- gsub("weighted_mean\\.", "", colnames(out))
  out
}

# This is much faster than using extract()
process_raster_zonal <- function(fnames) {

  d <- terra::rast(fnames)
  wts <- terra::cellSize(d, unit = "km") # For scaling

  layer_varnames <- terra::varnames(d) # vector of layers
  date_seq <- rep(terra::time(d)) # vector of dates

  # New colnames for use post zonal calculation below
  new_colnames <- c("sea", paste(layer_varnames, date_seq, sep = "_"))

  # Better colnames
  tdf_seas <- d |>
    terra::rotate() |>   # Convert 0 to 360 lon to -180 to +180 lon
    terra::zonal(unwrap(seas_zonal_wrapped), mean, na.rm = TRUE)
  colnames(tdf_seas) <- new_colnames

  # Reshape to long
  tdf_seas |>
    tidyr::pivot_longer(-sea,
                        names_to = c("measure", "date"),
                        values_to = "value",
                        names_pattern ="(.*)_(.*)") |>
    tidyr::pivot_wider(names_from = measure, values_from = value)

}
