# Load required libraries
setwd("")

library(raster)
library(tidyverse, lib = "Packages/")
library(trend, lib = "Packages/")
library(ncdf4)
library(doParallel, lib = "Packages/")
library(foreach, lib = "Packages/")
library(future, lib = "Packages/")
library(future.apply, lib = "Packages/")
library(modifiedmk, lib = "Packages/")

Create_New <- function(groups) {
  # Set working directory and get list of .nc files
  setwd("")
  nc_files <- list.files(pattern = ".nc")
  
  # Set up parallel processing
  num_cores <- as.integer(Sys.getenv('SLURM_CPUS_PER_TASK'))
  plan(multicore, workers = num_cores)
  
  # Month mapping and categorization
  months <- list(
    Jan = "001",
    Feb = "032",
    March = c("060", "061"),
    April = c("091", "092"),
    May = c("121", "122"),
    June = c("152", "153"),
    July = c("182", "183"),
    Aug = c("213", "214"),
    Sep = c("244", "245"),
    Oct = c("274", "275"),
    Nov = c("305", "306"),
    Dec = c("335", "336")
  )
  
  file_names <- substr(nc_files, 6, 8)
  month_mapping <- unlist(months)
  names(month_mapping) <- rep(names(months), sapply(months, length))
  file_month_mapping <- sapply(as.character(file_names), function(number) {
    month_name <- names(month_mapping)[month_mapping == number]
    return(month_name)
  })
  
  # Define seasons
  Summer <- which(file_month_mapping %in% c('Dec', 'Jan', 'Feb'))
  
  # Stack raster files for Summer season and crop to new extent
  setwd("")
  raster_stack <- stack(nc_files[Summer], varname = groups)
  new_extent <- extent(raster_stack)
  new_extent@ymin <- -90
  new_extent@ymax <- -30
  cropped_stack <- crop(raster_stack, new_extent)
  
  # Prepare data for Sens slope calculation
  xy <- data.frame(xyFromCell(cropped_stack[[1]], 1:ncell(cropped_stack[[1]])))
  lon <- as.vector(xy[, 1])
  lat <- as.vector(xy[, 2])
  
  # Transpose cropped_stack data for further processing
  tvalues_df <- as.data.frame(t(as.data.frame(cropped_stack)))
  
  # Sens slope function
  calculate_sens_slope <- function(column, index) {
    if (sum(!is.na(column)) >= 5) {
      slope <- modifiedmk::mmky(column)[[2]]
      return(list(result = slope, index = index))
    } else {
      return(list(result = NA, index = index))
    }
  }
  
  # Use future_mapply for parallel Sens slope calculation
  unique_column_names <- paste0(sprintf("%08d", 1:ncol(tvalues_df)))
  colnames(tvalues_df) <- unique_column_names
  results_list <- future_mapply(calculate_sens_slope, tvalues_df, index = colnames(tvalues_df))
  
  # Reset to sequential processing
  plan(sequential)
  
  # Organize results
  results_df <- do.call(rbind, results_list)
  results_df <- results_df[order(results_df[, "index"]), ]
  results <- as.numeric(results_df[, "result"])
  
  # Create output raster
  e <- raster(xmn = min(lon), xmx = max(lon), ymn = min(lat), ymx = max(lat),
              ncols = length(unique(lon)), nrows = length(unique(lat)))
  raster::values(e) <- results
  
  # Define dimensions and create NetCDF output file
  xdim <- ncdim_def("lon", units = "degrees_east", longname = "Longitude", as.double(unique(lon)))
  ydim <- ncdim_def("lat", units = "degrees_north", longname = "Latitude", as.double(unique(lat)))
  fillvalue <- 1e32
  tmp_def <- ncvar_def("Sen", "slope", list(xdim, ydim), fillvalue, "Sen", prec = "single")
  
  setwd("")
  ncfname <- paste0("P_", groups, "_Sen_Summer_OC2pval", ".nc")
  ncout <- nc_create(ncfname, list(tmp_def), force_v4 = TRUE)
  
  # Save data to NetCDF file
  ncvar_put(ncout, tmp_def, getValues(e))
  ncatt_put(ncout, "lon", "axis", "X")
  ncatt_put(ncout, "lat", "axis", "Y")
  
  nc_close(ncout)
}

# Apply Create_New function for each group
groups_list <- c("chlor_a")
for (group in groups_list) {
  Create_New(group)
}
