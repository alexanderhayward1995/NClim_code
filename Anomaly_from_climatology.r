# Load required packages
library(ncdf4, lib = "Packages/")
library(raster)
library(lubridate, lib = "Packages/")
library(tidyverse, lib = "Packages/")

# Set up command-line arguments for input
args = commandArgs(trailingOnly = TRUE)  
print(args)

# Extract date components from file path for labeling
d <- strsplit(args, "_")[[1]][1]
d <- strsplit(d, "/")[[1]][8]
d <- substr(d, 2, nchar(d))
e <- substr(d, 1, 7)

# Convert extracted date to a readable format
date <- as.Date(e, format = "%Y%j")  # Format as Date
nd <- format(date, "%Y%m%d")          # Format date as YYYYMMDD
Mnth <- substr(nd, 5, 6)              # Extract month as MM for lookup

# Define helper functions
get_month_names <- function(num) {
  # Get month name from month number (for file paths)
  months <- list(
    Jan = "01", Feb = "02", March = "03", April = "04", 
    May = "05", June = "06", July = "07", Aug = "08", 
    Sep = "09", Oct = "10", Nov = "11", Dec = "12"
  )
  month_names <- names(months)[unlist(which(sapply(months, function(x) num %in% x)))]
  return(month_names)
}

# Define data extraction function for netCDF
Nuts <- function(Y) {
  r <- raster(t(Y), xmn = min(lon), xmx = max(lon), ymn = min(lat), ymx = max(lat))
  xy <- data.frame(xyFromCell(r, 1:ncell(r)))
  v <- getValues(r)
  return(v)
}

# Initialize and open netCDF files for calculations
A <- nc_open(args)
lat <- ncvar_get(A, "lat")
lon <- ncvar_get(A, "lon")

# Define anomaly calculation and saving function
calculate_anomaly <- function(variable_name, variable_label) {
  # Open files for current and reference month
  B <- Nuts(ncvar_get(A, variable_name))  # Extract current data
  s2 <- get_month_names(Mnth)
  
  # Open reference month data
  C <- nc_open(paste0("/", s2, "/Complete/", variable_name, "_", s2, ".nc"))
  D <- Nuts(ncvar_get(C, variable_name))  # Extract reference data
  
  # Calculate anomaly as the difference
  nvals <- B - D
  
  # Create a new raster for storing the anomaly
  e <- raster(nrow = 2160, ncol = 4320, xmn = min(lon), xmx = max(lon), ymn = min(lat), ymx = max(lat))
  raster::values(e) <- nvals
  
  # Convert raster values to a format suitable for netCDF storage
  r <- e
  xy <- data.frame(xyFromCell(r, 1:ncell(r)))
  v <- getValues(r)
  
  # Define dimensions and variable details for netCDF
  xdim <- ncdim_def("lon", units = "degrees_east", longname = "Longitude", as.double(unique(xy[,1])))
  ydim <- ncdim_def("lat", units = "degrees_north", longname = "Latitude", as.double(unique(xy[,2])))
  fillvalue <- 1e32
  
  # Create netCDF variable and file for output
  tmp_def <- ncvar_def(variable_name, variable_label, list(xdim, ydim), fillvalue, variable_label, prec = "single")
  setwd(paste0("", variable_name, "/MA"))
  ncfname <- paste("P", d, "_", variable_name, "_Anom", ".nc", sep = "")
  
  ncout <- nc_create(ncfname, list(tmp_def), force_v4 = TRUE)
  ncvar_put(ncout, tmp_def, v)
  
  # Add coordinate axis attributes
  ncatt_put(ncout, "lon", "axis", "X")
  ncatt_put(ncout, "lat", "axis", "Y")
  
  nc_close(ncout)
}

# Calculate and save anomalies for each variable
variables <- list(
  PGreen = "Chla ug/L", PPelago = "Chla ug/L", 
  PDino = "Chla ug/L", Green = "Chla ug/L", 
  Pelago = "Chla ug/L", Dino = "Chla ug/L"
)

for (variable in names(variables)) {
  calculate_anomaly(variable, variables[[variable]])
}

# Close the main netCDF connection
nc_close(A)
setwd("")

