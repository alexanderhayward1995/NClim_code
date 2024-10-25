# Load necessary libraries
library(ncdf4, lib="Packages/")
library(raster)
library(lubridate, lib="Packages/")
library(tidyverse, lib="Packages/")
library(randomForest, lib="Packages/")
library(caret, lib="Packages/")

# Get command line arguments
args = commandArgs(trailingOnly = TRUE)  
print(args)

# Extract the date information from the command line argument
d <- strsplit(args, "_")[[1]][[1]]
d <- strsplit(d, "/")[[1]][[8]]
d <- substr(d, 2, nchar(d))  # Remove the leading character
e <- substr(d, 1, 7)          # Extract the date portion

# Convert the date to a Date object and format it
date <- as.Date(e, format = "%Y%j")  # Convert from YDDD format
nd <- format(date, "%Y%m%d")          # Reformat date
Mnth <- substr(nd, 5, 6)              # Extract month

# Load necessary NetCDF files
Diatom <- nc_open(args)
Green <- nc_open(paste("/nesi/project/niwa00020/Haywarda/Green/MO/", "P", d, "_Green.nc", sep=""))
Crypto <- nc_open(paste("/nesi/project/niwa00020/Haywarda/Crypto/MO/", "P", d, "_Crypto.nc", sep=""))
Syn <- nc_open(paste("/nesi/project/niwa00020/Haywarda/Syn/MO/", "P", d, "_Syn.nc", sep=""))
Pelago <- nc_open(paste("/nesi/project/niwa00020/Haywarda/Pelago/MO/", "P", d, "_Pelago.nc", sep=""))
Hapto <- nc_open(paste("/nesi/project/niwa00020/Haywarda/Hapto/MO/", "P", d, "_Hapto.nc", sep=""))
Dino <- nc_open(paste("/nesi/project/niwa00020/Haywarda/Dino/MO/", "P", d, "_Dino.nc", sep=""))

# Get latitude and longitude from one of the datasets
lat <- ncvar_get(Green, "lat")
lon <- ncvar_get(Green, "lon")

# Function to create a raster and return its values
createRasterValues <- function(Y, flip = FALSE) {
  r <- raster(t(Y), xmn=min(lon), xmx=max(lon), ymn=min(lat), ymx=max(lat))
  if (flip) {
    r <- flip(r, direction='y')  # Flip the raster if specified
  }
  return(getValues(r))  # Return the raster values
}

# Get values for various species
Diatom <- createRasterValues(ncvar_get(Diatom, "Diatom"))
Green <- createRasterValues(ncvar_get(Green, "Green"))
Crypto <- createRasterValues(ncvar_get(Crypto, "Crypto"))
Pelago <- createRasterValues(ncvar_get(Pelago, "Pelago"))
Hapto <- createRasterValues(ncvar_get(Hapto, "Hapto"))
Dino <- createRasterValues(ncvar_get(Dino, "Dino"))
Syn <- createRasterValues(ncvar_get(Syn, "Syn"))

# List of vectors for summation
vectors_list <- list(Diatom, Green, Crypto, Pelago, Hapto, Dino, Syn)

# Sum all vectors in the list
PSUM <- Reduce(+, vectors_list)

# Combine longitude, latitude, and sum values into a dataframe
s100 <- cbind(lon, lat, PSUM)

# Create a raster object for output
e <- raster(nrow=2160, ncol=4320, xmn=min(lon), xmx=max(lon), ymn=min(lat), ymx=max(lat))
values(e) <- PSUM  # Assign summed values to the raster

# Define dimensions for NetCDF
xdim <- ncdim_def("lon", units="degrees_east", longname="Longitude", as.double(unique(s100[, 1])))
ydim <- ncdim_def("lat", units="degrees_north", longname="Latitude", as.double(unique(s100[, 2])))

# Define fill value and variable for NetCDF
fillvalue <- 1e32
dlname <- "SUM"
tmp_def <- ncvar_def("SUM", "Chla ug/L", list(xdim, ydim), fillvalue, dlname, prec="single")

# Set working directory and create NetCDF file
setwd("/nesi/project/niwa00020/Haywarda/Sum")
ncfname <- paste("P", d, "_SUM", ".nc", sep="")
ncout <- nc_create(ncfname, list(tmp_def), force_v4=TRUE)

# Write data to NetCDF file
ncvar_put(ncout, tmp_def, values(e))

# Add metadata attributes
ncatt_put(ncout, "lon", "axis", "X")
ncatt_put(ncout, "lat", "axis", "Y")

# Close NetCDF file
nc_close(ncout)

# Repeat the process for different datasets (Hapto, Crypto, Syn, etc.)

# Function to handle processing for each species
processSpecies <- function(species_name, species_var) {
  # Open the corresponding NetCDF file
  species_data <- nc_open(paste0("/nesi/project/niwa00020/Haywarda/", species_name, "/MO/P", d, "_", species_name, ".nc"))
  
  # Get species values
  species_values <- createRasterValues(ncvar_get(species_data, species_var))
  
  # Calculate the proportion
  PSUM <- species_values / Sum  # Replace 'Sum' with the appropriate variable or data
  
  # Create a raster for output
  ks <- PSUM
  s100 <- cbind(lon, lat, ks)
  e <- raster(nrow=2160, ncol=4320, xmn=min(lon), xmx=max(lon), ymn=min(lat), ymx=max(lat))
  raster::values(e) <- ks
  
  # Define dimensions for NetCDF
  xdim <- ncdim_def("lon", units="degrees_east", longname="Longitude", as.double(unique(s100[, 1])))
  ydim <- ncdim_def("lat", units="degrees_north", longname="Latitude", as.double(unique(s100[, 2])))
  
  # Define fill value and variable for NetCDF
  fillvalue <- 1e32
  dlname <- paste0("P", species_name)
  tmp_def <- ncvar_def(dlname, "Chla ug/L", list(xdim, ydim), fillvalue, dlname, prec="single")
  
  # Set working directory for output
  setwd(paste0("/nesi/project/niwa00020/Haywarda/", dlname, "/MO/"))
  ncfname <- paste("P", d, "_", dlname, ".nc", sep="")
  
  # Create and write to NetCDF file
  ncout <- nc_create(ncfname, list(tmp_def), force_v4=TRUE)
  ncvar_put(ncout, tmp_def, values(e))
  ncatt_put(ncout, "lon", "axis", "X")
  ncatt_put(ncout, "lat", "axis", "Y")
  
  # Close NetCDF file
  nc_close(ncout)
}

# Process each species
species_list <- c("Hapto", "Crypto", "Syn", "Green", "Pelago", "Dino")  # Add other species as needed
for (species in species_list) {
  processSpecies(species, species)
}
