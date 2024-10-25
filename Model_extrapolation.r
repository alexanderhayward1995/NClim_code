# Load necessary libraries
library(ncdf4)       # For handling netCDF data
library(raster)      # For raster data processing
library(lubridate)   # For date handling
library(tidyverse)   # For data manipulation
library(randomForest) # For machine learning models
library(caret)       # For machine learning workflows

# Parse command-line arguments
args <- commandArgs(trailingOnly = TRUE)
print(args)

# Extract date from input arguments
date_str <- substr(strsplit(args, "_")[[1]][[1]], 2, 8)
date <- as.Date(date_str, format = "%Y%j")
formatted_date <- format(date, "%Y%m%d")
month <- substr(formatted_date, 5, 6)

# Load the initial model file
model_file <- "Diatom_Nature5.rds"
N <- readRDS(model_file)

# Open various netCDF files containing environmental variables
open_nc_var <- function(path, file_suffix) {
  nc_open(file.path("/nesi/project/niwa00020/Haywarda/Ecco", path, paste0("E", date_str, file_suffix)))
}

Chla <- nc_open(args)
Alk <- open_nc_var("ALK_9_nc", "_MO_ALK.nc")
FeT <- open_nc_var("FeT_9_nc", "_MO_FeT.nc")
SSS <- open_nc_var("SSS_9_nc", "_MO_SALTanom.nc")
NO3 <- open_nc_var("NO3_9_nc", "_MO_NO3.nc")
PO4 <- open_nc_var("PO4_9_nc", "_MO_PO4.nc")
pCO2 <- open_nc_var("pCO2_9_nc", "_MO_pCO2.nc")
MLD <- nc_open(paste("/nesi/project/niwa00020/Haywarda/MLD/MO/", "E", date_str, "_MO_MldBoyer.nc", sep=""))
SST <- nc_open(paste("/nesi/project/niwa00020/Haywarda/sstcci/MO/", "E", date_str, "_sst.nc", sep=""))
aph_490 <- nc_open(paste("/nesi/project/niwa00020/Haywarda/IOPs/aph_490/", "O", date_str, "_aph490.nc", sep=""))

# Define helper functions for processing environmental variables

# Converts array to raster, flips vertically, returns values for each cell
process_raster_values <- function(Y, fill_na = TRUE) {
  if (fill_na) Y[is.na(Y)] <- 0
  r <- raster(t(Y), xmn = min(lon), xmx = max(lon), ymn = min(lat), ymx = max(lat))
  v <- getValues(flip(r, direction = 'y'))
  return(v)
}

# Load environmental data and process into variables
lat <- ncvar_get(Chla, "lat")
lon <- ncvar_get(Chla, "lon")
aChla <- process_raster_values(ncvar_get(Chla, "chlor_a"))
Alka <- process_raster_values(ncvar_get(Alk, "ALK"))
NO3a <- process_raster_values(ncvar_get(NO3, "NO3"))
FeTa <- process_raster_values(ncvar_get(FeT, "FeT"))
SSSa <- process_raster_values(ncvar_get(SSS, "SALTanom"))
pCO2a <- process_raster_values(ncvar_get(pCO2, "pCO2"))
PO4a <- process_raster_values(ncvar_get(PO4, "PO4"))
MLDa <- process_raster_values(ncvar_get(MLD, "MldBoyer"))
ICEa <- process_raster_values(ncvar_get(SST, "sea_ice_fraction"), fill_na = FALSE)
SSTa <- process_raster_values(ncvar_get(SST, "analysed_sst"))

# Combine abiotic data into a single matrix
abiotic <- cbind(aChla, MLDa, SSSa, SSTa, PO4a, ICEa, Alka, NO3a, pCO2a, FeTa)

# Function to process and predict data using machine learning model, then output netCDF
create_nc_output <- function(N, dlname, folder_name) {
  df_clean <- na.omit(abiotic)
  colnames(df_clean) <- N$finalModel$x
  prediction <- predict(N$finalModel, as.data.frame(df_clean))
  
  # Create grid with predicted data
  km <- rep(NA, length(SSTa))
  km[!is.na(abiotic[, 1])] <- prediction
  e <- raster(nrow = 2160, ncol = 4320, xmn = min(lon), xmx = max(lon), ymn = min(lat), ymx = max(lat))
  raster::values(e) <- km
  
  # Define dimensions and variables for netCDF
  xdim <- ncdim_def("lon", units = "degrees_east", longname = "Longitude", as.double(unique(lon)))
  ydim <- ncdim_def("lat", units = "degrees_north", longname = "Latitude", as.double(unique(lat)))
  tmp_def <- ncvar_def(dlname, "Chla ug/L", list(xdim, ydim), fillvalue = 1e32, dlname, prec = "single")
  
  # Create and save netCDF file
  ncfname <- paste0("P", date_str, "_", dlname, ".nc")
  setwd(paste0("/nesi/project/niwa00020/Haywarda/", folder_name, "/MO/"))
  ncout <- nc_create(ncfname, list(tmp_def), force_v4 = TRUE)
  ncvar_put(ncout, tmp_def, getValues(e))
  
  # Add attributes
  ncatt_put(ncout, "lon", "axis", "X")
  ncatt_put(ncout, "lat", "axis", "Y")
  nc_close(ncout)
}

# Apply the function to each organism model and store in netCDF files
model_list <- list(
  list(file = "Hapto_Nature5.rds", label = "Hapto", folder = "Hapto"),
  list(file = "Green_Nature5.rds", label = "Green", folder = "Green"),
  list(file = "Crypto_Nature5.rds", label = "Crypto", folder = "Crypto"),
  list(file = "Pelago_Nature5.rds", label = "Pelago", folder = "Pelago"),
  list(file = "Dino_Nature5.rds", label = "Dino", folder = "Dino"),
  list(file = "Syn_Nature5.rds", label = "Syn", folder = "Syn")
)

for (model in model_list) {
  N <- readRDS(model$file)
  create_nc_output(N, model$label, model$folder)
}

# Set working directory back to initial location
setwd("/nesi/project/niwa03483/Satellite_functions/")
