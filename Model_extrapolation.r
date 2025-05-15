# Load necessary libraries
library(ncdf4,lib="Packages/")
library(raster)
library(lubridate,lib="Packages/")
library(tidyverse,lib="Packages/")
library(randomForest,lib="Packages/")
library(caret,lib="Packages/")

# Parse command-line arguments
args <- commandArgs(trailingOnly = TRUE)
print(args)

args = commandArgs(trailingOnly = TRUE)  
print(args)
d <- strsplit(args,"_")
d <- d[[1]][[1]]
d <- strsplit(d,"/")
d <- d[[1]][[8]]
e <- substr(d,2,nchar(d))
f <- substr(e,1,7)

# Extract date from input arguments
date_str <- e
date <- as.Date(f, format = "%Y%j")
formatted_date <- format(date, "%Y%m%d")
month <- substr(formatted_date, 5, 6)

# Helper function to open netCDF variables
open_nc_var <- function(path, file_suffix) {
  nc_open(file.path("", path, paste0("E", date_str, file_suffix)))
}

# Open various netCDF files containing environmental variables
Chla <- nc_open(args)
Alk <- open_nc_var("ALK_9_nc", "_MO_ALK.nc")
FeT <- open_nc_var("FeT_9_nc", "_MO_FeT.nc")
SSS <- open_nc_var("SSS_9_nc", "_MO_SALTanom.nc")
NO3 <- open_nc_var("NO3_9_nc", "_MO_NO3.nc")
PO4 <- open_nc_var("PO4_9_nc", "_MO_PO4.nc")
pCO2 <- open_nc_var("pCO2_9_nc", "_MO_pCO2.nc")
MLD <- nc_open(paste("", "E", date_str, "_MO_Mld.nc", sep=""))
SST <- nc_open(paste("", "E", date_str, "_sst.nc", sep=""))

# Define helper function for processing environmental variables
Nuts <- function(Y){
  r <- raster(t(Y),xmn=min(lon),xmx=max(lon),ymn=min(lat),ymx=max(lat))
  xy <- data.frame(xyFromCell(r, 1:ncell(r)))
  v <- getValues(r)
  return(v)
}

ICE <- function(Y){
  Y[is.na(Y)] <- 0
  r <- raster(t(Y),xmn=min(lon),xmx=max(lon),ymn=min(lat),ymx=max(lat))
  xy <- data.frame(xyFromCell(r, 1:ncell(r)))
  v <- getValues(r)
  return(v)
}
# Load environmental data and process into variables
lat <- ncvar_get(Chla, "lat")
lon <- ncvar_get(Chla, "lon")
aChla <- Nuts(ncvar_get(Chla,"chlor_a"))
Alka <- Nuts(ncvar_get(Alk,"ALK"))
NO3a <- Nuts(ncvar_get(NO3,"NO3"))
FeTa <- Nuts(ncvar_get(FeT,"FeT"))
SSSa <- Nuts(ncvar_get(SSS,"SALTanom"))
pCO2a <- Nuts(ncvar_get(pCO2,"pCO2"))
PO4a <- Nuts(ncvar_get(PO4,"PO4"))
MLDa <- Nuts(ncvar_get(MLD,"MldBoyer"))
ICEa <- ICE(ncvar_get(SST,"sea_ice_fraction"))
SSTa <- Nuts(ncvar_get(SST,"analysed_sst"))

# Combine abiotic data into a single matrix
abiotic <- cbind(aChla, MLDa, SSSa, SSTa, PO4a, ICEa, Alka, NO3a, pCO2a, FeTa)

# Function to process and predict data using a model and save netCDF
create_nc_output <- function(N, dlname, folder_name, model_name) {
  df_clean <- na.omit(abiotic)
  colnames(df_clean) <- c('N.Tchla','N.MLD','N.SSS','N.sst','N.PO4','N.ice','N.Alk','N.NO3','N.pCO2','N.FeT')
  prediction <- predict(N, as.data.frame(df_clean))
  Month <- rep(as.numeric(1),length(SSTa))
  # Create grid with predicted data
    df <- as.data.frame(abiotic)
    ks <- Month

    for (i in 1:ncol(df)){
      m <- which(is.na(df[,i]))
      ks[m] <- NA
    }
    km <- rep(NA,length(SSTa))
    km[!is.na(ks)] <- prediction
    s100 <- cbind(lon,lat,km)
    e <- raster(nrow=2160,ncol=4320,xmn=min(lon),xmx=max(lon),ymn=min(lat),ymx=max(lat))
    raster::values(e) <- km
    fillvalue <- 1e32
  # Define dimensions and variables for netCDF
  xdim <- ncdim_def("lon", units = "degrees_east", longname = "Longitude", as.double(unique(lon)))
  ydim <- ncdim_def("lat", units = "degrees_north", longname = "Latitude", as.double(unique(lat)))
  tmp_def <- ncvar_def(dlname, "Chla ug/L", list(xdim, ydim), fillvalue, dlname, prec = "single")
  
  # Create folder for model outputs
  output_folder <- file.path("", folder_name, model_name)
  dir.create(output_folder, recursive = TRUE, showWarnings = FALSE)
  
  # Save netCDF file
  ncfname <- file.path(output_folder, paste0("P", date_str, "_", dlname, ".nc"))
  ncout <- nc_create(ncfname, list(tmp_def), force_v4 = TRUE)
  ncvar_put(ncout, tmp_def, getValues(e))
  
  # Add attributes
  ncatt_put(ncout, "lon", "axis", "X")
  ncatt_put(ncout, "lat", "axis", "Y")
  nc_close(ncout)
}

# Apply function to all models
organism_models <- c("crypto")

for (organism in organism_models) {
  for (i in 1:11) {
    model_name <- paste0(organism, "_rf_model", i, ".rds")
    folder_name <- organism
    N <- readRDS(model_name)
    create_nc_output(N, organism, paste0(folder_name,'_unc'), paste0("Model_", i))
  }
}

# Return to initial working directory
setwd("")
