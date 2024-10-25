# Load required libraries from custom library paths if necessary
library(R.matlab, lib = "Packages/")
library(ncdf4, lib = "Packages/")
library(raster)
library(terra)
library(lubridate, lib = "Packages/")
library(tidyverse, lib = "Packages/")
library(hutilscpp, lib = "Packages/")
library(plyr, lib = "Packages/")

# Set working directory to location of NetCDF files
setwd("/nesi/project/niwa00020/Haywarda/Ecco/NO3_9_nc/")
# List all NetCDF files in directory
temp <- list.files(pattern = "*.nc")

# Extract date information from file names, convert to numeric, and sort files by date
split <- strsplit(temp, "_") 
split <- as.numeric(sapply(split, function(x) sub("E", "", x[1])))
myFiles.correct.order <- temp[order(split)]

# Convert Julian day format to standard date format (YYYYMMDD)
split <- substr(split, 1, 7)
date <- as.Date(split, format = "%Y%j")
split <- format(date, "%Y%m%d")

# Load metadata file and organize data by date
setwd("/nesi/project/niwa03483/Satellite_functions/")
MD <- read.csv("ODV_LatLon_Nov.csv", header = TRUE)
D <- as.numeric(MD$Date)
Date <- MD[[5]]
MD <- data.frame(MD)
MD2 <- MD[order(as.numeric(MD[[5]]), decreasing = FALSE), ]

# Extract Year and Month for grouping, add column to metadata
split <- substr(split, 1, 6)
D2 <- substr(MD2[, 5], 1, 6)
D2 <- as.vector(D2)
MD2$YearMonth <- D2

# Find and extract unique indices for each month/year combination
b <- list()
for (i in D2) {
  b[[length(b) + 1]] <- which(split == i)
}
b <- unique(unlist(b))

# Set working directory back to source files, open and load NetCDF data
setwd("/nesi/project/niwa00020/Haywarda/Ecco/NO3_9_nc/")
Res <- lapply(myFiles.correct.order[b], function(i) {
  nc_open(i)
})

# Extract latitude and longitude data from first file for later use
lat <- ncvar_get(Res[[1]], "lat")
lon <- ncvar_get(Res[[2]], "lon")

# Define functions for data extraction and processing
Biooptical <- function(Y) {
  Y <- ncvar_get(Y, "sea_ice_fraction")
  r <- raster(t(Y), xmn = min(lon), xmx = max(lon), ymn = min(lat), ymx = max(lat))
  r <- flip(r, direction = 'y')
  xy <- data.frame(xyFromCell(r, 1:ncell(r)))
  v <- getValues(r)
  v <- ifelse(v == 0, NA, v)
  cbind(xy, v)
}

Nuts <- function(Y) {
  lat <- ncvar_get(Y, "lat")
  lon <- ncvar_get(Y, "lon")
  Y <- ncvar_get(Y, "NO3")
  r <- raster(t(Y), xmn = min(lon), xmx = max(lon), ymn = min(lat), ymx = max(lat))
  crs(r) <- CRS("+init=epsg:4087")
  r <- projectExtent(r, CRS("+init=epsg:4087"))
  values(r) <- t(Y)
  xy <- data.frame(xyFromCell(r, 1:ncell(r)))
  v <- getValues(r)
  v <- ifelse(v == 0, NA, v)
  cbind(xy, v)
}

ICE <- function(Y) {
  Y <- ncvar_get(Y, "sea_ice_fraction")
  Y[is.na(Y)] <- 0
  r <- raster(t(Y), xmn = min(lon), xmx = max(lon), ymn = min(lat), ymx = max(lat))
  xy <- data.frame(xyFromCell(r, 1:ncell(r)))
  v <- getValues(r)
  v <- ifelse(v == 0, NA, v)
  cbind(xy, v)
}

# Apply the nutrient extraction function to the list of opened NetCDF data
Res <- lapply(Res, Nuts)

# Create lists of indices by year/month combination and subset metadata
P <- list()
for (i in split[b]) {
  P[[length(P) + 1]] <- which(MD2$YearMonth == i)
}

MD2a <- MD2[unlist(P), ]

# Create lists for chlorophyll data and location coordinates
chls <- list()
for (i in 1:length(P)) {
  chls[[length(chls) + 1]] <- MD2[P[[i]], ]
}

latlons <- list()
for (i in 1:length(chls)) {
  latlons[[length(latlons) + 1]] <- as.matrix(cbind(chls[[i]]$Lon, chls[[i]]$Lat))
}

# Define extractor functions for spatial data extraction and neighbor-based averaging
single_extractor <- function(Resd) {
  gr <- list()
  for (i in 1:length(Resd)) {
    e <- raster(nrow = 2160, ncol = 4320, xmn = min(Resd[[i]][, 1]), xmx = max(Resd[[i]][, 1]), ymn = min(Resd[[i]][, 2]), ymx = max(Resd[[i]][, 2]))
    values(e) <- Resd[[i]][, 3]
    gr[[length(gr) + 1]] <- terra::extract(e, latlons[[i]])
  }
  return(gr)
}

queen_extractor <- function(Resd) {
  gr <- list()
  ress <- list()
  for (i in 1:length(Resd)) {
    e <- raster(nrow = 2160, ncol = 4320, xmn = min(Resd[[i]][, 1]), xmx = max(Resd[[i]][, 1]), ymn = min(Resd[[i]][, 2]), ymx = max(Resd[[i]][, 2]))
    values(e) <- Resd[[i]][, 3]
    cells <- lapply(1:nrow(latlons[[i]]), function(j) cellFromXY(e, latlons[[i]][j, 1:2]))
    neighbors <- lapply(cells, function(cell) adjacent(e, cell, directions = 8, include = TRUE))
    ress[[length(ress) + 1]] <- sapply(neighbors, function(cell) mean(terra::extract(e, cell[, 2]), na.rm = TRUE))
  }
  return(ress)
}

# Perform single-cell and neighbor-averaged extractions, combine results
gr <- single_extractor(Res)
boxed <- queen_extractor(Res)

# Combine results with metadata and save as CSV
single <- unlist(gr)
boxed <- unlist(boxed)
testering <- cbind(MD2a, single, boxed)
setwd("/nesi/project/niwa03483/Satellite_functions/")
write.csv(testering, "New_NO3.csv")
