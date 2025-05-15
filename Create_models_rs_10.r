# Load necessary libraries
library(caret)       # For machine learning utilities
library(dplyr)       # For data manipulation
library(randomForest)  # For random forest models

# Set working directory
#omitted for security

# Load and preprocess data
data <- read.csv("RF_data_n2.csv") %>%
  na.omit() %>%                      # Remove rows with NA values
  filter(Depth < MLD)                # Filter rows where Depth is less than MLD

# Define target variables for modeling
variables <- c("Hapto", "Diatom", "Crypto", "Green", "Dino", "Pelago", "Syn")

# Define a vector of 10 distinct random seeds
# (Here we just use 1:10, but you could sample from a larger range if you like)
seeds <- 1:10

# Loop over each response variable
for (var in variables) {
  # Prepare the predictor + target dataset once
  model_data <- data %>%
    select(all_of(var), Tchla, MLD, SSS, sst, PO4, ice, Alk, NO3, FeT) %>%
    rename(target = all_of(var))
  
  # For each seed, train and save a separate model
  for (seed in seeds) {
    set.seed(seed)  # ensure reproducibility
    
    # Create stratified folds by Voyage:
    cvIndex <- groupKFold(data$Voyage, k = 10, returnTrain = TRUE)
    train_control <- trainControl(
      method    = "cv",
      number    = 10,
      index     = cvIndex,
      search    = "grid"
    )
    
    # Train the RF model
    model <- train(
      target ~ .,
      data      = model_data,
      method    = "rf",
      ntree     = 250,
      do.trace  = 10,
      importance= TRUE,
      trControl = train_control
    )
    
    # Build a filename that encodes both variable name and seed
    fname <- paste0(var, "_Nature5_seed", sprintf("%02d", seed), ".rds")
    saveRDS(model, file = fname)
  }
}
