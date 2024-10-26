# Load necessary libraries
library(caret)    # For machine learning utilities
library(dplyr)    # For data manipulation
library(randomForest)  # For random forest models

# Set working directory
setwd("C:/Users/haywarda/OneDrive - NIWA/Documents/Manuscripts_InProgress/Antarctic paper/")

# Load and preprocess data
data <- read.csv("RF_data_n2.csv") %>%
  na.omit() %>%                      # Remove rows with NA values
  filter(Depth < MLD)                # Filter rows where Depth is less than MLD

# Define target variables for modeling
variables <- c("Hapto", "Diatom", "Crypto", "Green", "Dino", "Pelago", "Syn")

# Initialize a list to store trained models
models <- list()

# Train a model for each variable
for (var in variables) {
  # Define cross-validation parameters
  cvIndex <- groupKFold(data$Voyage, k = 10,returnTrain = TRUE)
  train_control <- trainControl(method = "cv", number = 10, index = cvIndex)
  
  # Prepare dataset for modeling with the current target variable
  model_data <- data %>%
    select(all_of(var), Tchla, MLD, SSS, sst, PO4, ice, Alk, NO3, FeT) %>%
    rename(target = all_of(var))     # Rename target variable for formula compatibility
  
  # Train random forest model
  model <- train(
    target ~ ., data = model_data, method = 'rf', 
    ntree = 250, do.trace = 10, importance = TRUE, trControl = train_control
  )
  
  # Save the trained model
  saveRDS(model, file = paste0(var, "_Nature5.rds"))
  
  # Store the model in the models list
  models[[var]] <- model
}
