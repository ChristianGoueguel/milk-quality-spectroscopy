#' Read Sensor CSV Files from a Specified Directory
#' @author Christian L. Goueguel
#' @description
#' This function reads laboratory results and sensor data CSV files from a 
#' specified sensor folder within the dataset path. It adds a sensor identifier
#' column to both datasets for tracking purposes.
#' @param sensor_folder Character string. The name of the sensor folder containing
#'   the CSV files to be read. This folder should be located within the `dataset_path`
#'   directory and contain both "lab_results.csv" and "sensor.csv" files.
#' @return A named list containing two data frames:
#'   \describe{
#'     \item{lab_results}{Data frame containing laboratory results with an added 'sensor' column identifying the source sensor}
#'     \item{sensor_data}{Data frame containing sensor measurements with an added 'sensor' column identifying the source sensor}
#'   }
#' @details
#' The function expects the following file structure:
#' - dataset_path/
#'   - sensor_folder/
#'     - lab_results.csv
#'     - sensor.csv
#'     
#' Both CSV files are read with column type inference suppressed to avoid
#' parsing messages. The sensor folder name is appended as a new column to
#' both datasets for data provenance tracking.
#' @note
#' - Requires the \code{readr} package for the \code{read_csv} function
#' - The global variable \code{dataset_path} must be defined before calling this function
#' - Both CSV files must exist in the specified sensor folder
#' @export read_sensor_csvs
#' @importFrom readr read_csv
#'
read_sensor_csvs <- function(sensor_folder) {
  if (missing(sensor_folder)) {
    stop("Argument 'sensor_folder' is required")
  }
  if (!is.character(sensor_folder) || length(sensor_folder) != 1) {
    stop("Argument 'sensor_folder' must be a single character string")
  }
  if (!exists("dataset_path")) {
    stop("Global variable 'dataset_path' is not defined. ",
         "Please set it before calling this function.")
  }
  sensor_path <- file.path(dataset_path, sensor_folder)
  
  if (!dir.exists(sensor_path)) {
    stop(sprintf("Sensor folder does not exist: %s", sensor_path))
  }
  
  lab_results_file <- file.path(sensor_path, "lab_results.csv")
  sensor_data_file <- file.path(sensor_path, "sensor.csv")
  
  if (!file.exists(lab_results_file)) {
    stop(sprintf("Required file 'lab_results.csv' not found in: %s", sensor_path))
  }
  if (!file.exists(sensor_data_file)) {
    stop(sprintf("Required file 'sensor.csv' not found in: %s", sensor_path))
  }
  
  lab_results <- read_csv(
    lab_results_file, 
    show_col_types = FALSE
  )
  
  lab_results$sensor <- sensor_folder
  
  sensor_data <- read_csv(
    sensor_data_file, 
    show_col_types = FALSE
  )
  
  sensor_data$sensor <- sensor_folder
  
  if (nrow(lab_results) == 0) {
    warning(sprintf("Lab results file is empty for sensor: %s", sensor_folder))
  }
  if (nrow(sensor_data) == 0) {
    warning(sprintf("Sensor data file is empty for sensor: %s", sensor_folder))
  }
  
  return(
    list(
      lab_results = lab_results, 
      sensor_data = sensor_data
    )
  )
}