#' Read Sensor Parquet Files from a Specified Directory
#' @author Christian L. Goueguel
#' @description
#' This function reads dark spectra and tube spectra parquet files from a 
#' specified sensor folder within the dataset path. It processes multiple tube
#' files and adds sensor identifiers to all datasets for tracking purposes.
#' The function is designed to handle spectroscopic data organized in a 
#' standardized directory structure.
#' @param sensor_folder Character string. The name of the sensor folder containing
#'   the parquet files to be read. This folder should be located within the 
#'   dataset_path directory and follow the expected structure with dark_spectra
#'   and spectra subdirectories.
#'
#' @return A named list containing two elements:
#'   \describe{
#'     \item{dark_spectra}{Data frame containing dark spectra measurements with an added 'sensor' column identifying the source sensor}
#'     \item{tube_spectra}{Named list of data frames, where each element represents a tube's spectral data. Names correspond to the original file names (e.g., "tube_no_1.parquet"). Each data frame includes an added 'sensor' column.}
#'   }
#' @details
#' The function expects the following file structure:
#' - dataset_path/
#'   - sensor_folder/
#'     - dark_spectra/
#'       - dark_spectra.parquet
#'     - spectra/
#'       - tube_no_1.parquet
#'       - tube_no_2.parquet
#'       - ... (additional tube files following the pattern)
#' 
#' The function automatically detects all tube files matching the pattern
#' "tube_no_\\d+\\.parquet" in the spectra directory. Each file is read into
#' memory and stored in a named list for easy access and processing.
#'
#' @note
#' - Requires the \code{arrow} package for the \code{read_parquet} function
#' - Requires the \code{purrr} package for the \code{map} function
#' - The global variable \code{dataset_path} must be defined before calling this function
#' - The dark_spectra.parquet file must exist in the expected location
#' - At least one tube file should exist in the spectra directory for meaningful results
#' - Memory usage can be significant when processing many large parquet files
#' @export read_sensor_parquets
#' @importFrom arrow read_parquet
#' @importFrom purrr map
#'
read_sensor_parquets <- function(sensor_folder) {
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
  if (!requireNamespace("arrow", quietly = TRUE)) {
    stop("Package 'arrow' is required but not installed. ",
         "Please install it with: install.packages('arrow')")
  }
  if (!requireNamespace("purrr", quietly = TRUE)) {
    stop("Package 'purrr' is required but not installed. ",
         "Please install it with: install.packages('purrr')")
  }
  
  sensor_path <- file.path(dataset_path, sensor_folder)
  
  if (!dir.exists(sensor_path)) {
    stop(sprintf("Sensor folder does not exist: %s", sensor_path))
  }
  
  dark_spectra_path <- file.path(sensor_path, "dark_spectra", "dark_spectra.parquet")
  
  if (!file.exists(dark_spectra_path)) {
    stop(sprintf("Required dark spectra file not found: %s", dark_spectra_path))
  }
  
  spectra_path <- file.path(sensor_path, "spectra")
  
  if (!dir.exists(spectra_path)) {
    stop(sprintf("Required spectra directory not found: %s", spectra_path))
  }
  
  tryCatch({
    dark_spectra <- read_parquet(dark_spectra_path)
  }, error = function(e) {
    stop(sprintf("Failed to read dark spectra file: %s\nError: %s", dark_spectra_path, e$message))
  })
  
  if (!is.data.frame(dark_spectra)) {
    stop("Dark spectra file did not return a valid data frame")
  }
  if (nrow(dark_spectra) == 0) {
    warning(sprintf("Dark spectra file is empty for sensor: %s", sensor_folder))
  }
  
  dark_spectra$sensor <- sensor_folder
  
  tube_files <- list.files(
    spectra_path, 
    pattern = "tube_no_\\d+\\.parquet$", 
    full.names = TRUE
  )
  
  if (length(tube_files) == 0) {
    warning(sprintf("No tube files found in spectra directory: %s", spectra_path))
    return(list(
      dark_spectra = dark_spectra, 
      tube_spectra = list()
    ))
  }
  
  message(sprintf("Found %d tube file(s) in sensor folder: %s", length(tube_files), sensor_folder))
  
  tube_spectra_list <- map(tube_files, function(file_path) {
    tryCatch({
      read_parquet(file_path)
    }, error = function(e) {
      warning(sprintf("Failed to read tube file: %s\nError: %s", 
                      file_path, e$message))
      return(NULL)
    })
  })
  
  failed_reads <- sum(sapply(tube_spectra_list, is.null))
  if (failed_reads > 0) {
    warning(sprintf("Failed to read %d tube file(s)", failed_reads))
    tube_spectra_list <- tube_spectra_list[!sapply(tube_spectra_list, is.null)]
  }
  
  names(tube_spectra_list) <- basename(tube_files[!sapply(tube_spectra_list, is.null)])
  
  tube_spectra_list <- map(tube_spectra_list, function(tube_data) {
    if (!is.data.frame(tube_data)) {
      warning("Tube file did not return a valid data frame, skipping sensor assignment")
      return(tube_data)
    }
    tube_data$sensor <- sensor_folder
    return(tube_data)
  })
  
  if (length(tube_spectra_list) == 0 && nrow(dark_spectra) == 0) {
    warning(sprintf("No valid spectral data found for sensor: %s", sensor_folder))
  }
  
  return(
    list(
      dark_spectra = dark_spectra, 
      tube_spectra = tube_spectra_list
    )
  )
}