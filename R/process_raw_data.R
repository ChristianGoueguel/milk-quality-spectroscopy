#' Process Complete NIR Dataset from Multiple Sensors
#' @author Christian L. Goueguel
#' @description
#' This function provides a pipeline for processing Near-Infrared 
#' (NIR) spectroscopy data from multiple sensors. It orchestrates the complete 
#' workflow from raw data ingestion through baseline correction, handling both 
#' CSV metadata files and parquet spectral data files. The function processes 
#' data from 12 sensors, performs dark spectrum subtraction for baseline 
#' correction, and consolidates all results into organized data structures 
#' suitable for subsequent analysis.
#' @param verbose Logical value (default: FALSE). Controls the verbosity of 
#'   processing output. When FALSE, displays minimal progress indicators (dots) 
#'   for each sensor. When TRUE, provides detailed processing messages including 
#'   sensor-by-sensor progress and diagnostic information from subfunctions. 
#' @return A named list containing four consolidated data frames representing the complete processed dataset:
#'   \describe{
#'     \item{lab_results}{Combined laboratory results from all sensors. Includes 
#'     analytical reference measurements and metadata with sensor identification 
#'     columns for data provenance.}
#'     \item{sensor_data}{Combined sensor metadata from all sensors. Contains 
#'     operational parameters, timestamps, and measurement conditions for each sensor.}
#'     \item{dark_spectra}{Combined dark spectra from all sensors. Represents 
#'     baseline measurements taken without sample illumination, used for instrumental 
#'     noise characterization.}
#'     \item{corrected_spectra}{Combined baseline-corrected spectra from all tubes 
#'     across all sensors. Each spectrum has been corrected by subtracting the 
#'     corresponding dark spectrum, with tube numbers and sensor identifiers 
#'     preserved for traceability.}
#'   }
#' @details
#' The function executes a multi-stage processing pipeline that handles the 
#' complexity of multi-sensor NIR datasets. The processing workflow consists 
#' of several coordinated stages.
#' 
#' Data Structure Requirements: The function expects a specific directory 
#' structure under the global dataset_path variable, with folders named 
#' sensor_1 through sensor_12. Each sensor folder must contain standardized 
#' CSV files for metadata and parquet files for spectral data following the 
#' structure defined in read_sensor_csvs and read_sensor_parquets functions.
#' 
#' CSV Processing Stage: The function first reads laboratory results and sensor 
#' metadata from CSV files for all 12 sensors. These files typically contain 
#' reference measurements, sample identifications, and experimental conditions. 
#' The data from all sensors is consolidated using row binding operations that 
#' preserve sensor identification.
#' 
#' Parquet Processing Stage: Raw spectral data stored in parquet format is 
#' read for both dark spectra and tube measurements. The parquet format 
#' provides efficient storage for high-dimensional spectral data while 
#' maintaining numerical precision critical for spectroscopic analysis.
#' 
#' Baseline Correction Stage: For each sensor, the function applies dark 
#' spectrum subtraction to all tube measurements. This correction removes 
#' instrumental baseline effects, electronic noise, and systematic artifacts 
#' that would otherwise interfere with spectral analysis. The correction is 
#' performed tube by tube to maintain data integrity and traceability.
#' 
#' Data Consolidation: After processing, all corrected spectra are combined 
#' into a single data frame with proper tube number and sensor identification. 
#' This consolidated format facilitates subsequent multivariate analysis, 
#' calibration model development, and cross-sensor comparisons.
#' 
#' Memory Considerations: Processing data from 12 sensors with multiple tubes 
#' each can require substantial memory. The function loads all data into memory 
#' simultaneously, which may require systems with adequate RAM for large datasets.
#' @note
#' This function requires the global variable dataset_path to be defined before 
#' execution. The dataset_path should point to the root directory containing 
#' all sensor folders. The function assumes exactly 12 sensors numbered 1 through 12. 
#' All subsidiary functions (read_sensor_csvs, read_sensor_parquets, and 
#' subtract_dark_spectrum) must be loaded in the R environment before calling 
#' this function. Progress indicators are displayed during processing. In non-verbose mode, 
#' dots indicate sensor completion, while verbose mode provides detailed 
#' status updates throughout the pipeline.
#' @seealso
#' \code{\link{read_sensor_csvs}} for CSV data ingestion
#' \code{\link{read_sensor_parquets}} for parquet data ingestion
#' \code{\link{subtract_dark_spectrum}} for baseline correction
#' @export process_raw_data 
#' @importFrom purrr map map_dfr
process_raw_data <- function(verbose = FALSE) {
  if (!is.logical(verbose) || length(verbose) != 1) {
    stop("Argument 'verbose' must be a single logical value (TRUE or FALSE)")
  }
  if (!exists("dataset_path")) {
    stop("Global variable 'dataset_path' is not defined. ",
         "Please set it to the root directory containing sensor folders ",
         "before calling this function.")
  }
  if (!dir.exists(dataset_path)) {
    stop(sprintf("Dataset path does not exist or is not accessible: %s", 
                 dataset_path))
  }
  if (!requireNamespace("purrr", quietly = TRUE)) {
    stop("Package 'purrr' is required but not installed. ",
         "Please install it with: install.packages('purrr')")
  }
  
  required_functions <- c("read_sensor_csvs", "read_sensor_parquets", "subtract_dark_spectrum")
  missing_functions <- required_functions[!sapply(required_functions, exists)]
  if (length(missing_functions) > 0) {
    stop(sprintf("Required function(s) not found: %s. ",
                 paste(missing_functions, collapse = ", ")),
         "Please ensure all supporting functions are loaded.")
  }
  
  sensor_folders <- paste0("sensor_", 1:12)
  missing_sensors <- character()
  for (sensor in sensor_folders) {
    sensor_path <- file.path(dataset_path, sensor)
    if (!dir.exists(sensor_path)) {
      missing_sensors <- c(missing_sensors, sensor)
    }
  }
  
  if (length(missing_sensors) > 0) {
    warning(sprintf("Missing sensor folder(s): %s. ", paste(missing_sensors, collapse = ", ")), "Processing will continue with available sensors.")
    sensor_folders <- setdiff(sensor_folders, missing_sensors)
    if (length(sensor_folders) == 0) {
      stop("No sensor folders found for processing. ", "Please check the dataset structure.")
    }
  }
  
  cat(sprintf("Beginning NIR data processing for %d sensor(s)\n", length(sensor_folders)))
  if (verbose) {
    cat(sprintf("Sensors to process: %s\n", paste(sensor_folders, collapse = ", ")))
  }
  
  cat("Processing CSV files...\n")
  csv_data <- list()
  csv_read_errors <- character()
  
  for (sensor in sensor_folders) {
    tryCatch({
      csv_data[[sensor]] <- read_sensor_csvs(sensor)
      if (verbose) {
        cat(sprintf("  Successfully read CSV files for %s\n", sensor))
      }
    }, error = function(e) {
      csv_read_errors <<- c(csv_read_errors, sensor)
      warning(sprintf("Failed to read CSV files for %s: %s", 
                      sensor, e$message))
      csv_data[[sensor]] <<- list(lab_results = NULL, sensor_data = NULL)
    })
  }
  
  if (length(csv_read_errors) > 0) {
    warning(sprintf("Failed to read CSV files for %d sensor(s): %s", length(csv_read_errors), paste(csv_read_errors, collapse = ", ")))
  }
  
  names(csv_data) <- sensor_folders
  all_lab_results <- map_dfr(csv_data, function(x) {
    if (!is.null(x$lab_results)) x$lab_results else data.frame()
  })
  all_sensor_data <- map_dfr(csv_data, function(x) {
    if (!is.null(x$sensor_data)) x$sensor_data else data.frame()
  })
  
  if (nrow(all_lab_results) == 0) {
    warning("No laboratory results data was successfully read from CSV files")
  }
  if (nrow(all_sensor_data) == 0) {
    warning("No sensor data was successfully read from CSV files")
  }
  
  cat("Processing parquet files...\n")
  parquet_data <- list()
  parquet_read_errors <- character()
  
  for (sensor in sensor_folders) {
    tryCatch({
      parquet_data[[sensor]] <- read_sensor_parquets(sensor)
      if (verbose) {
        cat(sprintf("  Successfully read parquet files for %s\n", sensor))
      }
    }, error = function(e) {
      parquet_read_errors <<- c(parquet_read_errors, sensor)
      warning(sprintf("Failed to read parquet files for %s: %s", 
                      sensor, e$message))
      parquet_data[[sensor]] <<- list(dark_spectra = NULL, 
                                      tube_spectra = list())
    })
  }
  
  if (length(parquet_read_errors) > 0) {
    warning(sprintf("Failed to read parquet files for %d sensor(s): %s", length(parquet_read_errors), paste(parquet_read_errors, collapse = ", ")))
  }
  
  names(parquet_data) <- sensor_folders
  all_dark_spectra <- map_dfr(parquet_data, function(x) {
    if (!is.null(x$dark_spectra)) x$dark_spectra else data.frame()
  })
  
  if (nrow(all_dark_spectra) == 0) {
    warning("No dark spectra data was successfully read from parquet files")
  }
  
  corrected_spectra_list <- list()
  cat("Performing baseline correction...\n")
  
  for (sensor in sensor_folders) {
    if (verbose) {
      cat(paste("Processing", sensor, "...\n"))
    }
    if (is.null(parquet_data[[sensor]]$dark_spectra) || 
        length(parquet_data[[sensor]]$tube_spectra) == 0) {
      if (verbose) {
        cat(sprintf("  Skipping %s - insufficient spectral data\n", sensor))
      }
      corrected_spectra_list[[sensor]] <- list()
      next
    }
    dark_spectrum <- parquet_data[[sensor]]$dark_spectra
    tube_spectra <- parquet_data[[sensor]]$tube_spectra
    corrected_tubes <- map(
      tube_spectra,
      ~ subtract_dark_spectrum(.x, dark_spectrum, verbose = verbose)
    )
    for (i in seq_along(corrected_tubes)) {
      tube_file <- names(corrected_tubes)[i]
      tube_number <- gsub("tube_no_(\\d+)\\.parquet", "\\1", tube_file)
      corrected_tubes[[i]]$tube_number <- tube_number
    }
    corrected_spectra_list[[sensor]] <- corrected_tubes
    if (verbose) {
      cat(sprintf("  Processed %d tube(s) for %s\n", length(corrected_tubes), sensor))
    }
  }

  cat("\n")
  all_corrected_spectra <- map_dfr(corrected_spectra_list, ~ map_dfr(.x, identity))
  
  if (nrow(all_corrected_spectra) == 0) {
    warning("No corrected spectra were generated. ", "Check input data and processing steps.")
  } else {
    if (verbose) {
      cat(sprintf("Successfully processed %d corrected spectra records\n", nrow(all_corrected_spectra)))
    }
  }
  
  cat("Data processing complete!\n")
  if (verbose) {
    cat("\nProcessing Summary:\n")
    cat(sprintf("  Lab results records: %d\n", nrow(all_lab_results)))
    cat(sprintf("  Sensor data records: %d\n", nrow(all_sensor_data)))
    cat(sprintf("  Dark spectra records: %d\n", nrow(all_dark_spectra)))
    cat(sprintf("  Corrected spectra records: %d\n", nrow(all_corrected_spectra)))
    if (nrow(all_corrected_spectra) > 0) {
      n_sensors <- length(unique(all_corrected_spectra$sensor))
      n_tubes <- length(unique(paste(all_corrected_spectra$sensor, 
                                     all_corrected_spectra$tube_number)))
      cat(sprintf("  Unique sensors processed: %d\n", n_sensors))
      cat(sprintf("  Total tubes processed: %d\n", n_tubes))
    }
  }
  return(
    list(
      lab_results = all_lab_results,
      sensor_data = all_sensor_data,
      dark_spectra = all_dark_spectra,
      corrected_spectra = all_corrected_spectra
    )
  )
}