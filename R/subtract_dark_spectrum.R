#' Subtract Dark Spectrum from Tube Spectrum for Baseline Correction
#' @author Christian L. Goueguel
#' @description
#' This function performs dark spectrum subtraction from tube spectrum measurements
#' to correct for baseline noise and systematic sensor artifacts.
#' @param tube_spectrum Data frame containing tube spectrum measurements. Should 
#'   include numeric columns representing spectral intensities at different 
#'   wavelengths or channels. May also contain non-numeric metadata columns 
#'   which will be preserved in the output.
#' @param dark_spectrum Data frame containing dark spectrum measurements (baseline).
#'   Should have matching numeric columns with tube_spectrum for proper subtraction.
#'   The dark spectrum represents measurements taken without sample illumination
#'   and captures instrumental noise and offset.
#' @param verbose Logical value (default: FALSE). When TRUE, the function prints
#'   diagnostic messages about data processing steps, including row count 
#'   adjustments and column mismatches.
#' @return A data frame with the same structure as tube_spectrum, where numeric
#'   spectral columns have had the corresponding dark spectrum values subtracted.
#'   Non-numeric columns and metadata are preserved unchanged. The corrected
#'   spectrum represents the true signal with baseline artifacts removed.
#' @details
#' The function performs the following operations:
#' 
#' Row Alignment: If the input data frames have different numbers of rows, the
#' function automatically trims both to the minimum row count to ensure proper
#' element-wise subtraction. This handles cases where acquisition timing or 
#' sampling rates differ between measurements.
#' 
#' Column Matching: The function identifies all numeric columns in the tube 
#' spectrum (excluding the 'sensor' column if present) and attempts to subtract
#' corresponding columns from the dark spectrum. Only matching columns are 
#' processed, allowing for flexibility in data structure.
#' 
#' Baseline Correction: For each matching numeric column, the function performs
#' element-wise subtraction (tube - dark) to remove baseline effects.
#' 
#' Data Integrity: Non-numeric columns and any columns not present in both 
#' data frames are preserved unchanged, maintaining metadata and auxiliary 
#' information throughout the processing pipeline.
#' @note
#' This function assumes that dark spectrum measurements were taken under 
#' identical conditions as the tube measurements except for the absence of 
#' sample signal. The quality of baseline correction depends on the stability
#' of the instrument during dark spectrum acquisition. For optimal results, dark 
#' spectra should be acquired immediately before or
#' after tube measurements to minimize temporal drift in instrumental response.
#' @export subtract_dark_spectrum
subtract_dark_spectrum <- function(tube_spectrum, dark_spectrum, verbose = FALSE) {
  if (missing(tube_spectrum)) {
    stop("Argument 'tube_spectrum' is required")
  }
  if (!is.data.frame(tube_spectrum)) {
    stop("Argument 'tube_spectrum' must be a data frame")
  }
  # if (nrow(tube_spectrum) == 0) {
  #   stop("Tube spectrum data frame is empty")
  # }
  if (missing(dark_spectrum)) {
    stop("Argument 'dark_spectrum' is required")
  }
  if (!is.data.frame(dark_spectrum)) {
    stop("Argument 'dark_spectrum' must be a data frame")
  }
  if (nrow(dark_spectrum) == 0) {
    stop("Dark spectrum data frame is empty")
  }
  if (!is.logical(verbose) || length(verbose) != 1) {
    stop("Argument 'verbose' must be a single logical value (TRUE or FALSE)")
  }
  
  if (nrow(tube_spectrum) != nrow(dark_spectrum)) {
    min_rows <- min(nrow(tube_spectrum), nrow(dark_spectrum))
    if (verbose) {
      cat(
        "Row count mismatch - Tube:",
        nrow(tube_spectrum),
        "Dark:",
        nrow(dark_spectrum),
        "- Trimming to",
        min_rows,
        "rows\n"
      )
    }
    if (!verbose) {
      warning(sprintf(
        "Row count mismatch detected. Trimming from Tube: %d, Dark: %d to %d rows",
        nrow(tube_spectrum), nrow(dark_spectrum), min_rows
      ))
    }
    tube_spectrum <- tube_spectrum[1:min_rows, ]
    dark_spectrum <- dark_spectrum[1:min_rows, ]
  }
  
  spectral_cols <- names(tube_spectrum)[sapply(tube_spectrum, is.numeric)]
  spectral_cols <- spectral_cols[spectral_cols != "sensor"]
  
  if (length(spectral_cols) == 0) {
    warning("No numeric spectral columns found in tube spectrum for correction")
    return(tube_spectrum)  
  }
  
  if (verbose) {
    cat(sprintf("Processing %d spectral columns for dark subtraction\n", 
                length(spectral_cols)))
  }
  
  corrected_spectrum <- tube_spectrum
  processed_cols <- character(0)
  skipped_cols <- character(0)
  
  for (col in spectral_cols) {
    if (col %in% names(dark_spectrum)) {
      if (length(tube_spectrum[[col]]) == length(dark_spectrum[[col]])) {
        corrected_spectrum[[col]] <- tube_spectrum[[col]] - dark_spectrum[[col]]
        processed_cols <- c(processed_cols, col)
        if (any(corrected_spectrum[[col]] < 0, na.rm = TRUE)) {
          neg_count <- sum(corrected_spectrum[[col]] < 0, na.rm = TRUE)
          if (verbose) {
            cat(sprintf("Warning: Column '%s' has %d negative values after correction\n", col, neg_count))
          }
        }
      } else {
        skipped_cols <- c(skipped_cols, col)
        if (verbose) {
          cat("Column", col, "length mismatch, skipping\n")
        }
      }
    } else {
      skipped_cols <- c(skipped_cols, col)
      if (verbose) {
        cat(sprintf("Column '%s' not found in dark spectrum, skipping\n", col))
      }
    }
  }
  
  if (verbose) {
    cat(sprintf("\nDark spectrum subtraction complete:\n"))
    cat(sprintf("  - Columns processed: %d\n", length(processed_cols)))
    cat(sprintf("  - Columns skipped: %d\n", length(skipped_cols)))
    if (length(skipped_cols) > 0 && length(skipped_cols) <= 10) {
      cat(sprintf("  - Skipped columns: %s\n", paste(skipped_cols, collapse = ", ")))
    } else if (length(skipped_cols) > 10) {
      cat(sprintf("  - First 10 skipped columns: %s, ...\n", paste(skipped_cols[1:10], collapse = ", ")))
    }
  }
  
  if (ncol(corrected_spectrum) != ncol(tube_spectrum)) {
    warning("Output column count differs from input - data structure may be corrupted")
  }
  return(corrected_spectrum)
}