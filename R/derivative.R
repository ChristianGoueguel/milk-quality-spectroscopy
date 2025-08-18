#' Calculate Savitzky-Golay Derivatives for Spectral Data
#' @author Christian L. Goueguel
#' @param x Data frame containing spectral data and metadata
#' @param window_size Integer, filter window size (must be odd, >= 3). Default: 11
#' @param poly_order Integer, polynomial order for fitting (< window_size). Default: 2  
#' @param deriv Integer, derivative order (0 = smoothing, 1 = first derivative, 2 = second derivative). Default: 1
#' @param spectral_cols Character vector or selection helper for spectral columns. Default: starts_with("X")
#' @param preserve_metadata Logical, whether to keep non-spectral columns. Default: TRUE
#' @param prefix Character, prefix for output column names. Default: paste0("d", deriv, "_")
#' 
#' @return Data frame with derivative spectra and optionally preserved metadata
#' @examples
#' # First derivative with default parameters
#' first_deriv <- derivatives(spectra_data)
#' 
#' # Second derivative with larger window
#' second_deriv <- derivatives(spectra_data, window_size = 15, deriv = 2)
#' 
#' # Smoothing only (0th derivative)
#' smoothed <- derivatives(spectra_data, deriv = 0, prefix = "smooth_")

derivative <- function(
    x, 
    window_size = 11, 
    poly_order = 2, 
    deriv = 1,
    spectral_cols = starts_with("X"),
    preserve_metadata = TRUE,
    prefix = paste0("d", deriv, "_")
    ) {
  if (!requireNamespace("signal", quietly = TRUE)) {
    stop("Package 'signal' is required but not installed.")
  }
  library(signal, quietly = TRUE)
  
  if (!is.data.frame(x)) {
    stop("Input 'x' must be a data frame.")
  }
  
  if (window_size < 3 || window_size %% 2 == 0) {
    stop("window_size must be odd and >= 3")
  }
  
  if (poly_order >= window_size) {
    stop("poly_order must be less than window_size")
  }
  
  if (deriv < 0 || deriv > poly_order) {
    stop("deriv must be between 0 and poly_order")
  }
  
  tryCatch({
    spectral_data <- x %>%
      select({{ spectral_cols }})
  }, error = function(e) {
    stop("Could not select spectral columns. Check spectral_cols parameter.")
  })
  
  if (ncol(spectral_data) == 0) {
    stop("No spectral columns found with the specified selection.")
  }
  
  non_numeric <- spectral_data %>%
    select(!where(is.numeric)) %>%
    names()
  
  if (length(non_numeric) > 0) {
    warning("Non-numeric columns found in spectral data: ", 
            paste(non_numeric, collapse = ", "))
    spectral_data <- spectral_data %>% select(where(is.numeric))
  }
  
  x_matrix <- as.matrix(spectral_data)
  
  if (any(is.na(x_matrix))) {
    warning("Missing values detected in spectral data. Results may be unreliable.")
  }
  
  if (ncol(x_matrix) < window_size) {
    stop("Number of spectral channels (", ncol(x_matrix), 
         ") is less than window_size (", window_size, ")")
  }
  
  apply_savgol <- function(spectrum, deriv_order) {
    tryCatch({
      savgol(spectrum, fl = window_size, forder = poly_order, dorder = deriv_order)
    }, error = function(e) {
      warning("Error applying Savitzky-Golay filter to spectrum: ", e$message)
      rep(NA, length(spectrum))
    })
  }
  
  derivative_matrix <- t(apply(x_matrix, 1, function(spectrum) {
    apply_savgol(spectrum, deriv)
  }))
  
  original_names <- colnames(spectral_data)
  if (is.null(original_names)) {
    original_names <- paste0("X", seq_len(ncol(spectral_data)))
  }
  
  new_names <- original_names
  colnames(derivative_matrix) <- new_names
  derivative_df <- as.data.frame(derivative_matrix)
  
  if (preserve_metadata) {
    metadata <- x %>% 
      select(!{{ spectral_cols }})
    
    if (ncol(metadata) > 0) {
      result <- bind_cols(metadata, derivative_df)
    } else {
      result <- derivative_df
    }
  } else {
    result <- derivative_df
  }
  
  attr(result, "preprocessing") <- list(
    method = "savitzky_golay_derivative",
    window_size = window_size,
    poly_order = poly_order,
    derivative_order = deriv,
    timestamp = Sys.time()
  )
  
  return(result)
}

#' Apply Multiple Derivative Orders
#' @author Christian L. Goueguel
#' @param x Data frame containing spectral data
#' @param orders Integer vector of derivative orders to calculate
#' @param ... Additional parameters passed to derivatives()
#' 
#' @return Data frame with all requested derivatives
#' @examples
#' # Calculate 1st and 2nd derivatives
#' multi_deriv <- derivatives_multi(spectra_data, orders = c(1, 2))

derivatives_multi <- function(x, orders = c(1, 2), ...) {
  
  # Validate orders
  if (!all(orders >= 0)) {
    stop("All derivative orders must be >= 0")
  }
  
  result_list <- list()
  
  for (order in orders) {
    cat("\n--- Processing derivative order", order, "---\n")
    
    deriv_result <- derivatives(x, deriv = order, 
                                preserve_metadata = (order == min(orders)), ...)
    
    if (order == min(orders)) {
      # Keep metadata for first result
      result_list[[paste0("order_", order)]] <- deriv_result
    } else {
      # Only spectral columns for subsequent results
      spectral_cols <- deriv_result %>% 
        select(starts_with(paste0("d", order, "_")))
      result_list[[paste0("order_", order)]] <- spectral_cols
    }
  }
  
  # Combine all results
  final_result <- bind_cols(result_list)
  
  # Add comprehensive attributes
  attr(final_result, "preprocessing") <- list(
    method = "savitzky_golay_multiple_derivatives",
    derivative_orders = orders,
    timestamp = Sys.time()
  )
  
  return(final_result)
}

#' Check Preprocessing History
#'
#' @param x Data frame with preprocessing attributes
#' @return Prints preprocessing information

check_preprocessing <- function(x) {
  preproc_info <- attr(x, "preprocessing")
  
  if (is.null(preproc_info)) {
    cat("No preprocessing information found.\n")
  } else {
    cat("Preprocessing Information:\n")
    cat("Method:", preproc_info$method, "\n")
    
    if (!is.null(preproc_info$derivative_order)) {
      cat("Derivative order:", preproc_info$derivative_order, "\n")
    }
    if (!is.null(preproc_info$derivative_orders)) {
      cat("Derivative orders:", paste(preproc_info$derivative_orders, collapse = ", "), "\n")
    }
    if (!is.null(preproc_info$window_size)) {
      cat("Window size:", preproc_info$window_size, "\n")
    }
    if (!is.null(preproc_info$poly_order)) {
      cat("Polynomial order:", preproc_info$poly_order, "\n")
    }
    
    cat("Processed on:", as.character(preproc_info$timestamp), "\n")
  }
}