
-   [Milk Quality Spectroscopy](#milk-quality-spectroscopy)
-   [Dataset](#dataset)
    -   [Key data components](#key-data-components)
    -   [Important limitations](#important-limitations)
    -   [Data evolution](#data-evolution)

<!-- README.md is generated from README.Rmd. Please edit that file -->

# Milk Quality Spectroscopy

<!-- badges: start -->
<!-- badges: end -->

``` r
if (requireNamespace("ggplot2", quietly = TRUE)) {
  ggplot2::theme_set(ggplot2::theme_minimal(base_size = 12))
}
```

``` r
set.seed(123)
```

``` r
suppressPackageStartupMessages({
  library(tidyverse)
})
```

# Dataset

The data is structured by sensor/stall with each sensor directory
containing:

-   **Lab results** (CSV format with tube_no as primary key)
-   **Sensor configuration** (wavelengths, calibration coefficients)
-   **Spectral measurements** (Parquet files, one per milk sample)
-   **Dark spectra** (reference measurements)

## Key data components

**Lab Results**: Milk analysis including fat percentage, protein,
somatic cell count and, lactose linked to specific cows and milking
sessions.

**Spectral Data**: Raw 16-bit spectral arrays captured during milking,
with metadata like temperature, LED current, and integration time. Each
spectrum is classified as “dark,” “sample,” or “empty.”

**Sensor Info**: Each sensor has unique wavelength calibrations and
measurement parameters that aren’t standardized across sensors.

## Important limitations

-   Temperature and LED measurements are raw ADC values, not
    standardized between sensors
-   Each sensor measures different wavelengths
-   Some timing discrepancies remain due to clock source differences

## Data evolution

The dataset has evolved from initial CSV format to Parquet compression,
with spectral data consolidated into array columns rather than
individual wavelength columns for more efficient storage and processing.
