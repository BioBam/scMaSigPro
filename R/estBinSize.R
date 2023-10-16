#' estBinSize
#'
#' @description
#' This function calculates the optimal bin size for discretizing a continuous variable. It uses various bin_methods for bin size estimation such as "Freedman.Diaconis", "Sqrt", "Sturges", "Rice", "Doane", and "Scott.Normal".
#'
#' @param time_vector A numeric vector. The time series data points that need to be binned.
#' @param nPoints An integer. The total number of data points in time_vector.
#' @param drop_fac A numeric. A factor to adjust the calculated bin size. The estimated bin size is multiplied by this value. It helps in refining the bin size when the original bin size calculation results in too many empty bins.
#' @param bin_method A character string. The bin_method to estimate the bin size. Possible values include "Freedman.Diaconis", "Sqrt", "Sturges", "Rice", "Doane", and "Scott.Normal".
#'
#' @return
#' A numeric value representing the estimated bin size, adjusted by the drop_fac.
#'
#' @details
#' The function uses various rules for calculating the bin size:
#' - "Freedman.Diaconis": bin size is proportional to the interquartile range (IQR) and inversely proportional to the cube root of the number of data points.
#' - "Sqrt": bin size is proportional to the square root of the number of data points.
#' - "Sturges": bin size is proportional to the log (base 2) of the number of data points.
#' - "Rice": bin size is proportional to twice the cube root of the number of data points.
#' - "Doane": bin size accounts for data skewness in the calculation.
#' - "Scott.Normal": bin size is proportional to the standard deviation and inversely proportional to the cube root of the number of data points, assuming the data is nearly normal in distribution.
#'
#' After estimating the bin size, it is scaled down by a factor specified by 'drop_fac'.
#'
#' @examples
#' \dontrun{
#' estBinSize(time_vector = c(1, 2, 3, 4, 5), nPoints = 5, drop_fac = 0.5, bin_method = "Sturges")
#' }
#' @author Priyansh Srivastava \email{spriyansh29@@gmail.com}
#'
#' @keywords internal
estBinSize <- function(time_vector, nPoints, drop_fac, bin_method) {
  estBins <- switch(bin_method,
    "Freedman.Diaconis" = {
      # Freedman-Diaconis rule: bin size is proportional to the interquartile range (IQR)
      # and inversely proportional to the cube root of the number of data points.
      2 * IQR(time_vector) / nPoints^(1 / 3)
    },
    "Sqrt" = {
      # Square root rule: bin size is proportional to the square root of the number of data points.
      nPoints^(1 / 2)
    },
    "Sturges" = {
      # Sturges' rule: bin size is proportional to the log (base 2) of the number of data points.
      log2(nPoints) + 1
    },
    "Rice" = {
      # Rice Rule: bin size is proportional to twice the cube root of the number of data points.
      2 * nPoints^(1 / 3)
    },
    "Doane" = {
      # Doane's formula: accounts for data skewness in the calculation of bin size.
      sigma <- ((6 * (nPoints - 2)) / ((nPoints + 1) * (nPoints + 3)))^(1 / 2)
      sk <- skewness(time_vector)
      1 + log2(nPoints) + log2(1 + (abs(sk) / sigma))
    },
    "Scott.Normal" = {
      # Scott's normal reference rule: assumes the data is nearly normal in distribution.
      # Bin size is proportional to the standard deviation and inversely proportional to the cube root of the number of data points.
      3.49 * abs(sd(time_vector)) / nPoints^(1 / 3)
    },
    stop(paste("Invalid bin_method: ", bin_method, ". Please choose one of the following: 'Freedman.Diaconis', 'Sqrt', 'Sturges', 'Rice', 'Doane', 'Scott.Normal'"))
  )


  # Scale the estimated bin size by the drop factor.
  estBins <- round(drop_fac * estBins, 4)

  return(estBins)
}
