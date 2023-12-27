#################################################
## Author: Priyansh Srivastava ##################
## Year: 2023 ###################################
## Email: spriyansh29@gmail.com #################
#################################################
# Internal function for pseudo-bulking counts along Pseudotime. They are called
# by sc.squeeze()
###############################################################################
#' @title Calculate Bin Size Function
#'
#' @description
#' `calc_bin_size()` calculates the bin size based on the number of elements in
#' the "cluster.members" column of the input data frame.
#'
#' @param x A data frame containing the "cluster.members" column.
#'
#' @return A numeric value representing the size of the bin (number of elements
#' in the "cluster.members" column).
#'
#' @author Priyansh Srivastava \email{spriyansh29@@gmail.com}
#'
#' @importFrom stringr str_split
#'
#' @keywords internal
# Define a function 'calc_bin_size' which takes a data frame 'x' as input

calc_bin_size <- function(x, clus_mem_col = "scmp_cluster_members") {
  # Use the 'str_split' function from the 'stringr' package to split the 'cluster.members' column
  # of the input data frame 'x' by the '|' character.
  # This returns a list where each element is a vector of the split strings.
  # 'c()' is used to concatenate these vectors into a single vector.
  # Finally, 'length' is used to get the length of this vector (i.e., the number of split strings),
  # which is stored in the 'size' variable.
  size <- length(c(str_split(x[[clus_mem_col]], "\\|"))[[1]])

  # Convert the 'size' variable to a numeric value and return it as the result of the function
  return(as.numeric(size))
}
###############################################################################
#' @title Convert Vector Elements to Path Names (Internal)
#'
#' @description
#' `convert_to_path()` transforms a vector by renaming its unique elements
#' (excluding the "root" element) to a sequence named "Path1", "Path2", etc.
#'
#' @param vec A character vector where elements may be repeated and
#' might contain the value "root".
#' @param path_prefix Prefix used to annoate the paths, default is "Path".
#'
#' @return A character vector with the same length as the input where
#' unique elements, excluding "root", are renamed to "Path1", "Path2", etc.
#'
#' @importFrom stats setNames
#'
#' @keywords internal
convert_to_path <- function(vec, path_prefix, root_label) {
  # Exclude "root" from the transformation and get the unique values
  unique_vals <- unique(vec[vec != root_label])

  # Create a named vector for the mapping
  name_map <- setNames(paste0(path_prefix, 1:length(unique_vals)), unique_vals)

  # Map and replace the non-root elements
  vec[vec != root_label] <- name_map[vec[vec != root_label]]

  return(vec)
}
###############################################################################
#' @title Create Range Function
#'
#' @description
#' `create_range()` converts a factor column "bin" into a character vector, extracts
#' numeric range values from the character vector, and combines them with
#' additional columns "bin_size" and "binned_time" to return a numeric vector
#' representing the range information.
#'
#' @param x A data frame that should contain the following columns:
#'   \describe{
#'     \item{bin}{A factor column representing the bin intervals in the format "[x, y]".}
#'     \item{bin_size}{A numeric column representing the bin size.}
#'     \item{binned_time}{A numeric column representing the binned time.}
#'   }
#'
#' @return A numeric vector containing four elements:
#'   \describe{
#'     \item{lower_bound}{The lower bound of the bin interval.}
#'     \item{upper_bound}{The upper bound of the bin interval.}
#'     \item{bin_size}{The bin size.}
#'     \item{binned_time}{The binned time.}
#'   }
#'
#' @author Priyansh Srivastava \email{spriyansh29@@gmail.com}
#'
#' @importFrom stringr str_remove_all
#'
#' @keywords internal
create_range <- function(x, bin_size_colname = "scmp_bin_size",
                         bin_col = "scmp_bin", verbose = TRUE) {
  # Convert the factor column "bin" to character
  y <- as.character(x[[bin_col]])

  # Remove square and round brackets from the character string
  y <- y %>% str_remove_all(pattern = "\\[|\\]|\\(|\\)")

  # Split the character string by comma and extract the first element (lower bound of the range)
  y1 <- as.numeric(sapply(strsplit(y, ","), "[", 1))
  # Extract the first element
  # y1 <- as.numeric(vapply(strsplit(y, ","), function(x) x[1], numeric(1)))

  # Split the character string by comma and extract the second element (upper bound of the range)
  y2 <- as.numeric(sapply(strsplit(y, ","), "[", 2))
  # y2 <- as.numeric(vapply(strsplit(y, ","), function(x) x[2], numeric(1)))

  if (verbose) {
    message(paste0(
      "Lower Bound:", y1, ", Upper Bound:", y2, ", Number of cells:", x[[bin_size_colname]]
    ))
  }

  # Combine the lower bound, upper bound, bin size, and binned time into a numeric vector
  rangeVec <- c(y1, y2, x[[bin_size_colname]])

  # Return the numeric vector
  return(as.numeric(rangeVec))
}
###############################################################################
#' @title estBinSize
#'
#' @description
#' `estBinSize()` calculates the optimal bin size for discretizing a continuous variable.
#' It uses various bin_methods for bin size estimation such as "Freedman.Diaconis",
#' "Sqrt", "Sturges", "Rice", "Doane", and "Scott.Normal".
#'
#' @param time_vector A numeric vector. The time series data points that need to be binned.
#' @param nPoints An integer. The total number of data points in time_vector.
#' @param drop_fac A numeric. A factor to adjust the calculated bin size. The
#' estimated bin size is multiplied by this value. It helps in refining the bin
#' size when the original bin size calculation results in too many empty bins.
#' @param bin_method A character string. The bin_method to estimate the bin size.
#' Possible values include "Freedman.Diaconis", "Sqrt", "Sturges", "Rice", "Doane",
#' and "Scott.Normal". See details.
#'
#' @return
#' A numeric value representing the estimated bin size, adjusted by the drop_fac.
#'
#' @details
#' The function contains various rules for calculating the bin size:
#' \describe{
#'   \item{"Freedman.Diaconis"}{bin size is proportional to the interquartile range
#'   (IQR) and inversely proportional to the cube root of the number of data points.}
#'   \item{"Sqrt"}{bin size is proportional to the square root of the number of data points.}
#'   \item{"Sturges"}{bin size is proportional to the log (base 2) of the number of data points.}
#'   \item{"Rice"}{bin size is proportional to twice the cube root of the number of data points.}
#'   \item{"Doane"}{bin size accounts for data skewness in the calculation.}
#'   \item{"Scott.Normal"}{bin size is proportional to the standard deviation and
#'   inversely proportional to the cube root of the number of data points, assuming
#'   the data is nearly normal in distribution.}
#' }
#' After estimating the bin size, it is scaled down by a factor specified by 'drop_fac'.
#'
#' @author Priyansh Srivastava \email{spriyansh29@@gmail.com}
#'
#' @importFrom stats IQR sd
#' @importFrom e1071 skewness
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
  estBins <- drop_fac * estBins

  return(estBins)
}
###############################################################################
#' @title Extract Fitting
#'
#' @param reg Regression Result form `stats::glm()`
#' @param lmf Full model object generated from the `stats::glm()`
#' @param model.glm.0 Intercept model object generated from the `stats::glm()`
#' @param dis Dataframe for the predictors
#' @param family Family used for modelling the data.
#' @param name description
#' @param vars.in description
#' @param alfa Significance threshold (Q-Value)
#' @param influ.info description
#'
#' @keywords internal
extract_fitting <- function(reg, lmf, model.glm.0, dis, family, name, vars.in, alfa, influ.info) {
  sol <- coefficients <- group.coeffs <- t.score <- sig.profiles <- NULL
  y <- reg$y
  result <- summary(lmf)
  novar <- vars.in[!is.element(vars.in, names(result$coefficients[, 4]))]
  influ <- influence.measures(reg)$is.inf
  influ <- influ[, c(ncol(influ) - 3, ncol(influ) - 1)]
  influ1 <- which(apply(influ, 1, all))
  if (length(influ1) != 0) {
    paste.names <- function(a) {
      paste(names(a)[a], collapse = "/")
    }
    match <- match(rownames(dis), rownames(influ))
    influ <- as.data.frame(apply(influ, 1, paste.names))
    influ.info <- cbind(influ.info, influ[match, ])
    colnames(influ.info)[ncol(influ.info)] <- name
    influ.info <- as.matrix(influ.info)
  }
  result <- summary(reg)
  if ((!(result$aic == -Inf) & !is.na(result$aic) & family$family == "gaussian") | family$family != "gaussian") {
    if (family$family == "gaussian") {
      test <- anova(model.glm.0, reg, test = "F")
      p.value <- test[6][2, 1]
    } else {
      test <- anova(model.glm.0, reg, test = "Chisq")
      p.value <- test[5][2, 1]
    }
    # Computing goodness of fitting:
    bondad <- (reg$null.deviance - reg$deviance) / reg$null.deviance
    if (bondad < 0) {
      bondad <- 0
    }
    beta.coeff <- result$coefficients[, 1]
    beta.p.valor <- result$coefficients[, 4]
    coeff <- rep(0, (length(vars.in) + 1))
    if (length(novar) != 0) {
      for (m in 1:length(novar)) {
        coeff[position(dis, novar[m]) + 1] <- NA
      }
    }
    p.valor <- t <- as.numeric(rep(NA, (length(vars.in) + 1)))

    if (result$coefficients[, 4][rownames(result$coefficients) ==
      "(Intercept)"] < alfa) {
      coeff[1] <- result$coefficients[, 1][rownames(result$coefficients) ==
        "(Intercept)"]
      p.valor[1] <- result$coefficients[, 4][rownames(result$coefficients) ==
        "(Intercept)"]
      t[1] <- result$coefficients[, 3][rownames(result$coefficients) ==
        "(Intercept)"]
    }
    for (j in 2:length(coeff)) {
      if (is.element(vars.in[j - 1], rownames(result$coefficients))) {
        coeff[j] <- result$coefficients[, 1][rownames(result$coefficients) ==
          vars.in[j - 1]]
        p.valor[j] <- result$coefficients[, 4][rownames(result$coefficients) ==
          vars.in[j - 1]]
        t[j] <- result$coefficients[, 3][rownames(result$coefficients) ==
          vars.in[j - 1]]
      }
    }
    if (!all(is.na(p.valor))) {
      sol <- rbind(sol, as.numeric(c(
        p.value, bondad,
        p.valor
      )))
      coefficients <- rbind(coefficients, coeff)
      t.score <- rbind(t.score, t)
      sig.profiles <- rbind(sig.profiles, y)
      h <- nrow(sol)
      rownames(sol)[h] <- name
      rownames(coefficients)[h] <- name
      rownames(t.score)[h] <- name
      rownames(sig.profiles)[h] <- name
    }
  }
  # Return Calculation
  return(list(
    p_value = p.value,
    bondad = bondad,
    p_valor = p.valor,
    coeff = coeff,
    t = t,
    sig_profiles = y,
    sol = sol,
    influ.info = influ.info,
    feature_name = name
  ))
}

###############################################################################
#' @title Extract Intervals
#'
#' @param time.vector Numeric vector of Pseudotime values.
#' @param nBins Expected number of bins.
#' @param bin Column name for the bin column.
#' @param bin.size Column name for the bin size column.
#' @param lbond Column name for the lower bound column.
#' @param ubond Column name for the upper bound column.
#' @keywords internal

extract_interval <- function(time.vector, nBins = 1, bin, bin.size, lbound, ubound) {
  # Create Dataframe
  new_range_current <- as.data.frame(entropy::discretize(
    time.vector,
    numBins = nBins, r = range(time.vector)
  ))
  # Set columns
  colnames(new_range_current) <- c(bin, bin.size)

  # Current bin new table
  new_bin_table_current <- as.data.frame(
    t(as.data.frame(
      apply(
        new_range_current, 1, create_range,
        bin_size_colname = bin.size,
        bin_col = bin,
        verbose = FALSE
      )
    ))
  )
  # Set column names
  colnames(new_bin_table_current) <- c(lbound, ubound, bin.size)

  # Return
  return(new_bin_table_current)
}

###############################################################################
#' Select the Longer of Two Vectors
#'
#' This function compares the lengths of two vectors and returns the longer one.
#' If both vectors have the same length, the function returns 1.
#'
#' @param vector1 A numeric vector.
#' @param vector2 A numeric vector.
#' @param vector1_label A label (character string) for `vector1`.
#' @param vector2_label A label (character string) for `vector2`.
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{long_vec}: The longer vector of `vector1` and `vector2`.
#'   \item \code{long_vec_label}: The label corresponding to the longer vector.
#'   \item \code{short_vec}: The shorter vector of `vector1` and `vector2`.
#'   \item \code{short_vec_label}: The label corresponding to the shorter vector.
#' }
#' If both vectors are of the same length, the function returns 1.
#'
#' @author Priyansh Srivastava \email{spriyansh29@@gmail.com}
#' @keywords internal
select_longer_vector <- function(vector1, vector2,
                                 vector1_label, vector2_label) {
  if (length(vector1) > length(vector2)) {
    return(list(
      long_vec = vector1, long_vec_label = vector1_label,
      short_vec = vector2, short_vec_label = vector2_label
    ))
  } else if (length(vector2) > length(vector1)) {
    return(list(
      long_vec = vector2, long_vec_label = vector2_label,
      short_vec = vector1, short_vec_label = vector1_label
    ))
  } else {
    return(list(empty = 1))
  }
}

###############################################################################
#' @title Homogenize Bins of a Dataframe
#'
#' @description
#' Adjust the bins of a dataframe based on a specified maximum bin size threshold.
#' It checks if the bin size is too large relative to this threshold.
#' If a bin is too large, it will be split into two equal parts until all bins are
#' below the maximum allowed size. Optionally, bins too small can be merged with
#' adjacent bins based on the provided method.
#'
#' @param bin_table A dataframe containing the bin data.
#' @param max_allowed The maximum allowed size for a bin.
#' @param verbose Logical; if TRUE, prints detailed output.
#' @param time_vector A vector representing time or a sequence that bins refer to.
#' @param lbound The name of the lower bound column in `bin_table`.
#' @param ubound The name of the upper bound column in `bin_table`.
#' @param bin The name of the bin identifier column in `bin_table`.
#' @param bin_size The name of the bin size column in `bin_table`.
#' @param method The method for handling small bins: 'merge' to merge with previous or next bin,
#' 'drop' to remove small bins, or 'ignore' to leave small bins as they are.
#' @param drop The threshold below which a bin is considered too small and subject to the method.
#'
#' @return A dataframe with adjusted bins. The structure of the dataframe will be the same as `bin_table`
#' with updated values for the bin size and bounds.
#' @keywords internal
optimize_bin_max <- function(bin_table, max_allowed, verbose = TRUE,
                             time_vector, lbound, ubound, bin, bin.size, method,
                             drop) {
  # Initiate an empty uniform bin
  uniform_bin_df <- data.frame(matrix(NA, ncol = ncol(bin_table)))
  colnames(uniform_bin_df) <- colnames(bin_table)

  # Binning one-by one
  for (i in c(1:nrow(bin_table))) {
    # histo_bin
    current_bin_df <- c(unlist(bin_table[i, , drop = FALSE]))

    # Get size of the bins
    current_bin_size <- current_bin_df[[bin.size]]

    # Check if the bin size is big
    if (current_bin_size > max_allowed) {
      # get pseudotime
      pTime.for.big.interval <- time_vector[(time_vector >= current_bin_df[[lbound]] & time_vector <= current_bin_df[[ubound]])]

      # Estimate how many splits are required
      potential_splits <- ceiling(estBinSize(
        time_vector = pTime.for.big.interval, nPoints = length(pTime.for.big.interval),
        drop_fac = drop, bin_method = method
      ))

      # Run the extraction
      small_splitted_bin_df <- extract_interval(
        time.vector = pTime.for.big.interval,
        nBins = potential_splits,
        bin = bin, bin.size = bin.size, lbound = lbound, ubound = ubound
      )
      # Validation
      if (verbose) {
        message(paste("Splitting bin with", current_bin_size, "cells into", nrow(small_splitted_bin_df), "with sizes as", paste(small_splitted_bin_df[[bin.size]], collapse = ", ")))
      }

      # Add new Row
      uniform_bin_df <- rbind(uniform_bin_df, small_splitted_bin_df)
    } else if (current_bin_size <= max_allowed) {
      uniform_bin_df <- rbind(uniform_bin_df, current_bin_df)
      if (verbose) {
        if (current_bin_size != 0) {
          message(paste("Skipping bin with size", current_bin_size))
        }
      }
    }
  }
  uniform_bin_df <- uniform_bin_df[-1, ]

  # Correct the rows
  rownames(uniform_bin_df) <- NULL

  # Return the new data
  return(uniform_bin_df)
}

################################################################################
#' @title Create Pseduo-bulk Counts
#'
#' @description
#' `pb_counts()` creates a dataframe of pseudo bulk counts from single
#' cell counts. It does this by either taking the mean or sum of counts across
#' cells in each bin, depending on the specified method.
#'
#' @param scmpObj object of Class scMaSigPro. See \code{\link{ScMaSigPro}}
#' for more details.
#' @param bin_mem_col Column name in the Dense metadata storing information
#' about the members of the bins. (Default is 'scmp_bin_members').
#' @param bin_col  Column name in the Dense metadata storing information
#' about the bin labels. (Default is 'scmp_bin').
#' @param cluster_count_by A character string specifying the method to use to
#' aggregate counts within each cluster. Available options are 'mean' or 'sum'.
#' (Default = "sum").
#' @param assay_name Name of the Assay in the assay_name object from which
#' retrieve the counts. (Default = "counts").
#'
#' @return
#' A matrix. The matrix includes pseudo bulk counts with each row being
#' a gene and each column being a bin.
#'
#' @details
#' The function operates by iterating over each row of the pseudo_bulk_profile.
#' For each bin, it identifies the cells that belong to the bin and selects their
#' counts from the counts data frame. It then calculates the mean or sum of these
#' counts (depending on the specified method), and adds these to a new data frame
#' of pseudo bulk counts. The result is a pseudo bulk counts data frame where
#' each row is a gene and each column is a bin.
#'
#' @author Priyansh Srivastava \email{spriyansh29@@gmail.com}
#'
#' @keywords internal

pb_counts <- function(scmpObj,
                      bin_mem_col = scmpObj@Parameters@bin_mem_col,
                      bin_col = scmpObj@Parameters@bin_col,
                      assay_name = "counts",
                      cluster_count_by = "sum") {
  # Check Object Validity
  assert_that(is(scmpObj, "ScMaSigPro"),
    msg = "Please provide object of class 'scMaSigPro'."
  )

  # Count slot
  assert_that(
    all(
      assay_name %in% names(scmpObj@Sparse@assays@data@listData)
    ),
    msg = paste0("'", assay_name, "' ", "doesn't exit in scmpObj.")
  )

  # Get assay
  counts <- scmpObj@Sparse@assays@data@listData[[assay_name]]

  # Get Pseudobulk Profile
  pseudo_bulk_profile <- as.data.frame(colData(scmpObj@Dense))

  assert_that(bin_mem_col %in% colnames(pseudo_bulk_profile),
    msg = paste0("'", bin_mem_col, "' does not exist in level.meta.data")
  )
  assert_that(bin_col %in% colnames(pseudo_bulk_profile),
    msg = paste0("'", bin_col, "' does not exist in level.meta.data")
  )

  # Get the meta-information for pseudobulking
  meta.info <- pseudo_bulk_profile[, c(bin_mem_col, bin_col)]

  # Run mclapply
  pb.counts <- lapply(1:nrow(meta.info), function(i) {
    # Get the bin.info
    bin <- meta.info[i, , drop = FALSE]

    # Split the row
    cell.vector <- c(str_split(bin[1], "\\|"))[[1]]

    # Get col cells
    col_indices <- which(colnames(counts) %in% cell.vector)

    # Subset the matrix using these indices
    bin_matrix <- as.matrix(counts[, col_indices, drop = FALSE])

    # Get Pseudobulked-counts
    pb.vector <- switch(cluster_count_by,
      "mean" = as.matrix(round(rowMeans(bin_matrix))),
      "sum"  = as.matrix(rowSums(bin_matrix)),
      stop("Invalid cluster_count_by value. Please choose either 'mean' or 'sum'.")
    )

    # Return
    return(pb.vector)
  })

  # Convert the list output of mclapply to a matrix and set the row names
  pb.counts <- do.call(cbind, pb.counts)
  rownames(pb.counts) <- rownames(counts)
  colnames(pb.counts) <- meta.info[[bin_col]]

  # Return the counts
  scmpObj@Dense@assays@data@listData$bulk.counts <- as(pb.counts, "dgCMatrix")

  # return
  return(scmpObj)
}
