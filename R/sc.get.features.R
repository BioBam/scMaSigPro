#' Get Features from scmp Object
#'
#' This function extracts features from an object of class 'scmp'
#' based on a specified query and group. The function supports two types of
#' queries: 'unique' or 'union', and allows for the selection of a specific
#' group from the available groups in the object.
#'
#' @param scmpObj An object of class 'scmp'.
#' @param query A string specifying the type of query.
#'              Valid options are 'unique' or 'union'.
#'              Default is 'unique'.
#' @param unique.group A string specifying the group to consider.
#'              Default is 'Path1'.
#' @param union.ref A string specifying the reference group for 'union' query.
#' @param union.target A string specifying the target group for 'union' query.
#' @param rsq A numeric value specifying the R-squared threshold.
#' @param Q A numeric value, derived from the scmpObj object, for analysis.
#' @param vars A string specifying the variables to be used.
#' @param significant.intercept A string indicating the treatment of intercepts in the analysis.
#' @param term.Q A numeric value specifying the threshold for term significance.
#' @param includeInflu A boolean indicating whether to include influencers in the analysis.
#' @param unique.trend A string specifying the trend type for 'unique' query.
#' @param union.ref.trend A string specifying the trend type for the reference group in 'union' query.
#' @param union.target.trend A string specifying the trend type for the target group in 'union' query.
#'
#' @return Depending on the query, the function returns either a vector of
#'         unique features (when query is 'unique') or a vector of features
#'         common to all groups (when query is 'union').
#'
#' @importFrom assertthat assert_that
#'
#' @export
#'
sc.get.features <- function(scmpObj,
                            query = "unique",
                            unique.group = "Path1",
                            union.ref = "Path1",
                            union.target = "Path2vsPath1",
                            rsq = 0.7,
                            Q = scmpObj@param@Q,
                            vars = "groups",
                            significant.intercept = "dummy",
                            term.Q = 0.05,
                            includeInflu = TRUE,
                            unique.trend = "any",
                            union.ref.trend = "any",
                            union.target.trend = "any") {
  # Check Validity of the object
  assert_that(is(scmpObj, "scmp"),
    msg = "Please provide object of class 'scmp'"
  )

  # Check
  assert_that(
    all(
      query %in% c("unique", "union")
    ),
    msg = "Invalid selction for 'vars'"
  )

  # Extract avail groups
  avail_groups <- unique(scmpObj@estimate@path)

  # Check if group exist
  assert_that(
    all(
      unique.group %in% avail_groups
    ),
    msg = paste(unique.group, "does not exist. Please choose one of", paste(avail_groups, collapse = ", "))
  )

  # Override
  scmpObj <- sc.filter(scmpObj,
    rsq = rsq,
    vars = vars,
    significant.intercept = significant.intercept,
    Q = Q,
    term.Q = term.Q,
    includeInflu = includeInflu
  )

  # Create Intersection
  interserction <- plotIntersect(scmpObj)

  # Get intersection Data
  interserction.data <- interserction[[1]][["data"]][, c("gene", avail_groups),
    drop = FALSE
  ]

  # Check
  if (query == "unique") {
    res <- lapply(interserction.data[, -1], function(col) {
      allFalse <- rowSums(interserction.data[, -1] == TRUE) == 1
      interserction.data$gene[allFalse & col]
    })

    names(res) <- colnames(interserction.data[, -1])
    res <- res[[unique.group]]

    if (unique.trend != "any") {
      # Get all coeff
      coeff.df <- showCoeff(scmpObj)

      # Subset by the genes
      coeff.df.sub <- coeff.df[res, , drop = FALSE]

      # If up trends are requested
      if (unique.trend %in% c("up", "down")) {
        # Get path related columns
        if (unique.trend == "up") {
          res <- rownames(getPattern(
            frame = coeff.df.sub, trend = "up", group = unique.group,
            groups_vector = scmpObj@estimate@path
          ))
        } else if (unique.trend == "down") {
          res <- rownames(getPattern(
            frame = coeff.df.sub, trend = "down", group = unique.group,
            groups_vector = scmpObj@estimate@path
          ))
        }
      }
    }
  } else if (query == "union") {
    res <- interserction.data$gene[rowSums(
      interserction.data[, -1] == TRUE
    ) == ncol(interserction.data[, -1])]


    if ((union.target.trend != "any" | union.ref.trend != "any")) {
      # Get all coeff
      coeff.df <- showCoeff(scmpObj)

      # Subset by the genes
      coeff.df.sub <- coeff.df[res, , drop = FALSE]

      # Both up
      if ((union.target.trend == "up" & union.ref.trend == "up")) {
        # Get list
        res.list <- lapply(avail_groups, function(group_i) {
          res.frame <- getPattern(
            frame = coeff.df.sub,
            trend = "up", group = group_i,
            groups_vector = scmpObj@estimate@path
          )
          return(rownames(res.frame))
        })

        # Get shared elements
        res <- Reduce(intersect, res.list)
      } else if ((union.target.trend == "down" & union.ref.trend == "down")) {
        # Get list
        res.list <- lapply(avail_groups, function(group_i) {
          res.frame <- getPattern(
            frame = coeff.df.sub,
            trend = "down", group = group_i,
            groups_vector = scmpObj@estimate@path
          )

          return(rownames(res.frame))
        })

        # Get shared elements
        res <- Reduce(intersect, res.list)
      } else if ((union.target.trend == "up" & union.ref.trend == "down")) {
        # Get list
        res.list <- lapply(avail_groups, function(group_i) {
          if (group_i == union.ref) {
            trend_value <- union.ref.trend
          } else if (group_i == union.target) {
            trend_value <- union.target.trend
          }

          res.frame <- getPattern(
            frame = coeff.df.sub,
            trend = trend_value, group = group_i,
            groups_vector = scmpObj@estimate@path
          )

          return(rownames(res.frame))
        })

        # Get shared elements
        res <- Reduce(intersect, res.list)
      } else if ((union.target.trend == "down" & union.ref.trend == "up")) {
        # Get list
        res.list <- lapply(avail_groups, function(group_i) {
          if (group_i == union.ref) {
            trend_value <- union.ref.trend
          } else if (group_i == union.target) {
            trend_value <- union.target.trend
          }

          res.frame <- getPattern(
            frame = coeff.df.sub,
            trend = trend_value, group = group_i,
            groups_vector = scmpObj@estimate@path
          )

          return(rownames(res.frame))
        })

        # Get shared elements
        res <- Reduce(intersect, res.list)
      } else if ((union.target.trend == "any" & union.ref.trend == "down")) {
        # Get list
        res.list <- lapply(avail_groups, function(group_i) {
          if (group_i == union.ref) {
            trend_value <- union.ref.trend
          } else if (group_i == union.target) {
            trend_value <- union.target.trend
          }

          res.frame <- getPattern(
            frame = coeff.df.sub,
            trend = trend_value, group = group_i,
            groups_vector = scmpObj@estimate@path
          )

          return(rownames(res.frame))
        })

        # Get shared elements
        res <- Reduce(intersect, res.list)
      } else if ((union.target.trend == "any" & union.ref.trend == "up")) {
        # Get list
        res.list <- lapply(avail_groups, function(group_i) {
          if (group_i == union.ref) {
            trend_value <- union.ref.trend
          } else if (group_i == union.target) {
            trend_value <- union.target.trend
          }

          res.frame <- getPattern(
            frame = coeff.df.sub,
            trend = trend_value, group = group_i,
            groups_vector = scmpObj@estimate@path
          )

          return(rownames(res.frame))
        })

        # Get shared elements
        res <- Reduce(intersect, res.list)
      } else if ((union.target.trend == "down" & union.ref.trend == "any")) {
        # Get list
        res.list <- lapply(avail_groups, function(group_i) {
          if (group_i == union.ref) {
            trend_value <- union.ref.trend
          } else if (group_i == union.target) {
            trend_value <- union.target.trend
          }

          res.frame <- getPattern(
            frame = coeff.df.sub,
            trend = trend_value, group = group_i,
            groups_vector = scmpObj@estimate@path
          )

          return(rownames(res.frame))
        })

        # Get shared elements
        res <- Reduce(intersect, res.list)
      } else if ((union.target.trend == "up" & union.ref.trend == "any")) {
        # Get list
        res.list <- lapply(avail_groups, function(group_i) {
          if (group_i == union.ref) {
            trend_value <- union.ref.trend
          } else if (group_i == union.target) {
            trend_value <- union.target.trend
          }

          res.frame <- getPattern(
            frame = coeff.df.sub,
            trend = trend_value, group = group_i,
            groups_vector = scmpObj@estimate@path
          )

          return(rownames(res.frame))
        })

        # Get shared elements
        res <- Reduce(intersect, res.list)
      }
    }
  }

  return(res)
}
