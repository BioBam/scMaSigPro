#' Cluster the counts or coefficients
#'
#' This function clusters the counts or coefficients
#'
#' @param scmpObj object of class scmpObj
#' @param cluster_by description
#' @param geneSet description
#' @param cluster_method description
#' @param hclust.agglo.method description
#' @param hclust.distance description
#' @param includeInflu description
#' @param k description
#' @param verbose description
#' @param kmeans.iter.max description
#' @param mclust.k description
#' @param fill_na_by description
#' @param fill_dim description
#'
#' @import ggplot2
#' @importFrom stats complete.cases cutree hclust
#' @importFrom RColorConesa getConesaColors
#' @importFrom mclust Mclust
#' @importFrom stats as.dist cor kmeans
#' @return Generates a plot.
#' @export
sc.cluster.trend <- function(scmpObj,
                             geneSet,
                             cluster_by = "coeff",
                             cluster_method = "hclust",
                             hclust.agglo.method = "ward.D",
                             kmeans.iter.max = 500,
                             hclust.distance = "cor",
                             includeInflu = TRUE,
                             mclust.k = FALSE,
                             k = 9,
                             verbose = FALSE,
                             fill_na_by = "zero",
                             fill_dim = "col") {
  # Global vars
  scmp_clusters <- "scmp_clusters"

  # Check if the gene set exists
  assert_that(any(geneSet %in% c(names(scmpObj@sig.genes@sig.genes), "intersect", "union")),
    msg = paste(
      paste0("'", geneSet, "'"), "does not exist. Please use one of",
      paste(c(names(scmpObj@sig.genes@sig.genes), "intersect", "union"), collapse = ", ")
    )
  )
  assert_that(any(cluster_by %in% c("coeff", "counts")),
    msg = paste(
      paste0("'", cluster_by, "'"), "is not a valid option. Please use one of",
      paste(c("coeff", "counts"), collapse = ", ")
    )
  )
  assert_that(any(fill_dim %in% c("row", "col")),
    msg = paste(
      paste0("'", fill_dim, "'"), "is not a valid option. Please use one of",
      paste(c("row", "col"), collapse = ", ")
    )
  )
  assert_that(any(fill_na_by %in% c("mean", "median", "zero")),
    msg = paste(
      paste0("'", fill_na_by, "'"), "is not a valid option. Please use one of",
      paste(c("mean", "median", "zero"), collapse = ", ")
    )
  )
  assert_that(any(cluster_method %in% c("hclust", "kmeans", "Mclust")),
    msg = paste(
      paste0("'", cluster_method, "'"), "is not a valid method. Please use one of",
      paste(c("hclust", "kmeans", "Mclust"), collapse = ", ")
    )
  )

  # Get gene set vector
  if (geneSet == "intersect") {
    gene_set_vector <- Reduce(intersect, scmpObj@sig.genes@sig.genes)
  } else if (geneSet == "union") {
    gene_set_vector <- Reduce(union, scmpObj@sig.genes@sig.genes)
  } else {
    gene_set_vector <- scmpObj@sig.genes@sig.genes[[geneSet]]
  }

  # Extract data based on 'cluster_by'
  if (cluster_by == "counts") {
    # Extract bulk counts
    cluster_matrix_input <- as.matrix(showSigProf(scmpObj, includeInflu = includeInflu))
    cluster_matrix_input <- cluster_matrix_input[rownames(cluster_matrix_input) %in% gene_set_vector, , drop = FALSE]
  } else if (cluster_by == "coeff") {
    # Extract coefficients
    cluster_matrix_input <- as.matrix(showCoeff(scmpObj, includeInflu = includeInflu))
    cluster_matrix_input <- cluster_matrix_input[rownames(cluster_matrix_input) %in% gene_set_vector, , drop = FALSE]
  }

  # Check fill by
  if (fill_dim == "row") {
    use_dim <- 1
  } else if (
    fill_dim == "col"
  ) {
    use_dim <- 2
  }

  # Replace NAs with column means
  if (fill_na_by == "mean") {
    cluster_matrix_input <- apply(cluster_matrix_input, use_dim, function(x) {
      x[is.na(x)] <- mean(x, na.rm = TRUE)
      return(x)
    })
  } else if (fill_na_by == "median") {
    cluster_matrix_input <- apply(cluster_matrix_input, use_dim, function(x) {
      x[is.na(x)] <- median(x, na.rm = TRUE)
      return(x)
    })
  } else if (fill_na_by == "zero") {
    cluster_matrix_input <- apply(cluster_matrix_input, use_dim, function(x) {
      x[is.na(x)] <- 0
      return(x)
    })
  }

  # Transpose if neccesary
  if (fill_dim == "col") {
    # Transpose the matrix
    cluster_data <- cluster_matrix_input
  } else if(fill_dim == "row") {
    cluster_data <- t(cluster_matrix_input)
  }
  
  # Perform Clustering
  if (cluster_method == "hclust") {
      if (hclust.distance == "cor") {
          # Compute correlation-based distance for rows
          correlation_matrix <- cor(t(cluster_data), use = "pairwise.complete.obs")
          dcorrel <- as.dist(1 - correlation_matrix)
          clust <- hclust(dcorrel, method = hclust.agglo.method)
      } else {
          # Compute distance based on the specified method for rows
          clust <- hclust(dist(t(cluster_data), method = hclust.distance), method = hclust.agglo.method)
      }
      # Cut the dendrogram to get cluster assignments
      cluster_results <- cutree(clust, k = k)
  }  else if (cluster_method == "kmeans") {
      kmeans_result <- kmeans(cluster_data, centers = k, iter.max = kmeans.iter.max)
      cluster_results <- kmeans_result$cluster
  }else if (cluster_method == "Mclust") {
      if (mclust.k) {
          mclust_res <- Mclust(cluster_data,verbose = FALSE)
          k <- mclust_res$G
      } else {
          mclust_res <- Mclust(cluster_data, G = k, verbose = FALSE)
      }
      cluster_results <- mclust_res$classification
  }
  
  # Convert to df
  clusters_df <- as.data.frame(cluster_results)
  colnames(clusters_df) <- scmp_clusters
  clusters_df[["feature_id"]] <- rownames(clusters_df)

  cluster.vector <- clusters_df[["scmp_clusters"]]
  names(cluster.vector) <- clusters_df[["feature_id"]]
  scmpObj@sig.genes@feature.clusters <- as.list(cluster.vector)
  
  # Update Parameters
  scmpObj@param@cluster.method <- cluster_method
  scmpObj@param@cluster.fill.dim <- fill_dim
  scmpObj@param@cluster.fill.na <- fill_na_by
  
  return(scmpObj)
}
