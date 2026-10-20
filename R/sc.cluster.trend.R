#' @title Cluster the counts or coefficients.
#'
#' @description
#' This function clusters the counts or coefficients to visualize the collective
#' trends later.
#'
#' @importFrom stats cutree hclust
#' @importFrom stats as.dist cor kmeans
#'
#' @param scmpObj An object of class \code{\link{ScMaSigPro}}.
#' @param cluster_by Whether to use counts or coefficients for clustering.
#' @param geneSet Specify the gene set to be used for clustering.
#' (Default is "intersect")
#' @param cluster_method Clustering method for data partioning. Currently "hclust",
#' "kmeans" and "Mclust" are supported. "Mclust" needs the suggested package
#' 'mclust'.
#' @param hclust.agglo.method Aggregation Method. (Default is "ward.D")
#' @param hclust.distance Distance measurement.(Default is "cor")
#' @param k Number of clusters for data partioning.
#' @param kmeans.iter.max Maximum number of iterations when cluster.method is
#' kmeans
#' @param mclust.k TRUE for computing the optimal number of clusters with
#' Mclust algorithm. (Default is FALSE)
#' @param fill_na Fill the NAs (if Present), by zero, mean or median.
#' @param use_dim Whether to use rows or columns for taking mean or median while
#' filling the NAs with `fill_na`.
#' @param includeInflu Whether to include genes with influential observations.
#' @param verbose Print detailed output in the console. (Default is TRUE)
#'
#' @return An object of class \code{\link{ScMaSigPro}}, with updated
#' `Significant` slot with clusters.
#'
#' @seealso `maSigPro::see.genes()`
#'
#' @references{Conesa, A., Nueda M.J., Alberto Ferrer, A., Talon, T. 2006.
#' maSigPro: a Method to Identify Significant Differential Expression Profiles
#' in Time-Course Microarray Experiments. Bioinformatics 22, 1096-1102}
#'
#' @author Priyansh Srivastava \email{spriyansh29@@gmail.com}, Ana Conesa and
#' Maria Jose Nueda, \email{mj.nueda@@ua.es}
#'
#' @export
sc.cluster.trend <- function(scmpObj,
                             geneSet = "intersect",
                             cluster_by = "coeff",
                             cluster_method = "hclust",
                             hclust.agglo.method = "ward.D",
                             kmeans.iter.max = 500,
                             hclust.distance = "cor",
                             includeInflu = TRUE,
                             mclust.k = FALSE,
                             k = 9,
                             verbose = FALSE,
                             fill_na = "zero",
                             use_dim = "col") {
    # Global vars
    scmp_clusters <- "scmp_clusters"

    # Check Object Validity
    .check_scmp(scmpObj)

    # The gene sets come from sc.filter()
    .check(!isEmpty(scmpObj@Significant@genes),
        msg = "'Significant' slot is empty, please run 'sc.filter()' first."
    )

    # Check if the gene set exists
    geneSet <- .check_choice(geneSet,
        c(names(scmpObj@Significant@genes), "intersect", "union"),
        msg = paste(
            paste0("'", geneSet, "'"), "does not exist. Please use one of",
            paste(
                c(names(scmpObj@Significant@genes), "intersect", "union"),
                collapse = ", "
            )
        )
    )
    cluster_by <- .check_choice(cluster_by, c("coeff", "counts"),
        msg = paste(
            paste0("'", cluster_by, "'"),
            "is not a valid option. Please use one of",
            paste(c("coeff", "counts"), collapse = ", ")
        )
    )
    use_dim <- .check_choice(use_dim, c("row", "col"),
        msg = paste(
            paste0("'", use_dim, "'"),
            "is not a valid option. Please use one of",
            paste(c("row", "col"), collapse = ", ")
        )
    )
    fill_na <- .check_choice(fill_na, c("mean", "median", "zero"),
        msg = paste(
            paste0("'", fill_na, "'"),
            "is not a valid option. Please use one of",
            paste(c("mean", "median", "zero"), collapse = ", ")
        )
    )
    cluster_method <- .check_choice(cluster_method,
        c("hclust", "kmeans", "Mclust"),
        msg = paste(
            paste0("'", cluster_method, "'"),
            "is not a valid method. Please use one of",
            paste(c("hclust", "kmeans", "Mclust"), collapse = ", ")
        )
    )

    # Get gene set vector
    if (geneSet == "intersect") {
        gene_set_vector <- Reduce(intersect, scmpObj@Significant@genes)
    } else if (geneSet == "union") {
        gene_set_vector <- Reduce(union, scmpObj@Significant@genes)
    } else {
        gene_set_vector <- scmpObj@Significant@genes[[geneSet]]
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

    # Enough genes for the clusters: hclust() needs two genes and at most one
    # cluster per gene, kmeans() and Mclust() more genes than clusters (with
    # one gene, the NA filling below would return a vector)
    n_genes <- nrow(cluster_matrix_input)
    uses_k <- cluster_method != "Mclust" || !mclust.k
    whole_k <- is.numeric(k) && length(k) == 1L && !is.na(k) && k == round(k)
    min_genes <- 2
    if (uses_k && whole_k) {
        min_genes <- max(2, if (cluster_method == "hclust") k else k + 1)
    }
    .check(n_genes >= min_genes,
        msg = paste0(
            "The '", geneSet, "' gene set has ", n_genes, " gene(s), too few ",
            "for cluster_method = '", cluster_method, "'",
            if (uses_k) paste0(" with k = ", k), ": it needs at least ",
            min_genes, "."
        )
    )

    # Check fill by
    if (use_dim == "row") {
        dim <- 1
    } else if (
        use_dim == "col"
    ) {
        dim <- 2
    }

    # Replace NAs with column means
    if (fill_na == "mean") {
        cluster_matrix_input <- apply(cluster_matrix_input, dim, function(x) {
            x[is.na(x)] <- mean(x, na.rm = TRUE)
            return(x)
        })
    } else if (fill_na == "median") {
        cluster_matrix_input <- apply(cluster_matrix_input, dim, function(x) {
            x[is.na(x)] <- median(x, na.rm = TRUE)
            return(x)
        })
    } else if (fill_na == "zero") {
        cluster_matrix_input <- apply(cluster_matrix_input, dim, function(x) {
            x[is.na(x)] <- 0
            return(x)
        })
    }

    # Transpose if neccesary
    if (use_dim == "col") {
        # Transpose the matrix
        cluster_data <- cluster_matrix_input
    } else if (use_dim == "row") {
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
            clust <- hclust(dist(cluster_data, method = hclust.distance),
                method = hclust.agglo.method
            )
        }
        # Cut the dendrogram to get cluster assignments
        cluster_results <- cutree(clust, k = k)
    } else if (cluster_method == "kmeans") {
        kmeans_result <- kmeans(cluster_data, centers = k, iter.max = kmeans.iter.max)
        cluster_results <- kmeans_result$cluster
    } else if (cluster_method == "Mclust") {
        .check(requireNamespace("mclust", quietly = TRUE),
            msg = paste(
                "cluster_method = 'Mclust' needs the 'mclust' package;",
                "install it with install.packages('mclust')."
            )
        )
        # Mclust() evaluates mclustBIC() in the frame it is called from, so
        # call it from a frame enclosed by the mclust namespace: then it
        # works without library(mclust)
        mclust_frame <- new.env(parent = asNamespace("mclust"))
        mclust_frame$cluster_data <- cluster_data
        mclust_frame$k <- k
        if (mclust.k) {
            mclust_res <- eval(
                quote(Mclust(cluster_data, verbose = FALSE)), mclust_frame
            )
            k <- mclust_res$G
        } else {
            mclust_res <- eval(
                quote(Mclust(cluster_data, G = k, verbose = FALSE)),
                mclust_frame
            )
        }
        cluster_results <- mclust_res$classification
    }

    # One cluster per gene, named after it: with one component, Mclust()
    # returns its classification without names
    names(cluster_results) <- rownames(cluster_data)

    # Convert to df
    clusters_df <- as.data.frame(cluster_results)
    colnames(clusters_df) <- scmp_clusters
    clusters_df[["feature_id"]] <- rownames(clusters_df)

    cluster.vector <- clusters_df[["scmp_clusters"]]
    names(cluster.vector) <- clusters_df[["feature_id"]]
    scmpObj@Significant@clusters <- as.list(cluster.vector)

    # Update Parameters
    scmpObj@Parameters@cluster_method <- cluster_method
    scmpObj@Parameters@use_dim <- use_dim
    scmpObj@Parameters@fill_na <- fill_na
    validObject(scmpObj)

    return(scmpObj)
}
