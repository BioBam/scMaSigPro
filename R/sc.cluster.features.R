#' @export
sc.cluster.features <- function(scmpObj, includeInflu = FALSE,
                                cluster.method = "hclust",
                                distance = "cor",
                                k = 9, k.mclust=FALSE,
                                cor.methods_covariance = "pairwise.complete.obs",
                                cor.method = "pearson",
                                hclust.agglo_method = "ward.D"){
    
    # Check Object Validity
    assert_that(is(scmpObj, "scMaSigProClass"),
                msg = "Please provide object of class 'scMaSigPro'"
    )
    
    # Check if compressed exist
    assert_that(
        all(length(scmpObj@compress.sce@colData$Time) > 1 &
                length(scmpObj@compress.sce@colData$Group) >=2 &
                nrow(scmpObj@compress.sce@assays@data@listData$bulk.counts > 1)),
        msg = "Please run 'squeeze', sc.pvector' & 'sc.Tfit' before"
    )
    
    # Extract Data
    sigCounts <- showSigProf(scmpObj = scmpObj, includeInflu = includeInflu)
    sigCounts[is.na(sigCounts)] <- 0
    sigCoeff <- showCoeff(scmpObj = scmpObj, includeInflu = includeInflu)
    sigCoeff[is.na(sigCoeff)] <- 0
    
    # Create list for 
    matrix.list <- list(sigCounts = sigCounts,
                        sigCoeff = sigCoeff)
    
    # Run Clustering
    cluster.list <- lapply(matrix.list, function(sig.element, cluster_method = cluster.method,
                                                 dis = distance, hclust_agglo_method = hclust.agglo_method,
                                                 cor_method = cor.method,
                                                 numClus = k, cor_methods_covariance = cor.methods_covariance){
        
        
        # Hclust
        if(cluster_method == "hclust"){
            
            # Calculate correlation-based dissimilarity matrix
            dcorrel <- matrix(rep(1, nrow(sig.element)^2), nrow(sig.element)) - cor(t(sig.element), use = cor_methods_covariance,
                                                                                    method = cor_method)
            # Run hclust
            clustering.result <- hclust(as.dist(dcorrel), method = hclust_agglo_method)
            
            # Extract Vector
            cluster.vector <- as.vector(cutree(clustering.result, numClus))
            names(cluster.vector) <- rownames(sig.element)
        }
        else if (cluster_method == "kmeans"){
            
            # Standardizing the data can be important for k-means
            standardized_data <- scale(sig.element) 
            
            # Run k-means clustering
            clustering.result <- kmeans(standardized_data, centers = numClus, nstart = 25)
            # nstart parameter is used to set the number of random sets chosen
            
            # Extract Vector
            cluster.vector <- as.vector(clustering.result$cluster)
            names(cluster.vector) <- rownames(sig.element)
        }else if (cluster_method == "mclust"){
            
            # Run Mclust
            # Mclust automatically determines the optimal number of clusters
            # You can specify the number of clusters if desired using the G parameter
            clustering.result <- Mclust(sig.element, G = numClus)
            
            # Extract the best model and get cluster assignments
            bestModel <- summary(clustering.result, parameters = TRUE)
            cluster.vector <- as.vector(as.integer(bestModel$classification))
            
            # Optionally, you can extract the BIC values and other model details
            # bic_values <- clustering.result$BIC
            
            # Assign row names
            names(cluster.vector) <- rownames(sig.element)
        }
        return(cluster.vector)
    })
    
    # Add information to the compressed Data
    scmpObj@sig.genes@feature.clusters <- list(
        sigCounts = cluster.list[["sigCounts"]],
        sigCoeff = cluster.list[["sigCoeff"]]
    )
    
    # return Object
    return(scmpObj)
}