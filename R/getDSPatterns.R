getDSPatterns <- function(tableDS, Cluster.Major, Cluster.minor) {
  rownames(tableDS$IsoClusters)[tableDS$IsoClusters[, 1] == Cluster.Major & tableDS$IsoClusters[, 2] == Cluster.minor]
}
