# make-multi-lin-sce.R
#
# Provenance script for the bundled dataset `multi.lin.sce` (see ?multi.lin.sce).
# This script is NOT run during package build or check; it documents how the
# object in data/multi.lin.sce.RData was produced so it can be regenerated.
# Requires: monocle3 (GitHub), magrittr, igraph, SingleCellExperiment, and the
# tradeSeq simulation RDS downloaded from the URL in the comments below.
# Run from the package root: Rscript inst/scripts/make-multi-lin-sce.R

# Load Libraries
library(monocle3)
library(magrittr)

# Download dataset of multi furcating trajectory from tradeSeq Article
# wget "https://github.com/statOmics/tradeSeqPaper/raw/master/simulation/sim2_dyntoy_multifurcating_4/multifurcating_4.rds"

multi_ob <- readRDS(file = "data/multifurcating_4.rds")

# Counts
raw_counts <- as.matrix(t(multi_ob[["counts"]]))

# Cell Metadata
cell_metadata_data <- as.data.frame(multi_ob[["cell_info"]])
rownames(cell_metadata_data) <- cell_metadata_data$cell_id

# Gene Metadata
gene_metadata_data <- as.data.frame(multi_ob[["feature_info"]])
rownames(gene_metadata_data) <- gene_metadata_data$feature_id
gene_metadata_data[["gene_short_name"]] <- gene_metadata_data$feature_id

# Convert to Monocle3 CDS
cds <- new_cell_data_set(
    expression_data = raw_counts,
    cell_metadata = cell_metadata_data,
    gene_metadata = gene_metadata_data
)

# Basic Steps
## Normalize
cds <- preprocess_cds(cds, norm_method = "log",
                      method = "PCA",
                      num_dim = 20,
                      pseudo_count = 1,
                      scaling = TRUE,
                      verbose = FALSE)

## Reduce Dimensions
set.seed(123)
cds <- reduce_dimension(cds, reduction_method = "tSNE", verbose = FALSE,
                        preprocess_method = "PCA",
                        cores = 1)
cds <- reduce_dimension(cds, reduction_method = "UMAP", verbose = FALSE,
                        preprocess_method = "PCA",
                        cores = 1)

# Overwrite UMAP Slots with tSNE as learn_graph only works on UMAP
reducedDims(cds)[["UMAP"]] <- reducedDims(cds)[["tSNE"]]
plot_cells(cds) + labs(title = "tSNE", x = "tSNE 1", y = "tSNE 2")

## Compute Clusters
cds <- cluster_cells(cds, verbose = FALSE, random_seed = 123, resolution = 0.8)
plot_cells(cds, color_cells_by = "cluster", cell_size = 3)  +
    labs(title = "tSNE", xlab = "tSNE 1", ylab = "tSNE 2")
plot_cells(cds, color_cells_by = "partition", cell_size = 3) +
    labs(title = "tSNE", x = "tSNE 1", y = "tSNE 2")

# Learn Graph
cds <- learn_graph(cds, verbose = FALSE,
                   learn_graph_control = list(minimal_branch_len = 15,
                                              prune_graph=TRUE, ncenter=200))
plot_cells(cds, color_cells_by = "cluster", cell_size = 3,
           label_principal_points = TRUE) + labs(title = "tSNE", x = "tSNE 1", y = "tSNE 2")

# Infer Pseudotime
cds <- order_cells(cds, root_pr_nodes = "Y_51", verbose = FALSE)
p <- plot_cells(cds, color_cells_by = "pseudotime", cell_size = 2,trajectory_graph_color = "red",trajectory_graph_segment_size = 2,
           label_principal_points = TRUE) + scale_color_viridis_c() +
    labs(title = "Simulated Multifurcating Trajectory",subtitle = "Simulation: Dyntoy | Latent Space: t-SNE", x = "tSNE-1", y = "tSNE-2", color = "Monocle3 Pseudotime")+
    theme(legend.position = "bottom") + geom_point(inherit.aes = TRUE, alpha = 0.5, cex = 0)
save(p, file = "extdata/multifurcating_trajectory.RData")

## Follow Steps from
## https://statomics.github.io/tradeSeq/articles/Monocle.html#extracting-the-pseudotimes-and-cell-weights-for-tradeseq-1

# Get the closest vertice for every cell
y_to_cells <-  principal_graph_aux(cds)$UMAP$pr_graph_cell_proj_closest_vertex %>%
    as.data.frame()
y_to_cells$cells <- rownames(y_to_cells)
y_to_cells$Y <- y_to_cells$V1

# Get the root vertices
# It is the same node as above
root <- cds@principal_graph_aux$UMAP$root_pr_nodes

# Extract Mst
mst <- principal_graph(cds)$UMAP

# Get the other endpoints
endpoints <- names(which(igraph::degree(mst) == 1))
endpoints <- endpoints[!endpoints %in% root]

# For each endpoint
cellWeights <- lapply(endpoints, function(endpoint) {
    # We find the path between the endpoint and the root
    path <- igraph::shortest_paths(mst, root, endpoint)$vpath[[1]]
    path <- as.character(path)
    # We find the cells that map along that path
    df <- y_to_cells[y_to_cells$Y %in% path, ]
    df <- data.frame(weights = as.numeric(colnames(cds) %in% df$cells))
    colnames(df) <- endpoint
    return(df)
}) %>% do.call(what = 'cbind', args = .) %>%
    as.matrix()
rownames(cellWeights) <- colnames(cds)
colnames(cellWeights) <- paste("path",colnames(cellWeights), sep = "_")

# Subset for 3 paths
cellWeights <- cellWeights[, c("path_Y_18", "path_Y_52", "path_Y_15"), drop=FALSE]
cellWeights <- cellWeights[rowSums(cellWeights) != 0, ]

# Create Cell Data
cellData <- data.frame(
    cell_id = rownames(cellWeights),
    row.names = rownames(cellWeights)
)

# Create Cell Metadata
cellData[["group"]] <- apply(cellWeights, 1, FUN = function(x) {

    npath <- length(names(x[x == 1]))

    if(npath == 3){
        return("root")
    }else if(npath == 2){
        return(paste(names(x[x == 1]), collapse = "|"))
    }else{
        return(names(x[x == 1]))
    }
})

# Get counts and Pseudotime
counts <- as.matrix(cds@assays@data@listData$counts)
counts <- counts[, rownames(cellData), drop=FALSE]

# Get Pseudotime
pseudotime_vector <- pseudotime(cds)
cellData[["Monocle3_Pseudotime"]] <- pseudotime_vector[rownames(cellData)]

# Create SingleCellExperiment Object
multi.lin.sce <- SingleCellExperiment(assays = list(counts = counts),
                            colData = cellData)

# Add dimensionality reduction
redDim <- reducedDims(cds)[["tSNE"]]
redDim <- redDim[rownames(cellData), , drop=FALSE]
reducedDims(multi.lin.sce)[["TSNE"]] <- redDim

# Save Object
save(multi.lin.sce, file = "data/multi.lin.sce.RData")
tools::resaveRdaFiles("data/multi.lin.sce.RData")
