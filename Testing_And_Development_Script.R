# This document is temporary and intended for testing and debugging
# This will be removed from the final package

# Set seed
set.seed(123)

# Load ScMaSigpro
library(scMaSigPro)
library(ggplot2)

# Step-1: Load a dataset for testing
data("Sim2Path", package = "scMaSigPro")

# Step-2: Convert to ScMaSigpro Object
scmp <- as_scmp(object = sim.sce,
                from = "sce", 
                additional_params = list(labels_exist = TRUE,
                                         existing_pseudotime_colname = "Step",
                                         existing_path_colname = "Group")
                )

# Step-3-A: Perform Binning
scmp <- entropy_discretize(scmp,drop.fac = 1,
                           verbose = T,
                           binning = "individual",
                           additional_params = list(use_unique_time_points = TRUE))

# Step-3-B: Create Pseudo-Bulk Cell Metadata
scmp <- make.pseudobulk.design(scmp,
                               verbose= T)

# Step-3-C: Pseudo-bulk the counts
scmp <- make.pseudobulk.counts(scmp)

# Step-4: Make Design-Matrix
scmp <- sc.make.design.matrix(scmp, poly_degree = 2)

# Step-5: Run P-vector
scmp <- sc.p.vector(scmp, parallel = T, family = gaussian())

# Step-6: Run T.fit
scmp <- sc.T.fit(scmp, parallel = T, verbose = T)

# Step-7: Select with R2 
scmp <- sc.get.siggenes(scmpObj = scmp,
                        vars = "all",
                        significant.intercept = "dummy")

# Step-8: Plot Gene Trends
sc.PlotGroups(scmpObj = scmp,
              feature_id = "Gene10", smoothness = 0.1)


# Developing Methods for monocle3

# Step-1: Load data and Monocle3 like function
load("extdata/rep1_processed.RData")
suppressPackageStartupMessages(library(shiny))
suppressPackageStartupMessages(library(plotly))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(igraph))
suppressPackageStartupMessages(library(assertthat))
suppressPackageStartupMessages(library(SingleCellExperiment))

# Step-2: Select the Paths
selectPath(obj = cds, annotation = "predicted.celltype.l2")


# Create SCMP Object
scmp.cds <-  as_scmp(cds, "cds")

# Select the Path
scmp.cds <- selectPath(scmp.cds,
           sel.path = c("Y_3", "Y_34"),
           pathCol = "Path",balance_paths = T, 
           pTimeCol = "Pseudotime", 
           plot_paths = T, verbose = T)

scmp.cds <- entropy_discretize(scmp.cds,drop.fac = 1,
                           verbose = T,
                           binning = "individual",
                           additional_params = list(use_unique_time_points = TRUE))

scmp.cds <- make.pseudobulk.design(scmp.cds,
                               verbose= T)
scmp.cds <- make.pseudobulk.counts(scmp.cds)

# Step-4: Make Design-Matrix
scmp.cds <- sc.make.design.matrix(scmp.cds, poly_degree = 2)

# Step-5: Run P-vector
scmp.cds <- sc.p.vector(scmp.cds, parallel = T)



# Plot the trajectory
pTime <- plot_cells(cds = cds,
           color_cells_by = "pseudotime",
           label_roots = T, cell_size = 2,
           trajectory_graph_segment_size = 2,
           label_principal_points = T) +
    theme(legend.position = "bottom")
cell.anno <- plot_cells(cds = cds,
                    color_cells_by = "predicted.celltype.l2",
                    label_roots = T, cell_size = 2,
                    trajectory_graph_segment_size = 2,
                    label_principal_points = T) +
    theme(legend.position = "bottom")

tr.graph <- plot_cells(cds = cds,
                       label_branch_points = T, 
                       label_leaves = T,
                       group_label_size = 3,
                        color_cells_by = "predicted.celltype.l2",
                        label_roots = T, cell_size = 0,
                        trajectory_graph_segment_size = 2,
                        label_principal_points = T) +
    theme(legend.position = "bottom")

partitions <- plot_cells(cds = cds,
                       color_cells_by = "partition",
                       label_roots = T, cell_size = 2,
                       trajectory_graph_segment_size = 1,
                       label_principal_points = T) +
    theme(legend.position = "bottom")

#Plot all together
ggpubr::ggarrange(pTime, cell.anno,
                  tr.graph, partitions,
                  ncol= 2, nrow =2)

# Step-3: Extract fan like structures
suppressPackageStartupMessages(library(igraph))

# Extract graph from monocle3
pgraph <- cds@principal_graph@listData$UMAP
end_nodes <- V(pgraph)[degree(pgraph) == 1]$name

extract_fans <- function(node) {
    node_index <- which(node_names == node)
    neighbors_indices <- neighbors(pgraph, node_index)
    neighbors_names <- node_names[neighbors_indices]
    
    valid_neighbors <- neighbors_names[neighbors_names %in% end_nodes | neighbors_names %in% candidate_nodes]
    
    # Check if there are enough valid_neighbors to form a combination
    if (length(valid_neighbors) < 3) {
        return(list(center=paste0(node, " (center)"), edges=NULL))
    }
    
    fans <- combn(valid_neighbors, 3, function(triplet) {
        edges_named <- sapply(triplet, function(edge_node) paste0(edge_node, " (edge-", which(triplet == edge_node), ")"))
        list(center=paste0(node, " (center)"), edges=edges_named)
    }, simplify = FALSE)
    
    unlist(fans, recursive = FALSE)
}



selectPath(obj = cds)


sample_node <- "Y_6"
# sample_fans <- extract_fans(sample_node)

cat("Central node:", sample_fans$center, "\n")
cat("Edges:", sample_fans$edges, "\n")


candidate_nodes <- V(pgraph)[degree(pgraph) >= 3]

extract_fans <- function(node) {
    node_index <- which(node_names == node)
    neighbors_indices <- neighbors(pgraph, node_index)
    neighbors_names <- node_names[neighbors_indices]
    
    fans <- combn(neighbors_names, 3, function(triplet) {
        edges_named <- sapply(triplet, function(edge_node) paste0(edge_node, " (edge-", which(triplet == edge_node), ")"))
        list(center=paste0(node, " (center)"), edges=edges_named)
    }, simplify = FALSE)
    
    unlist(fans, recursive = FALSE)
}


sample_node <- "Y_6"
sample_fans <- extract_fans(sample_node)

cat("Central node:", sample_fans$center, "\n")
cat("Edges:", sample_fans$edges, "\n")


Loop through the filtered fans and print
for (fan in filtered_fans) {
    cat("Central node:", fan$center, "\n")
    cat("Edges:", fan$edges, "\n\n")
}



### Create test datasets for scMaSigPro

## Test Reproducibility with MaSigPro

# Load
library(maSigPro)
library(scMaSigPro)

# Load data
data("data.abiotic")
data("edesign.abiotic")

# MaSigPro: Make Design Matrix
design <- make.design.matrix(edesign = edesign.abiotic, degree = 4)

# Convert Design
add_col <- function(x){
    return(names(x[x==1]))
}
edesign.abiotic.in <- as.data.frame(edesign.abiotic[, 1, drop = F])
Group <- c(unlist(apply(edesign.abiotic[, c(3:6)], 1, add_col, simplify = F), use.names = F))
edesign.abiotic.in$Group <- Group

# Create Design with scMaSigPro
test.scmp <- create_scmpObj(counts = data.abiotic,
            cell_data = edesign.abiotic.in,
            pseudotime_colname = "Time",
            path_colname = "Group",
            use_as_bin = T)

# Make Design
test.scmp <- sc.make.design.matrix(test.scmp,poly_degree = 4)

# Set -Pvector
gc <- capture.output(
    fit <- p.vector(data.abiotic, design, Q = 0.05, MT.adjust = "BH", min.obs = 20)
    )
gc <- NULL

test.scmp <- sc.p.vector(test.scmp,
                         min.obs = 20,
                         offset = F,
                         parallel = T,
                         family = gaussian())

# T step
gc <- capture.output(
    tstep <- T.fit(fit, step.method = "backward", alfa = 0.05)
)
gc <- NULL

test.scmp <- sc.T.fit(test.scmp,
                         parallel = T, verbose = F, offset = F)




# Sig Genes
sigs <- get.siggenes(tstep, rsq = 0.6, vars = "groups")

# Sc sig genes
test.scmp <- sc.get.siggenes(test.scmp, rsq = 0.6, vars = "groups")

# See Genes
see.genes(sigs$sig.genes$ColdvsControl, show.fit = T, dis =design$dis,
          cluster.method="hclust" ,cluster.data = 1, k = 9, )

STMDE66 <- data.abiotic[rownames(data.abiotic)=="STMDE66", ]
PlotGroups (STMDE66, edesign = edesign.abiotic)

PlotGroups (STMDE66, edesign = edesign.abiotic, show.fit = T,
            dis = design$dis, groups.vector = design$groups.vector)

PlotProfiles(data = data.abiotic
             )

atestthat::expect_equal(expected = sigs$summary,
                       object = test.scmp@sig.genes@summary)

testthat::expect_equal(expected = sigs$sig.genes,
                       object = test.scmp@sig.genes@sig.genes)
