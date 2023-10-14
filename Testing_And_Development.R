##########################################################
## Author: Priyansh Srivastava ###########################
## Email: spriyansh29@gmail.com ##########################
## Script: scMaSigPro Application ########################
##########################################################

set.seed(007)

# Load scMaSigPro
library(scMaSigPro)
library(ggplot2)

# Load Objects
load("../scMaSigPro_Supp/benchmarks/11_RealDataSmall/data/results/monocle3_inferred_pseudotime.RData")
load("../scMaSigPro_Supp/benchmarks/01_Sparsity/data/simulated/sce/sparsity_30.RData")

# SCE
sce.anno <- annotate_sce(sce = sim.sce,
             existing_pseudotime_colname = "Step",
             existing_path_colname = "Group",
             path_prefix = "Path",
             path_colname = "Path",
             overwrite_labels = T,
             )
View(as.data.frame(SingleCellExperiment::colData(sce.anno)))

# Convert the ScMaSigPro Object
scmp.obj.sce <- as_scmp(sim.sce, from = "sce",
                        path_colname = "a",
                        pseudotime_colname = "Pseudotime",
                        additional_params = list(overwrite_labels = T,
                                                 existing_pseudotime_colname = "Step",
                                                 existing_path_colname = "Group")
                        )
View(as.data.frame(SingleCellExperiment::colData(scmp.obj.sce@sce)))




scmp.obj.cds <- as_scmp(cds, from = "cds")
View(as.data.frame(SingleCellExperiment::colData(scmp.obj.cds@sce)))


cds <- annotate_monocle3_cds(cds,root_label = "Progenitor",pseudotime_colname = "Goam")
View(as.data.frame(SingleCellExperiment::colData(cds)))

# Compress
scmp.obj <- squeeze(
    scmp.ob = scmp.obj,
    time.col = "Step",
    path.col = "Group",
    method = "Sturges",
    drop.fac = 0.6,
    verbose = T,
    cluster.count.by = "sum"
)

# Make Design
scmp.obj <- sc.make.design.matrix(scmp.obj,
                                  degree = 2,
                                  time.col = "binnedTime",
                                  path.col = "path"
)

# Run p-vector
scmp.obj <- sc.p.vector(
        scmpObj = scmp.obj, verbose = F, min.obs = 6,
        counts = T, theta = 1,
        offset = T, epsilon = 0.00001
        )
    
# Run-Step-2
scmp.obj <- sc.T.fit(
        data = scmp.obj, verbose = T,
        step.method = "backward",
        family = scmp.obj@scPVector@family,
        offset = T
)

# Extract the genes
sig.gene <- sc.get.siggenes(scmpObj = scmp.obj, rsq = 0.7, vars = "groups")

sc.PlotGroups(scmpObj = scmp.obj,feature_id = "GATA1", dis = scmp.obj@scTFit@dis,
                        edesign =  scmp.obj@scTFit@edesign,
                        groups.vector = scmp.obj@scTFit@groups.vector)
sc.plot.bins(scmpObj = scmp.obj)

saveRDS(scmp.obj, "../scMaSigPro_Supp/Testing_And_Development/scmp.obj.latest.RDS")
scmp.obj <- readRDS("../scMaSigPro_Supp/Testing_And_Development/scmp.obj.latest.RDS")


# Create path Segments
compress.meta <- as.data.frame(SingleCellExperiment::colData(scmp.obj@compress.sce))
expand.meta <- as.data.frame(SingleCellExperiment::colData(scmp.obj@sce))
expand.meta <- expand.meta[, c("Pseudotime", "Path", "PrincipalPoints", "predicted.celltype.l2")]
expand.meta$cluster.members <- rownames(expand.meta)

# Gene count
gene_count <- data.frame(raw_count = scmp.obj@sce@assays@data@listData$counts["GATA1",],
           cluster.members = colnames(scmp.obj@sce@assays@data@listData$counts))

# Extract the reuired coulms
compress.meta <- compress.meta[, c("binnedTime", "cluster.members", "bin", "path")]
compress.meta <- compress.meta %>% 
    separate_rows(cluster.members, sep = "\\|")

plot.data <- left_join(expand.meta, compress.meta,
                       by = c("cluster.members"))

plot.data <- left_join(plot.data, gene_count,
                       by = c("cluster.members"))


plot.data <- plot.data[order(plot.data$binnedTime),]

segments_df <- data.frame(
    x = plot.data$binnedTime[-nrow(plot.data)],
    xend = plot.data$binnedTime[-1],
    y = plot.data$raw_count[-nrow(plot.data)],
    yend = plot.data$raw_count[-1],
    col = plot.data$predicted.celltype.l2[-nrow(plot.data)]
)

# Plot
# Plot
p <- ggplot(plot.data, aes(x = binnedTime, y = raw_count)) +
    geom_point(aes(color = path)) +
    geom_segment(data = segments_df, aes(x = x, xend = xend, y = y, yend = yend)) +
    theme_minimal() +
    labs(title = "Raw Count vs. Binned Time", x = "Binned Time", y = "Raw Count", color = "Path")

print(p)



