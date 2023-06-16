# ### Create Test STAT GLM Implementataions
#
# ### Create Testable objects for testthat()
# # Load the library
# suppressPackageStartupMessages(library(maSigPro))
# suppressPackageStartupMessages(library(MASS))
#
# # Load object from the /data
# sce.test <- readRDS("extdata/sce.test.RDS")
#
# # Create a simple design file
# design.file <- as.data.frame(colData(sce.test))
# design.file$Path1 <- ifelse(design.file$Group == "Path1", 1, 0)
# design.file$Path2 <- ifelse(design.file$Group == "Path2", 1, 0)
# design.file$Rep <- c(1:nrow(design.file))
# design.file <- design.file[, c(5, 8, 6, 7), ]
#
# # Extract the counts
# gCounts <- as.matrix(sce.test@assays@data$counts)
# all(colnames(gCounts) == rownames(design.file$edesign))
#
# # Create Design file
# x <- capture.output(new.design.file <- make.design.matrix(
#     edesign = design.file, degree = 2,
#     time.col = 1, repl.col = 2
# ))
#
# # Running Pvector and TFit
# x <- capture.output(p.ob <- p.vector(
#     data = as.matrix(gCounts), Q = 0.4,
#     design = new.design.file,
#     counts = T, family = negative.binomial(10)
# ))
#
# # Save the object for testing later
# saveRDS(p.ob, "extdata/pvector_maSigPro.RDS")
#
# x <- capture.output(tfit <- T.fit(p.ob, step.method = "backward", alfa = 0.05))
#
# # Save the object for testing later
# saveRDS(tfit, "extdata/tfit_maSigPro.RDS")
#
# # Get the significant genes
# sigs <- get.siggenes(tfit, rsq = 0.05, vars = "groups")
#
# # Save object for testing later
# saveRDS(sigs, "extdata/sigs_maSigPro.RDS")
#
# ###############################
# # Sc MaSigPro Test ############
# ###############################
#
# # Load the libraray
# suppressPackageStartupMessages(library(scMaSigPro))
#
# # Load Object
# sce.test <- readRDS("data/sce.test.RDS")
#
# # Create Design
# scmp.obj <- sc_makeDesign(sce.object = sce.test)
# saveRDS(scmp.obj, file = "extdata/sc_makeDesign().RDS")
#
# # Run P vector
# scmp.obj <- sc_pVector(scmp.obj, p.vector.sig = 0.4, model.type = "nb")
# saveRDS(scmp.obj, file = "extdata/sc_pVector().RDS")
#
# # Compuet T-step
# scmp.obj <- sc_tFit(scmp.obj, covar.sig = 0.05)
# saveRDS(scmp.obj, file = "extdata/sc_tFit().RDS")
#
# # Extract Significant Genes
# scmp.obj <- sc.sel.features(scmp.obj,
#                             filter = "rsq",
#                             Q =  0.05
# )
# # SaveObject
# saveRDS(scmp.obj, file = "extdata/scmp_sigGenes_sc.sel.features().RDS")
