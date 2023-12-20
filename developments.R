suppressPackageStartupMessages(library(scMaSigPro))
suppressPackageStartupMessages(library(maSigPro))


# Step-1: Load Data
data("data.abiotic")
data("edesign.abiotic")

# Step-2: Set-up data for scMaSigPro
count <- as.matrix(data.abiotic)
cell_metadata <- as.data.frame(edesign.abiotic)

# Step-2.1: Add group column
cell_metadata$Group <- apply(cell_metadata[, c(3:6)], 1, FUN = function(x) {
  return(names(x[x == 1]))
})

# Step-2.2: Remove Binary Columns
cell_metadata <- cell_metadata[, !(colnames(cell_metadata) %in% c("Control", "Cold", "Heat", "Salt")), drop = FALSE]

# Step-3: Create scmp Object
test.scmp <- create.scmp(
  counts = count,
  cell_data = cell_metadata,
  pseudotime_colname = "Time",
  path_colname = "Group",
  use_as_bin = T
)

# Step-5: Set polynomial and make design
test.scmp.2 <- sc.set.poly(test.scmp, poly_degree = 2)


# Step-7: Run sc.p.vector
test.scmp.2 <- sc.p.vector(test.scmp.2,
  min.na = 20, verbose = TRUE,
  offset = FALSE, parallel = TRUE, max_it = 25,
  epsilon = 0.00001, family = gaussian()
)
