# as_scmp rejects an unknown 'from'

    Code
      eval(cases[[nm]])
    Condition
      Error:
      ! 'from' must be "sce": as_scmp() converts SingleCellExperiment objects only. Monocle3 'cell_data_set' objects (from = "cds") are no longer converted; pass the counts and a cell table with pseudotime and path columns to create_scmp(), as shown in the 'Common Root Cells' vignette.

# as_scmp no longer converts Monocle3 objects (from = 'cds')

    Code
      eval(cases[[nm]])
    Condition
      Error:
      ! 'from' must be "sce": as_scmp() converts SingleCellExperiment objects only. Monocle3 'cell_data_set' objects (from = "cds") are no longer converted; pass the counts and a cell table with pseudotime and path columns to create_scmp(), as shown in the 'Common Root Cells' vignette.

# as_scmp rejects an object that is not a SingleCellExperiment

    Code
      eval(cases[[nm]])
    Condition
      Error:
      ! 'object' must be a SingleCellExperiment. To start from a count matrix and a table of cell metadata, use create_scmp().

# sc.squeeze rejects a non-ScMaSigPro input

    Code
      eval(cases[[nm]])
    Condition
      Error:
      ! Please provide object of class 'ScMaSigPro'.

# sc.squeeze rejects an unknown pseudotime column

    Code
      eval(cases[[nm]])
    Condition
      Error:
      ! 'nope' does not exist in cell.level.metadata. Please review the 'ptime_col' parameter.

# sc.squeeze rejects an unknown path column

    Code
      eval(cases[[nm]])
    Condition
      Error:
      ! 'nope' does not exist in cell.level.metadata. Please review the 'path_col' parameter.

# sc.squeeze rejects an unknown bin method

    Code
      eval(cases[[nm]])
    Condition
      Error:
      ! Available binning methods are 'Freedman.Diaconis', 'Sqrt', 'Sturges', 'Rice', 'Doane', and 'Scott.Normal'

# sc.squeeze rejects drop_fac below 0.3

    Code
      eval(cases[[nm]])
    Condition
      Error:
      ! Invalid value for 'drop_fac'. It should be between 0.3 and 1.

# sc.set.poly rejects a non-ScMaSigPro input

    Code
      eval(cases[[nm]])
    Condition
      Error:
      ! Please provide object of class 'ScMaSigPro'

# sc.p.vector rejects an unknown link

    Code
      eval(cases[[nm]])
    Condition
      Error:
      ! link function should be either 'log' or 'identity'

# sc.p.vector rejects min_na larger than the number of bins

    Code
      eval(cases[[nm]])
    Condition
      Error:
      ! 1000 for 'min_na' is too high. Try lowering the threshold.

# sc.t.fit rejects an unknown selection method

    Code
      eval(cases[[nm]])
    Condition
      Error in `sc.t.fit()`:
      ! stepwise method must be one of backward, forward, two.ways.backward, two.ways.forward

# sc.filter rejects an unknown vars

    Code
      eval(cases[[nm]])
    Condition
      Error:
      ! Invalid selection for 'vars'

# showSol errors before sc.t.fit

    Code
      eval(cases[[nm]])
    Condition
      Error:
      ! 'significance_matrix' is not computed yet

# showCoeff errors before sc.t.fit

    Code
      eval(cases[[nm]])
    Condition
      Error:
      ! 'coefficient_matrix' is not computed yet

# eSparse rejects an unknown assay

    Code
      eval(cases[[nm]])
    Condition
      Error in `eSparse()`:
      ! Assay 'nope' not found in the object.

# eDense rejects an unknown assay

    Code
      eval(cases[[nm]])
    Condition
      Error in `eDense()`:
      ! Assay 'nope' not found in the object.

# plotTrend rejects an unknown feature

    Code
      eval(cases[[nm]])
    Condition
      Error:
      ! 'NotAGene' doesn't exist please select another one.

# plotTrend refuses a non-significant feature by default

    Code
      eval(cases[[nm]])
    Condition
      Error:
      ! Feature Id didn't pass the R2 threshold, please re-run sc.filter() with a lower 'rsq' value or set 'significant' to 'FALSE'

# plotTrend needs at least one layer

    Code
      eval(cases[[nm]])
    Condition
      Error:
      ! At least one of 'curves', 'lines', or 'points' must be TRUE.

# plotBinTile errors before sc.squeeze

    Code
      eval(cases[[nm]])
    Condition
      Error:
      ! No binning information found. Please run 'sc.squeeze()', first.

# plotIntersect errors before sc.filter

    Code
      eval(cases[[nm]])
    Condition
      Error:
      ! 'Significant' slot is empty, please run 'sc.filter()' first.

# plotTrendCluster errors before sc.cluster.trend

    Code
      eval(cases[[nm]])
    Condition
      Error:
      ! Please run 'sc.cluster.trend', before plotting cluster trends

# sc.cluster.trend rejects an unknown gene set

    Code
      eval(cases[[nm]])
    Condition
      Error:
      ! 'nope' does not exist. Please use one of Path2vsPath1, Path1, intersect, union

# sc.cluster.trend rejects an unknown method

    Code
      eval(cases[[nm]])
    Condition
      Error:
      ! 'dbscan' is not a valid method. Please use one of hclust, kmeans, Mclust

# queryCoeff rejects an unknown query

    Code
      eval(cases[[nm]])
    Condition
      Error:
      ! Invalid query, please select one of 'pseudotime', 'pseudotime_path' or 'path'

# queryCoeff rejects an unknown change

    Code
      eval(cases[[nm]])
    Condition
      Error:
      ! Invalid change, please select one of 'increasing' or 'decreasing'

# create_scmp rejects mismatched cells

    Code
      eval(cases[[nm]])
    Condition
      Error:
      ! Number of cells in raw-counts and cell-level-metadata are different.

# sc.restruct errors before sc.squeeze

    Code
      eval(cases[[nm]])
    Condition
      Error:
      ! No binning information found. Please run 'sc.squeeze()', first.

