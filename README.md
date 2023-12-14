# ScMaSigPro

Implementation of MaSigPro for scRNA-Seq Data

[![Lint Code Base](https://github.com/spriyansh/scMaSigPro/actions/workflows/super-linter.yml/badge.svg)](https://github.com/spriyansh/scMaSigPro/actions/workflows/super-linter.yml) 

[![R-CMD-check](https://github.com/spriyansh/scMaSigPro/actions/workflows/r.yml/badge.svg)](https://github.com/spriyansh/scMaSigPro/actions/workflows/r.yml)

---

## Introduction

`scMaSigPro` is an R package designed for analyzing single-cell RNA-seq data over pseudotime. Building on the [maSigPro](https://www.bioconductor.org/packages/release/bioc/html/maSigPro.html) package, it identifies genes with significant expression changes across branching paths in a pseudotime-ordered dataset. This guide provides a step-by-step workflow for ScMaSigPro, making it accessible for users.

## Installation
To install `scMaSigPro` from GitHub, use the following R code:

```
# Install devtools if not already installed
if (!requireNamespace("devtools", quietly = TRUE)) {
    install.packages("devtools")
}

# Install scMaSigPro
devtools::install_github("spriyansh/scMaSigPro",
                         ref = "dev",
                         auth_token = "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",
                         build_vignettes = FALSE,
                         build_manual = TRUE,
                         upgrade = "never",
                         force = TRUE,
                         quiet = TRUE)
```

---

## Basic Usage
The basic workflow of `scMaSigPro` involves the following steps:

### 1. Load Package and Dataset
```
set.seed(123)
library(scMaSigPro)
data("splat.sim", package = "scMaSigPro")
```

### 2. Create `scMaSigPro` Object

```
# Helper Function to convert annotated SCE object to scmpObject
scmp.ob <- as_scmp(
  object = splat.sim, from = "sce",
  align_pseudotime = FALSE,
  verbose = TRUE,
  additional_params = list(
    labels_exist = TRUE,
    existing_pseudotime_colname = "Step",
    existing_path_colname = "Group"
  )
)
```

### 3. Pseudo-bulking with `squeeze()`

This function discretizes a continuous pseudotime column into bins:

```
scmp.sce <- squeeze(scmp.sce, ...)
```

### 4. Setting up the Polynomial Model

```
scmp.sce <- sc.make.design.matrix(scmp.sce, poly_degree = 2)
```

### 5. Detecting Genes with Non-Flat Profiles

```
scmp.sce <- sc.p.vector(scmp.sce, ...)
```

### 6. Model Refinement

```
scmp.sce <- sc.T.fit(scmp.sce, ...)
```

### 7. Gene selection with $R^2$

```
scmp.sce <- sc.get.siggenes(scmpObj = scmp.sce, ...)
```

### 8. Gene Trend Visualization

```
scmp.ob <- sc.cluster.features(scmp.ob, ...)
sc.PlotProfiles(scmp.ob, groupBy = "feature")
```

For detailed instructions, please refer to the quick start vignette included in the package, which contains comprehensive guidance and explanations for each step in the analysis process.

## Contributing
Contributions, including bug reports, suggestions, and pull requests, are welcome.


## License
This project is licensed under the MIT - see the LICENSE.md file for details.

## Citation
If you use `scMaSigPro` in your research, please cite:

---
