# ScMaSigPro

Implementation of MaSigPro for scRNA-Seq Data

[![Lint Code Base](https://github.com/spriyansh/scMaSigPro/actions/workflows/super-linter.yml/badge.svg)](https://github.com/spriyansh/scMaSigPro/actions/workflows/super-linter.yml) 

[![R-CMD-check](https://github.com/spriyansh/scMaSigPro/actions/workflows/r.yml/badge.svg)](https://github.com/spriyansh/scMaSigPro/actions/workflows/r.yml)


## Introduction

ScMaSigPro is an R package designed to analyze single-cell trajectory data, acting as an adaptation of [MaSigPro](https://www.bioconductor.org/packages/release/bioc/html/maSigPro.html). This package incorporates pseudotime as a temporal covariate into the polynomial Generalized Linear Model (GLM) used by MaSigPro. This feature empowers users to investigate differentially expressed genes across trajectory branches over pseudotime.

## Installation

### Step-1: Installation
```
# Install devtools if not already installed
if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")

# Install scMaSigPro (token will be removed once the directory is public)
devtools::install_github("spriyansh/scMaSigPro",
                         ref = "main",
                         auth_token = "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"
                         )
```
## Basic Usage
The basic workflow of scMaSigPro involves the following steps:

### Set Seed for Reproducibility
```
set.seed(123)
```

### Load the scMaSigPro Package and Dataset

```
library(scMaSigPro)
data("Sim2Path", package = "scMaSigPro")
```

### Convert to scMaSigPro Object

```
scmp.sce <- as_scmp(object = sim.sce, from = "sce", ...)
```

### Pseudo-Bulk Processing

```
scmp.sce <- squeeze(scmp.sce, ...)
```

Design Matrix Creation

```
scmp.sce <- sc.make.design.matrix(scmp.sce, poly_degree = 2)
```

Run P-vector Calculation

```
scmp.sce <- sc.p.vector(scmp.sce, ...)
```

T-fit Calculation

```
scmp.sce <- sc.T.fit(scmp.sce, ...)
```

Gene Selection with R²

```
scmp.sce <- sc.get.siggenes(scmpObj = scmp.sce, ...)
```

Gene Trend Visualization

```
sc.PlotGroups(scmpObj = scmp.sce, feature_id = "Gene138", ...)
```

For detailed instructions, please refer to the vignette included in the package, which contains comprehensive guidance and explanations for each step in the analysis process.

# ... continue with the workflow steps as outlined above ...
For the full example and explanation of each function, please see the vignette.

Contributing
We welcome contributions from the community, including bug reports, suggestions, and pull requests. Please read CONTRIBUTING.md for details on our code of conduct, and the process for submitting pull requests to us.

License
This project is licensed under the MIT - see the LICENSE.md file for details.

Citation
If you use scMaSigPro in your research, please cite:

Priyansh Srivastava (2023).


