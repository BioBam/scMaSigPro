# ScMaSigPro

Implementation of MaSigPro for scRNA-Seq Data

[![Lint Code Base](https://github.com/spriyansh/scMaSigPro/actions/workflows/super-linter.yml/badge.svg)](https://github.com/spriyansh/scMaSigPro/actions/workflows/super-linter.yml) 

[![R-CMD-check](https://github.com/spriyansh/scMaSigPro/actions/workflows/r.yml/badge.svg)](https://github.com/spriyansh/scMaSigPro/actions/workflows/r.yml)


## Introduction

ScMaSigPro is an R package designed to analyze single-cell trajectory data, acting as an adaptation of [MaSigPro](https://www.bioconductor.org/packages/release/bioc/html/maSigPro.html). This package incorporates pseudotime as a temporal covariate into the polynomial Generalized Linear Model (GLM) used by MaSigPro. This feature empowers users to investigate differentially expressed genes across trajectory branches over pseudotime.

## Installation

### Step-1: Install Devtools
```
install.packages("devtools")
```

### Step-2: Install ScMaSigPro
```
devtools::install_github("spriyansh/scMaSigPro", ref = "document", auth_token = "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")
```

---
