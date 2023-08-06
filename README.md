---

ScMaSigPro

---

Implementation of MaSigPro for scRNA-Seq Data

[![Lint Code Base](https://github.com/spriyansh/scMaSigPro/actions/workflows/super-linter.yml/badge.svg)](https://github.com/spriyansh/scMaSigPro/actions/workflows/super-linter.yml)

[![R](https://github.com/spriyansh/scMaSigPro/actions/workflows/r.yml/badge.svg)](https://github.com/spriyansh/scMaSigPro/actions/workflows/r.yml)


## Introduction

ScMaSigPro is an adaptation of [MaSigPro](https://www.bioconductor.org/packages/release/bioc/html/maSigPro.html) written R, specifically engineered to analyze single-cell trajectory data. It utilizes the inferred Pseudotime as a temporal covariate in the polynomial Generalized Linear Model (GLM) of MaSigPro, enabling the exploration of differential genes across trajectory branches over pseudotime.

## Installation

### Step-1: Install Devtools
```
install.packages("devtools")
```

### Step-2: Install ScMaSigPro
```
devtools::install_github("spriyansh/scMaSigPro")
```
