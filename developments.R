library(tidyverse)
library(assertthat)
library(scMaSigPro)
library(shiny)
library(igraph)
library(plotly)

# cds <- readRDS("/supp_data/Analysis_Public_Data/rep1/rep1_cds.RDS")


m3_select_path(cds,
  redDim = "umap",
  annotation_col = "cell_type",
  path_col = "Path",
  pseudotime_col = "Pseudotime",
  use_shiny = F,
  plot_purity = T
  # m3_pp = list(
  #     root_pp = c("Y_31", "Y_71"),
  #         path1_pp = c("Y_82", "Y_48"),
  #     path2_pp = c("Y_22", "Y_149")
  #                  )
)
