# Use while debugging
load_package <- function(){
    suppressPackageStartupMessages(library(assertthat))
    suppressPackageStartupMessages(library(SingleCellExperiment))
    suppressPackageStartupMessages(library(entropy))
    suppressPackageStartupMessages(library(tidyverse))
}