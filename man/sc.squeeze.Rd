% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sc.squeeze.R
\name{sc.squeeze}
\alias{sc.squeeze}
\title{Pseudo-bulking with optimal number of Pseudotime based bins}
\usage{
sc.squeeze(
  scmpObj,
  ptime_col = scmpObj@Parameters@ptime_col,
  path_col = scmpObj@Parameters@path_col,
  bin_method = "Sturges",
  drop_fac = 1,
  verbose = FALSE,
  bin_mem_col = "scmp_bin_members",
  bin_col = "scmp_bin",
  bin_size_col = "scmp_bin_size",
  bin_ptime_col = "scmp_binned_pseudotime",
  split_bins = FALSE,
  prune_bins = FALSE,
  assay_name = "counts",
  drop_trails = FALSE,
  aggregate = "sum",
  fill_gaps = FALSE,
  additional_params = list(use_unique_time_points = FALSE)
)
}
\arguments{
\item{scmpObj}{An object of class \code{\link{ScMaSigPro}}.}

\item{ptime_col}{A character string representing the column name
for inferred Pseudotime values in 'Sparse' data. (Default is "Pseudotime").}

\item{path_col}{A character string representing the column name for branching
path assignment in 'Sparse' or 'Dense' data. (Default is `path_prefix`).}

\item{bin_method}{A character string representing the algorithm used for
binning. Available options: "Freedman.Diaconis",
"Sqrt", "Sturges", "Rice", "Doane", and "Scott.Normal". (Default = "Sturges")}

\item{drop_fac}{A numeric value specifying the factor by which to adjust the
number of bins if the initial binning results in too many/few bins.
(Default = 1).}

\item{verbose}{Print detailed output in the console. (Default is TRUE)}

\item{bin_mem_col}{A character string representing the name of the column in
which cells per bin are stored. (Default is "scmp_bin_members").}

\item{bin_col}{A character string representing the name of the column in which
bin labels are stored. (Default is "scmp_bin").}

\item{bin_size_col}{A character string representing the name of the column in
which bin sizes per bin are stored. (Default is "scmp_bin_size").}

\item{bin_ptime_col}{A character string representing the column name
for binned Pseudotime values in 'Dense' data.
(Default is "scmp_binned_pseudotime").}

\item{split_bins}{If bin sizes are greater than mean + sd, split the bin into
smaller bins by re-running the sc.squeeze() function. (Default = FALSE).}

\item{prune_bins}{If bin sizes are smaller than mean - sd, remove the bin.
(Default = FALSE).}

\item{assay_name}{Name of the Assay in sparse data from which the counts are
used. (Default = "counts").}

\item{drop_trails}{If the paths have different lengths of the binned pseudotime,
drop the bins from the path with more bins. (Default = FALSE).}

\item{aggregate}{A character string specifying the method to aggregate counts
within each cluster. Available options are 'mean' or 'sum'. (Default = "sum").}

\item{fill_gaps}{If corresponding bin is missing for a time-point, pull the
successive bins and fill the gaps.}

\item{additional_params}{Pass additional parameters as a named list.
See examples}
}
\value{
An object of class \code{\link{ScMaSigPro}}, with updated `Dense`
slot.
}
\description{
`sc.squeeze()` discretizes a continuous time series column into bins
of equal size using entropy-based binning method. It automatically calculates
the optimal number of bins using one of the supported methods.
}
\seealso{
\code{\link{estBinSize}}, \code{\link{discretize}},
\code{\link{create_range}}
}
\author{
Priyansh Srivastava \email{spriyansh29@gmail.com}
}
