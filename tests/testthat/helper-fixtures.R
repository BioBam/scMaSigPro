# Fixture builders shared by the test files and by fixtures/make-golden.R.
# Every object is built once per session and cached, so tests can call the
# builders freely without re-running the pipeline.

.fixture_cache <- new.env(parent = emptyenv())

fixture <- function(name, builder) {
  if (!exists(name, envir = .fixture_cache, inherits = FALSE)) {
    assign(name, builder(), envir = .fixture_cache)
  }
  get(name, envir = .fixture_cache, inherits = FALSE)
}

load_dataset <- function(name) {
  e <- new.env()
  utils::data(list = name, package = "scMaSigPro", envir = e)
  get(name, envir = e)
}

# Evaluate an expression, returning the error message instead of failing. Used
# by the golden snapshot so that options which currently error are frozen as
# errors (and a fix that makes them work shows up as a deliberate golden change).
freeze <- function(expr) {
  tryCatch(expr, error = function(e) list(.error = conditionMessage(e)))
}

# ---- splat.sim: bifurcation, 200 cells x 100 genes --------------------------

splat_scmp <- function() {
  fixture("splat_scmp", function() {
    as_scmp(
      object = load_dataset("splat.sim"), from = "sce",
      verbose = FALSE,
      additional_params = list(
        labels_exist = TRUE,
        exist_ptime_col = "Step",
        exist_path_col = "Group"
      )
    )
  })
}

splat_binned <- function() {
  fixture("splat_binned", function() sc.squeeze(splat_scmp()))
}

splat_design <- function(degree = 2) {
  fixture(paste0("splat_design_", degree), function() {
    sc.set.poly(splat_binned(), poly_degree = degree)
  })
}

# Three fitting variants frozen by the goldens (roadmap bug ledger, item h):
#   default   offset = TRUE, log_offset = FALSE  (the package default)
#   vignette  offset = TRUE, log_offset = TRUE   (what the Quick Start uses)
#   no_offset offset = FALSE
splat_fitted <- function(variant = c("default", "vignette", "no_offset")) {
  variant <- match.arg(variant)
  fixture(paste0("splat_fitted_", variant), function() {
    ob <- splat_design(2)
    ob <- switch(variant,
      default = sc.p.vector(ob, verbose = FALSE),
      vignette = sc.p.vector(ob, verbose = FALSE, offset = TRUE, log_offset = TRUE),
      no_offset = sc.p.vector(ob, verbose = FALSE, offset = FALSE)
    )
    sc.t.fit(ob, verbose = FALSE)
  })
}

splat_filtered <- function() {
  fixture("splat_filtered", function() {
    sc.filter(splat_fitted("default"),
      rsq = 0.7, vars = "groups", intercept = "dummy", includeInflu = TRUE
    )
  })
}

splat_clustered <- function() {
  fixture("splat_clustered", function() {
    sc.cluster.trend(splat_filtered(), geneSet = "union", k = 4)
  })
}

# A significant gene and a non-significant gene of the default fit, chosen
# deterministically.
splat_sig_gene <- function() {
  sort(unique(unlist(splat_filtered()@Significant@genes)))[1]
}
splat_nonsig_gene <- function() {
  ob <- splat_filtered()
  sort(setdiff(rownames(eDense(ob)), unlist(ob@Significant@genes)))[1]
}

# ---- multi.lin.sce: multifurcation with a shared root, 1977 cells x 501 genes

multi_cell_data <- function() {
  sce <- load_dataset("multi.lin.sce")
  cd <- as.data.frame(SingleCellExperiment::colData(sce))
  cd$group <- gsub("path_", "", cd$group, fixed = TRUE)
  cd$group <- gsub("|", "_links_", cd$group, fixed = TRUE)
  # Capitalised so it sorts first in every locale and becomes the reference
  # path (see the Common Root Cells vignette and roadmap item L1).
  cd$group[cd$group == "root"] <- "Root"
  list(
    counts = as.matrix(SingleCellExperiment::counts(sce)),
    cell_data = cd
  )
}

multi_scmp <- function() {
  fixture("multi_scmp", function() {
    x <- multi_cell_data()
    create_scmp(
      counts = x$counts, cell_data = x$cell_data,
      ptime_col = "Monocle3_Pseudotime", path_col = "group",
      use_as_bin = FALSE
    )
  })
}

multi_binned <- function() {
  fixture("multi_binned", function() sc.squeeze(multi_scmp()))
}

multi_restruct <- function() {
  fixture("multi_restruct", function() {
    sc.restruct(multi_binned(),
      end_node_list = list("Y_15", "Y_18", "Y_52"),
      root_node = "Root", link_node_list = list("Y_18_links_Y_15"),
      verbose = FALSE, link_sep = "_links_",
      assay_name = "counts", aggregate = "sum"
    )
  })
}

# ---- Golden snapshot ---------------------------------------------------------
# Everything a refactor must not change. fixtures/make-golden.R saves the
# result of build_golden(); test-golden.R compares a fresh build_golden()
# against the saved file.

sorted_df <- function(df) df[order(rownames(df)), , drop = FALSE]

dense_snapshot <- function(ob) {
  list(bins = as.data.frame(cDense(ob)), bulk = eDense(ob))
}

fit_snapshot <- function(ob) {
  list(
    p_values = ob@Profile@p_values,
    adj_p_values = ob@Profile@adj_p_values,
    fdr = ob@Profile@fdr,
    non_flat = ob@Profile@non_flat,
    offset = ob@Design@offset,
    sol = ob@Estimate@significance_matrix,
    coefficients = ob@Estimate@coefficient_matrix,
    path_coefficients = ob@Estimate@path_coefficient_matrix,
    t_scores = ob@Estimate@t_score_matrix,
    influential = ob@Estimate@influential,
    path = ob@Estimate@path
  )
}

layer_snapshot <- function(p) lapply(p$layers, function(l) l$data)

build_golden <- function() {
  binned <- splat_binned()

  methods <- c("Sturges", "Sqrt", "Rice", "Doane", "Freedman.Diaconis", "Scott.Normal")
  by_method <- lapply(methods, function(m) {
    freeze(dense_snapshot(sc.squeeze(splat_scmp(), bin_method = m)))
  })
  names(by_method) <- methods

  options <- list(
    split_bins = freeze(dense_snapshot(sc.squeeze(splat_scmp(), split_bins = TRUE))),
    prune_bins = freeze(dense_snapshot(sc.squeeze(splat_scmp(), prune_bins = TRUE))),
    # fill_gaps and drop_trails alone are no-ops on splat.sim; pair them with
    # an option that creates a gap (prune_bins) or unequal paths (Doane: 9/8).
    fill_gaps = freeze(dense_snapshot(sc.squeeze(splat_scmp(), prune_bins = TRUE, fill_gaps = TRUE))),
    drop_trails = freeze(dense_snapshot(sc.squeeze(splat_scmp(), bin_method = "Doane", drop_trails = TRUE))),
    mean = freeze(dense_snapshot(sc.squeeze(splat_scmp(), aggregate = "mean")))
  )

  designs <- lapply(1:3, function(d) {
    ob <- splat_design(d)
    list(
      predictors = predictors(ob),
      assignment = pathAssign(ob),
      groups = ob@Design@groups.vector
    )
  })
  names(designs) <- paste0("degree_", 1:3)

  fits <- list(
    default = fit_snapshot(splat_fitted("default")),
    vignette = fit_snapshot(splat_fitted("vignette")),
    no_offset = fit_snapshot(splat_fitted("no_offset"))
  )

  filters <- list()
  for (v in c("all", "each", "groups")) {
    for (i in c("none", "dummy", "all")) {
      filters[[paste(v, i, sep = "_")]] <- freeze(
        sc.filter(splat_fitted("default"),
          rsq = 0.7, vars = v, intercept = i, includeInflu = TRUE
        )@Significant@genes
      )
    }
  }

  clustered <- splat_clustered()
  plots <- list(
    bintile = plotBinTile(binned)$data,
    trend = layer_snapshot(plotTrend(clustered, feature_id = splat_sig_gene(), logs = TRUE)),
    trend_cluster = layer_snapshot(plotTrendCluster(clustered, verbose = FALSE)),
    intersect = plotIntersect(clustered, return = TRUE)
  )

  list(
    meta = list(
      frozen_version = "0.0.4",
      generated = format(Sys.time(), "%Y-%m-%d"),
      r_version = R.version.string,
      # Frozen as they are today, known to be wrong; these move in Phase 5.
      known_wrong = c(
        "squeeze/by_method/Freedman.Diaconis" = "ledger g: bin width used as count",
        "squeeze/by_method/Scott.Normal" = "ledger g: bin width used as count",
        "fits/default" = "ledger h: raw size factors as log-link offset",
        "fits/vignette" = "ledger h: offsets recomputed in sc.t.fit",
        "squeeze" = "ledger S3: breaks parsed from rounded cut() labels"
      )
    ),
    sparse_cells = sorted_df(cSparse(binned)),
    squeeze = list(
      default = dense_snapshot(binned),
      by_method = by_method,
      options = options
    ),
    designs = designs,
    fits = fits,
    filters = filters,
    clusters = unlist(clustered@Significant@clusters),
    multi = list(
      squeeze = dense_snapshot(multi_binned()),
      restruct = dense_snapshot(multi_restruct())
    ),
    plots = plots
  )
}
