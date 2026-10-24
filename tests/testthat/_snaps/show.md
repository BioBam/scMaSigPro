# show() output is stable at every stage

    Code
      print(splat_scmp())
    Output
      Class: ScMaSigProClass
      nCells: 200
      nFeatures: 100
      Pseudotime Range: 1 100
      Branching Paths: Path2, Path1

---

    Code
      print(splat_binned())
    Output
      Class: ScMaSigProClass
      nCells: 200
      nFeatures: 100
      Pseudotime Range: 1 100
      Branching Paths: Path2, Path1
      Binned Pseudotime: 1-8(Range), 4.5(Mean), 
      Number of bins-> Path1: 8 Path2: 8
      Average bin Size-> Path1: 12 Path2: 13

---

    Code
      print(splat_design(2))
    Output
      Class: ScMaSigProClass
      nCells: 200
      nFeatures: 100
      Pseudotime Range: 1 100
      Branching Paths: Path2, Path1
      Binned Pseudotime: 1-8(Range), 4.5(Mean), 
      Number of bins-> Path1: 8 Path2: 8
      Average bin Size-> Path1: 12 Path2: 13
      Polynomial Order: 2

---

    Code
      print(splat_fitted("default"))
    Output
      Class: ScMaSigProClass
      nCells: 200
      nFeatures: 100
      Pseudotime Range: 1 100
      Branching Paths: Path2, Path1
      Binned Pseudotime: 1-8(Range), 4.5(Mean), 
      Number of bins-> Path1: 8 Path2: 8
      Average bin Size-> Path1: 12 Path2: 13
      Polynomial Order: 2
      No. of Significant Profiles: 52
      No. of Influential Features: 10

---

    Code
      print(splat_clustered())
    Output
      Class: ScMaSigProClass
      nCells: 200
      nFeatures: 100
      Pseudotime Range: 1 100
      Branching Paths: Path2, Path1
      Binned Pseudotime: 1-8(Range), 4.5(Mean), 
      Number of bins-> Path1: 8 Path2: 8
      Average bin Size-> Path1: 12 Path2: 13
      Polynomial Order: 2
      No. of Significant Profiles: 52
      No. of Influential Features: 10

---

    Code
      print(multi_restruct())
    Output
      Class: ScMaSigProClass
      nCells: 1977
      nFeatures: 501
      Pseudotime Range: 0 62.806
      Branching Paths: Y_15, Y_52, Y_18
      Binned Pseudotime: 1-28(Range), 13.39(Mean), 
      Number of bins-> Y_15: 23 Y_18: 24 Y_52: 19
      Average bin Size-> Y_15: 35 Y_18: 49 Y_52: 50

# showPoly() prints the full model formula

    Code
      showPoly(splat_design(1))
    Output
      [1] "beta0 + beta1*Path2vsPath1 + beta2*scmp_binned_pseudotime + beta3*scmp_binned_pseudotimexPath2"

---

    Code
      showPoly(splat_design(2))
    Output
      [1] "beta0 + beta1*Path2vsPath1 + beta2*scmp_binned_pseudotime + beta3*scmp_binned_pseudotimexPath2 + beta4*scmp_binned_pseudotime2 + beta5*scmp_binned_pseudotime2xPath2"

---

    Code
      showPoly(splat_design(3))
    Output
      [1] "beta0 + beta1*Path2vsPath1 + beta2*scmp_binned_pseudotime + beta3*scmp_binned_pseudotimexPath2 + beta4*scmp_binned_pseudotime2 + beta5*scmp_binned_pseudotime2xPath2 + beta6*scmp_binned_pseudotime3 + beta7*scmp_binned_pseudotime3xPath2"

# showParams() lists every parameter used so far

    Code
      params
    Output
                 parameters                  value
      1       bin_ptime_col scmp_binned_pseudotime
      2         path_prefix                   Path
      3          root_label                   root
      4           ptime_col             Pseudotime
      5    selection_method               backward
      6          bin_method                Sturges
      7            path_col                   Path
      8             bin_col               scmp_bin
      9        bin_size_col          scmp_bin_size
      10        bin_mem_col       scmp_bin_members
      11           anno_col              cell_type
      12                  g                    100
      13            p_value                   0.05
      14             min_na                      6
      15      mt_correction                     BH
      16               link                    log
      17            epsilon                  1e-08
      18             offset                   TRUE
      19             max_it                    100
      20        poly_degree                      2
      21     cluster_method                 hclust
      22            use_dim                    col
      23            fill_na                   zero
      24         min_counts                      0
      25 min_bins_expressed                      0
      26       distribution  Negative Binomial(10)

