getDS <- function(Model, vars = "all", rsq = 0.4) {
  data <- Model$data
  gen <- Model$gen
  selected.genes <- Model$DSG
  Tfit2 <- Model$Tfit

  data2 <- data[gen %in% selected.genes, ]
  gen2 <- gen[gen %in% selected.genes]

  #--------------------------------------------------------------------------------------
  get2 <- get.siggenes(Tfit2, vars = vars, rsq = rsq)
  sig.iso2 <- get2$summary
  gen.sig.iso2 <- as.character(gen2[rownames(data2) %in% sig.iso2])
  NumIso.by.gene <- tapply(sig.iso2, gen.sig.iso2, length) # Num of Transcripts from selected Isoforms
  DSG_distributed_by_number_of_DETs <- NumIso.by.gene
  T.iso2 <- table(DSG_distributed_by_number_of_DETs, useNA = "ifany")

  #--------------------------------------------------------------------------------------
  print(paste(length(selected.genes), " DSG selected"))
  print(paste(length(gen.sig.iso2), " DETs selected"))
  print(T.iso2)

  List0 <- setdiff(selected.genes, gen.sig.iso2) # list of DSG without any DET

  out <- list(Model, get2, selected.genes, sig.iso2, List0, NumIso.by.gene)
  names(out) <- c("Model", "get2", "DSG", "DET", "List0", "NumIso.by.gene")
  out
}
