calculate_TME_signature_score <- function(input_file, output_file) {
  cellMarker <- data.table::fread("TME_signature_marker.csv", data.table = FALSE)
  colnames(cellMarker)[2] <- "celltype"
  cellMarker <- lapply(split(cellMarker, cellMarker$celltype), function(x) {
    dd = x$Metagene
    unique(dd)
  })
  expr0 <- data.table::fread(input_file, data.table = FALSE)
  expr = dplyr::distinct(expr0, V1, .keep_all = TRUE)
  rownames(expr) <- expr[,1]
  expr <- expr[,-1]
  expr <- as.matrix(expr)
  library(GSVA)
  gsva_data <- gsva(expr, cellMarker, method = "gsva", min.sz = 1)
  write.csv(gsva_data, file = output_file)
}
