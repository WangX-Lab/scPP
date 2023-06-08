#' Title
#'
#' @param sc_dataset A seurat object of single cell RNA sequencing data.
#' @param geneList A gene list correlated with interested features.
#' @param probs Cutoff value of ScPP, default is 0.2.
#'
#' @return Meta data with information of ScPP selected cells.
#' @export
#'
#' @examples
ScPP = function(sc_dataset, geneList, probs = 0.2){
  if(length(geneList) != 2){
    stop("This gene list do not have enough information correlated with interested feature.")
  }

  if (missing(sc_dataset) || class(sc_dataset) != "Seurat")
    stop("'sc_dataset' is missing or not a seurat object.")

  library(AUCell)
  cellrankings = AUCell_buildRankings(sc_dataset@assays$RNA@data,plotStats = FALSE)
  cellAUC = AUCell_calcAUC(geneList,cellrankings)

  metadata = as.data.frame(sc_dataset@meta.data)
  metadata$AUCup <- as.numeric(getAUC(cellAUC)["gene_pos", ])
  metadata$AUCdown <- as.numeric(getAUC(cellAUC)["gene_neg", ])
  
  downcells1 = rownames(metadata)[which(metadata$AUCup <= quantile(metadata$AUCup,probs = probs))]
  upcells1 = rownames(metadata)[which(metadata$AUCup >= quantile(metadata$AUCup,probs = (1-probs)))]
  downcells2 = rownames(metadata)[which(metadata$AUCdown >= quantile(metadata$AUCdown,probs = (1-probs)))]
  upcells2 = rownames(metadata)[which(metadata$AUCdown <= quantile(metadata$AUCdown,probs = probs))]
  
  ScPP_neg = intersect(downcells1,downcells2)
  ScPP_pos = intersect(upcells1,upcells2)
  
  metadata$ScPP <- ifelse(rownames(metadata) %in% ScPP_pos, "Phenotype+", "Background")
  metadata$ScPP <- ifelse(rownames(metadata) %in% ScPP_neg, "Phenotype-", metadata$ScPP)
  
  return(metadata)
}
