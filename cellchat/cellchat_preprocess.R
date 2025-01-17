load_and_subset <- function(file_path, idents_value) {
  BM <- readRDS(file_path)
  Idents(BM) <- "group"
  BM_subset <- subset(BM, idents = idents_value)
  Idents(BM_subset) <- "new.cell"
  return(BM_subset)
}

process_cellchat <- function(BM, filename, in_nPatterns=in_nPatterns,out_nPatterns=out_nPatterns) {
  data.input <- GetAssayData(BM, assay = "RNA", slot = "data")
  labels <- Idents(BM)
  meta <- data.frame(group = labels, row.names = names(labels))
  cc <- createCellChat(object = BM)
  
  cc@DB <- CellChatDB.mouse
  cc <- subsetData(cc)
  cc <- identifyOverExpressedGenes(cc)
  cc <- identifyOverExpressedInteractions(cc)
  
  cc <- projectData(cc, PPI.mouse)
  cc <- computeCommunProb(cc, raw.use = FALSE)
  cc <- filterCommunication(cc, min.cells = 10)
  cc <- computeCommunProbPathway(cc)
  cc <- aggregateNet(cc)
  cc <- netAnalysis_computeCentrality(cc, slot.name = "netP")
  
  selectK(cc, pattern = "incoming")
  ggsave(paste0("selectk_", filename, "_in.png"))
  nPatterns = in_nPatterns
  cc <- identifyCommunicationPatterns(cc, pattern = "incoming", k = nPatterns)
  

  selectK(cc, pattern = "outgoing")
  ggsave(paste0("selectk_", filename, "_out.png"))
  nPatterns = out_nPatterns
  cc <- identifyCommunicationPatterns(cc, pattern = "outgoing", k = nPatterns)
  
 
  cc <- computeNetSimilarity(cc, type = "functional")
  saveRDS(cc, paste0(filename, "_cellchat.RDS"))
}
