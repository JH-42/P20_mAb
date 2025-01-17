update_barcodes <- function(tcr_subset, pattern) {
  parts <- str_split(tcr_subset[["barcode"]], pattern, simplify = TRUE)
  tcr_subset[["barcode"]] <- as.data.frame(parts)$V2
  return(tcr_subset)
}

process_column <- function(data, column_name, target_field, add_before = "", add_after = "", remove_field = TRUE) {
  if (!column_name %in% colnames(data)) {
    stop(paste("Column", column_name, "not found in the data!"))
  }
  
  if (!any(str_detect(data[[column_name]], target_field))) {
    stop(paste("Target field", target_field, "not found in column", column_name, "!"))
  }
  
  data[[column_name]] <- str_replace(
    data[[column_name]], 
    paste0(target_field), 
    paste0(add_before)
  )
  
  data[[column_name]] <- paste0(data[[column_name]], add_after)
  
  return(data)
}

plot_and_save_clonal_network <- function(immexp_subset, identity, file_prefix, manualcolors) {
  for (i in c("CD8-Tn","CD8-Tprolif","CD8-Teff","Treg","Th17","Tfh","NKT")) {
    plot <- clonalNetwork(immexp_subset, reduction = "PHATE", group.by = identity,
                          filter.identity = i, cloneCall = "aa") + 
            scale_color_manual(values = manualcolors)
    ggsave(filename = paste0(file_prefix, i, "_colonalNetwork.pdf"), plot)
  }
}

plot_and_save_clonal_network_Bcells <- function(immexp_subset, identity, file_prefix, manualcolors,
  cells=c("CD8-Tn","CD8-Tprolif","CD8-Teff","Treg","Th17","Tfh","NKT")) {
  for (i in cells) {
    plot <- clonalNetwork(immexp_subset, reduction = "PHATE", group.by = identity,
                          filter.identity = i, cloneCall = "aa") + 
            scale_color_manual(values = manualcolors)
    ggsave(filename = paste0(file_prefix, i, "_colonalNetwork.pdf"), plot)
  }
}

plot_occupied_scRepertoire <- function(immexp_subset, file_name) {
  immexp_subset$new.cell <- factor(immexp_subset$new.cell, levels = levels(immexp_subset))
  plot <- clonalOccupy(immexp_subset, x.axis = "new.cell") + RotatedAxis()
  ggsave(file_name, plot)
}

plot_alluvial_clonotypes <- function(immexp_subset, file_name,color=c("#a4a839","#3e3d1a","#e94f9e","#705da9","#231d65","#d6a0c9","#c8ba8b")) {
  plot <- alluvialClones(immexp_subset, cloneCall = "gene", 
                             y.axes = c("group", "ident", "cloneSize"), 
                             color = "ident") +
          scale_fill_manual(values = color)
  ggsave(file_name, plot, width = 9, height = 6)
}

plot_chord_diagram <- function(immexp_subset, file_name, group_by, grid_cols) {
  circles <- getCirclize(immexp_subset, group.by = group_by)
  pdf(file = file_name)
  circlize::chordDiagram(circles, self.link = 1, grid.col = grid_cols, direction.type = "arrows")
  dev.off()
}