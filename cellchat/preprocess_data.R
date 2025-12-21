preprocess_data <- function(sample1_dir, sample2_dir, project_name, n_features = 2000,
                           dims_for_cluster = 1:30, group_var = "stim",
                            resolution = 0.5, plot_file = "../Figure/part1/") {
  
  set.seed(123)

  setwd("C:/data/Dafei/raw/")
  
  # Read data from Sample 1 and Sample 2
  cc <- Read10X(sample1_dir)
  tt <- Read10X(sample2_dir)
  
  # Create Seurat objects
  c1 <- CreateSeuratObject(counts = cc, project = project_name)
  t1 <- CreateSeuratObject(counts = tt, project = project_name)
  
  # Set the "stim" metadata
  c1$stim <- "Ctrl"
  t1$stim <- "Treatment"
  
  # Merge the two Seurat objects
  ss <- merge(x = c1, y = t1)
  
  ss1<-as.SingleCellExperiment(ss)
  # Run scDblFinder
  ss2 <- scDblFinder(ss1,samples="stim",BPPARAM=SerialParam(RNGseed=123))
  table(ss2$scDblFinder.class)
  ss$scDblFinder.class <- ss2$scDblFinder.class
  
  ss<-subset(ss,subset=scDblFinder.class %in% "singlet")

  # Calculate percent.mt
  ss[["percent.mt"]] <- PercentageFeatureSet(ss, pattern = "^mt-")
  
  # Create QC plot and save it
  VlnPlot(ss, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 4)
  ggsave(file = paste0(plot_file,"1_QC_plot.pdf"),width=8,height=8)
  
  # Subset cells based on QC criteria
  ss <- subset(ss, subset = nFeature_RNA > 200 & nFeature_RNA < 7500 & percent.mt < 10 )
  
  # Normalize, find variable features, scale, and run PCA

  #ss <- SCTransform(ss, method = "glmGamPoi", vars.to.regress = "percent.mt", verbose = FALSE)
  ss <- NormalizeData(ss)
  ss <- FindVariableFeatures(ss, selection.method = "vst", nfeatures = n_features)
  ss <- ScaleData(ss)
  ss <- RunPCA(ss)
  
  # Run Harmony
  ss <- RunHarmony(ss, group.by.vars = group_var, reduction = "pca")
  
  # Run UMAP
  ss <- RunUMAP(ss, dims = dims_for_pca, reduction = "harmony")
  
  # Find neighbors and clusters
  ss <- FindNeighbors(ss, reduction = "harmony", dims = dims_for_pca)
  ss <- FindClusters(ss, resolution = resolution, algorithm=3)

  DimPlot(ss, reduction = "umap", label = TRUE, pt.size = 0.5) +scale_color_d3("category20")+ NoLegend()
  ggsave(filename=paste0(plot_file,"umap_test.pdf"),width=6,height=5)

  module_scores <- list(
  B_cells. = c("Cd19", "Cd22", "Cd79a", "Cd79b"),
  CD3_T_cells. = c("Cd3e", "Trac", "Cd3d", "Cd3e", "Cd3g"),
  CD4_T_cells. = c("Cd3e", "Cd3g", "Cd4", "Ctla4"),
  CD8_T_cells. = c("Cd3e", "Cd8a", "Cd8b1"),
  Tregs. = c("Cd3e", "Cd4", "Foxp3"),
  Monocytes. = c("Itgam", "Csf1r", "Gsr", "Ly6c2", "Spn"),
  Macrophages. = c("Itgam", "Csf1r", "Adgre1", "Mrc1", "Cd14", "Fcgr3"),
  NK_cells. = c("Ncr1", "Klrb1c", "Gzma", "Klra4", "Nkg7"),
  NKT_cells. = c("Cd3e", "Ncr1", "Klrb1c", "Nkg7"),
  DCs. = c("Itgam", "Itgax"),
  Neutrophils. = c("S100a8", "S100a9", "Mmp9", "Csf3r", "Cxcl3"),
  RBCs. = c("Hba-a1"),
  MDSCs. = c("Ly6g", "Ly6a"))

 for (module in names(module_scores)) {
   ss <- AddModuleScore(ss, features = module_scores[[module]], name = module)
   FeaturePlot(ss, features = paste0(module, "1"), min.cutoff = "q5")
   ggsave(filename=paste0(plot_file,"celltype_score/",module,".pdf"),width=6,height=5)
  }


  unique_features <- unique(c("Cd45","Cd19", "Cd22", "Cd79a","Cd79b","Cd3e", "Trac", "Cd3d", "Cd3e", "Cd3g","Cd3e", "Cd3g","Cd4", "Ctla4",
  "Cd3e", "Cd8a", "Cd8b1","Cd3e","Cd4", "Foxp3","Itgam", "Csf1r", "Gsr", "Ly6c2" , "Spn","Itgam", "Csf1r", "Adgre1", "Mrc1","Cd14",
  "Fcgr3","Ncr1", "Klrb1c", "Gzma", "Klra4", "Nkg7","Cd3e", "Ncr1", "Klrb1c","Nkg7","Itgam", "Itgax","S100a8", "S100a9", "Mmp9", 
  "Csf3r","Cxcl3","Hba-a1","Mki67", "Ly6a","Foxp3"))

  f1<-FeaturePlot(ss,features=unique_features,min.cutoff="q5")
  ggsave(filename=paste0(plot_file,"geneexpression.png"),f1,width=25,height=40)
  ggsave(filename=paste0(plot_file,"geneexpression.pdf"),f1,width=25,height=40)

  f2<-DotPlot(ss,features=unique_features)+RotatedAxis()
  ggsave(filename=paste0(plot_file,"geneexpression_dotplot.png"),f2,width=10,height=5)
  ggsave(filename=paste0(plot_file,"geneexpression_dotplot.pdf"),f2,width=10,height=5)



  #new.cluster.ids <- c("BCs", "Macro", "Neutrophils", "CD8-T.1", "CD4-T", "NK","CD8-T.2", "Fibroblast", 
  #                   "B-like", "DCs","CD8-T.3", "Granulocyte")

  new.cluster.ids <- c("BCs", "Neutrophils",  "CD8_TCs1", "Macro_1","NK", "Macro_2", "CD4_TCs1", "CD4_TCs2","CD8_TCs2",
                       "NK_like","TIICs","Monocyte","DCs","B_like", "IFIt+_BCs", "rest_BCs")

  names(new.cluster.ids) <- levels(ss)
  ss <- RenameIdents(ss, new.cluster.ids)
  ss$celltypes<-Idents(ss)

  Idents(ss)<-"celltypes"
  #levels(ss)<-c("BCs", "B-like",  "CD8-T.1", "CD8-T.2","CD8-T.3", "CD4-T", "NK", "Fibroblast", 
  #                   "DCs","Macro", "Neutrophils", "Granulocyte")
  levels(ss)<-c("BCs", "B_like", "IFIt+_BCs", "rest_BCs",  "CD4_TCs1", "CD4_TCs2","CD8_TCs1","CD8_TCs2", 
                "NK_like","NK", "DCs","TIICs","Neutrophils","Macro_1","Macro_2", "Monocyte")

  cluster.levels<-c("BCs", "B_like", "IFIt+_BCs", "rest_BCs",  "CD4_TCs1", "CD4_TCs2","CD8_TCs1","CD8_TCs2", 
                "NK_like","NK", "DCs","TIICs","Neutrophils","Macro_1","Macro_2", "Monocyte")



  ss$celltypes<-Idents(ss)


  DimPlot(ss)+scale_color_d3("category20")
  ggsave(filename=paste0(plot_file,"2_UMAP_cluster_nolabel.pdf"),width=6,height=5)

  DimPlot(ss,label=T)+scale_color_d3("category20")
  ggsave(filename=paste0(plot_file,"2_UMAP_cluster_label.pdf"),width=6,height=5)

  DimPlot(ss,group.by="stim")+scale_color_d3("category20")
  ggsave(filename=paste0(plot_file,"3_UMAP_stim.pdf"),width=6,height=5)

  DimPlot(ss,split.by="stim")+scale_color_d3("category20")
  ggsave(filename=paste0(plot_file,"4_UMAP_split_nolabel.pdf"),width=12,height=6)

  DimPlot(ss,split.by="stim",label=T)+scale_color_d3("category20")
  ggsave(filename=paste0(plot_file,"4_UMAP_split_label.pdf"),width=12,height=6)

  
 # levels(ss@meta.data$celltypes)<-c("BCs", "B-like",  "CD8-T.1", "CD8-T.2","CD8-T.3", "CD4-T", "NK", "Fibroblast", 
 #                     "DCs","Macro", "Neutrophils", "Granulocyte")

  cell.prop<- as.data.frame(prop.table(table(ss@meta.data$celltypes,ss@meta.data$stim)))
  cell.prop$Freq<-cell.prop$Freq*100  
  colnames(cell.prop)<-c("clusters","type","Percentages")
  levels(cell.prop$clusters)<-cluster.levels
  
  ct1<-ggplot(cell.prop,aes(type,Percentages,fill=clusters))+
  geom_bar(stat = "identity",position="fill")+
  ggtitle("")+
  theme_bw()+
  theme(axis.ticks.length=unit(0.5,"cm"))+
  guides(fill=guide_legend(title=NULL))+scale_fill_d3("category20")
  ggsave(filename=paste0(plot_file,"5_cellproportion.pdf"),ct1,width=4,height=5)

  ct2<-ggplot(cell.prop,aes(clusters,Percentages,fill=type))+
  geom_bar(stat = "identity",position="fill")+
  ggtitle("")+
  theme_bw()+
  theme(axis.ticks.length=unit(0.5,"cm"))+
  guides(fill=guide_legend(title=NULL))+scale_fill_npg()+ RotatedAxis()
  ggsave(filename=paste0(plot_file,"5_cellproportion_2.pdf"),ct2,width=8,height=6)

  ct3 <- ggplot(cell.prop, aes(clusters, Percentages, fill = type)) +
    geom_bar(stat = "identity", position = position_dodge()) +
    ggtitle("") +
    theme_bw() +
    RotatedAxis() +
    theme(axis.ticks.length = unit(0.5, "cm"),panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank()) +
    guides(fill = guide_legend(title = NULL)) +
    scale_fill_npg()
  ggsave(filename=paste0(plot_file,"5_cellproportion_3.pdf"),ct3,width=8,height=6)

  color_palette <- pal_d3("category20")(20)
  

  ct4 <-  plot_stat(ss, plot_type = "group_count",group_by="celltypes",pal_setup=color_palette)
   ggsave(filename=paste0(plot_file,"5_cellproportion_4.pdf"),ct4,width=8,height=6)


  dir.create(paste0(plot_file,"/celltype_score/"))


  cell_types_DEG<-FindAllMarkers(ss, logfc.threshold = 0.3, test.use = "wilcox", min.pct=0.2, only.pos=T,slot = "data")
  
  write.csv(cell_types_DEG,paste0(plot_file,"DEG.csv"))

  cell_types_DEG %>%
    group_by(cluster) %>%subset(subset=avg_log2FC>=1&p_val_adj<0.05)%>%
    top_n(n = 100, wt = avg_log2FC) -> topMarkerGenes

  DoHeatmap(ss, label=F, features = topMarkerGenes$gene, group.bar=T, group.colors=eval(pal_d3("category20")(15))) +
                                             scale_fill_gradient2(low = 'navy', mid = 'white', high = 'darkred') +
                                             theme(axis.text.y = element_blank())

  ggsave(filename=paste0(plot_file,"DEheatmap.png"),width=5,height=6)
  ggsave(filename=paste0(plot_file,"DEheatmap.pdf"),width=5,height=6)

 


  ##################enrichment analysis

  dir.create(paste0(plot_file,"/enrichment/"))

  ##GO
  ids=bitr(topMarkerGenes$gene,"SYMBOL","ENTREZID","org.Mm.eg.db")
  topDEG<-merge(topMarkerGenes,ids,by.x="gene",by.y="SYMBOL")
  gcSample=split(topDEG$ENTREZID,topDEG$cluster)
  gcSample<-gcSample[cluster.levels]
  xx<-compareCluster(gcSample,fun="enrichGO",OrgDb="org.Mm.eg.db",pvalueCutoff=0.05)
  go<-simplify(xx)
  clusterProfiler::dotplot(go, showCategory = 2, font.size = 8)+RotatedAxis()
  ggsave(filename=paste0(plot_file,"/enrichment/","GO.png"),width=6,height=8)
  ggsave(filename=paste0(plot_file,"/enrichment/","GO.pdf"),width=6,height=8)


  ##KEGG
  msig = msigdbr::msigdbr(species = "Mus musculus", category = "C2")
  msig = data.frame(msig$gs_name[which(msig$gs_subcat == "CP:KEGG")], msig$gene_symbol[which(msig$gs_subcat == "CP:KEGG")])
  gcKEGG=split(topMarkerGenes$gene,topMarkerGenes$cluster)
  gcKEGG<-gcKEGG[cluster.levels]
  cp = clusterProfiler::compareCluster(geneCluster = gcKEGG, fun = "enricher", TERM2GENE = msig)
  clusterProfiler::dotplot(cp, showCategory = 2, font.size = 8)+RotatedAxis()
  ggsave(filename=paste0(plot_file,"/enrichment/","KEGG.png"),width=6,height=8)
  ggsave(filename=paste0(plot_file,"/enrichment/","KEGG.pdf"),width=6,height=8)

   
  results<-list(ss,cell_types_DEG)

  return(results)
}

