# Figure 6 - SLE data
library(Seurat)
library(qs)
library(SCpubr)
library(ggpubr)
library(sccomp)

dc_lupus = qs::qread('lupus_DC.qsave')
# re-analysis of DC clusters
all.genes <- rownames(dc_lupus)
dc_lupus <- ScaleData(dc_lupus, features = all.genes)
dc_lupus <- FindVariableFeatures(dc_lupus)
dc_lupus <- RunPCA(dc_lupus, features = VariableFeatures(object = dc_lupus))
ElbowPlot(dc_lupus)
dc_lupus <- RunUMAP(dc_lupus, dims = 1:10)
Idents(dc_lupus) = dc_lupus$cell_type
DimPlot(dc_lupus, reduction = "umap")
dc_lupus <- FindNeighbors(dc_lupus, dims = 1:10)
dc_lupus <- FindClusters(dc_lupus, resolution = 0.5)
dc_lupus <- RunUMAP(dc_lupus, dims = 1:10)
DimPlot(dc_lupus, reduction = "umap")
dc_lupus <- RenameIdents(dc_lupus, '7' = 'pDC', '2' = 'pDC', '6' = "cDC1","0" = "cDC2","1" = "cDC2","3" = "cDC2","4" = "cDC2","5"="cDC2","8"="cDC2", "9"="cDC2")
dc_lupus$CellType = Idents(dc_lupus)

FeaturePlot(dc_lupus, features = c("CLEC9A","GZMB"))
genes <- list("cDC1" = c("IDO1","CLEC9A","CADM1","CLNK"), 
              "pDC"= c("GZMB","LILRA4","TCF4","IRF7"),
              "cDC2"=c("CD1C","CLEC10A","FCER1A"))

p <- SCpubr::do_DotPlot(sample = dc_lupus,cluster = TRUE, 
                        features = genes,dot.scale = 8,   legend.length = 10)
p
ggplot2::ggsave(plot = p, 'markers.pdf',width = 5, height= 4)

p <- SCpubr::do_DimPlot(sample = dc_lupus,plot.axes = TRUE, raster = TRUE, raster.dpi = 3048, pt.size = 9,
                        colors.use = c("cDC1" = "#e9d8a6",
                                       "pDC" = "#005f73",
                                       "cDC2" = "#bb3e03"))
ggplot2::ggsave(
  filename = "/Users/cohn/DendriticProject/Figures/Paper/RAW/6/umaps.pdf",
  plot = p,
  width = 5,
  height = 4,
  device = cairo_pdf,
  bg = "white"
)
qs::qsave(dc_lupus,'lupus_DC_process.qsave')


# PLD4 expression in pDCs
p=VlnPlot(dc_lupus, features ="ENSG00000166428", split.by='disease', sort='increasing',group.by='cell_type', slot = 'counts')
plot.data <- p$data
colnames(plot.data) <- c('expr','cell.type','group')
plot.data.filter = dplyr::filter(plot.data, cell.type == "plasmacytoid dendritic cell")
plot.data.filter$group = factor(plot.data.filter$group, levels = c("normal","systemic lupus erythematosus"))
p <- ggplot(plot.data.filter, aes(group, expr,fill = group)) +
  geom_boxplot(width = 0.5, outlier.shape = NA, color = "black", alpha = 0.8) +
  stat_summary(fun=median, geom="point",shape = 95,  size = 10, color="darkred")+
  scale_fill_manual(values=c("grey90", "#a2c5db")) +
  theme_bw() +
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        axis.text.y= element_text(size= 10, color = "black"),
        axis.text.x= element_text(size= 12, color = "black"),
        strip.placement = "outside",legend.position = "none",axis.title.x=element_blank(),
        strip.background = element_rect(fill = "gray95",colour="white")) +
  ylab("Expression (log1p)") +
  stat_compare_means(label.y = 2,method='wilcox.test') 
p  

# cell proportions
refine_metadata_levels <- function(seurat_data){
  for (i in base::colnames(seurat_data@meta.data)){
    if (base::is.factor(seurat_data@meta.data[[i]])){
      base::print(base::paste("Re-evaluating levels for a factor column", i))
      base::print(
        base::paste(
          "before:", base::paste(base::levels(seurat_data@meta.data[[i]]), collapse=", ")
        )
      )
      seurat_data@meta.data[[i]] <- base::droplevels(seurat_data@meta.data[[i]])  # need to drop levels of the removed values
      base::print(
        base::paste(
          "after:", base::paste(base::levels(seurat_data@meta.data[[i]]), collapse=", ")
        )
      )
    }
  }
  return (seurat_data)
}
dc_lupus = refine_metadata_levels(dc_lupus)
dc_lupus$disease <- relevel(dc_lupus$disease, "normal")

sccomp_result <- 
  dc_lupus |>
  sccomp_estimate( 
    formula_composition = ~ disease, 
    sample = "donor_id", 
    cell_group = "cell_type", 
    bimodal_mean_variability_association = TRUE,
    cores = 1 
  ) |> 
  sccomp_test()

sccomp_res_sig <- sccomp_result %>%
  mutate(signif = ifelse(c_FDR < 0.05, "FDR < 0.05", "FDR >= 0.05")) %>% arrange(c_effect)

plot = sccomp_result |> 
  sccomp_boxplot(factor = "disease")
library(ggplot2)
p = plot +   scale_fill_manual(values=c("#c3d9de"))+ 
  theme(plot.background = element_blank(),
        axis.line.x = element_line(color = 'black'),
        strip.text.x = element_text(size = 10),
        axis.ticks.length = unit(0.25, "cm"),
        panel.grid = element_blank(),
        axis.text =element_text(color = 'black', size = 15),
        panel.border = element_rect(colour = "black", fill=NA, size=0.25),
        axis.text.x =  element_text(angle=90, hjust = 1, size = 15),
        axis.text.y =  element_text(size = 10),
        plot.margin = unit(c(0,0,0,0), "cm"))
p


# milo analysis
sce <- as.SingleCellExperiment(dc_lupus)
sce$donor_id = as.character(sce$donor_id)
sce$disease <- factor(ifelse(sce$disease == "normal", "normal", "systemic lupus erythematosus"), levels=c("normal", "systemic lupus erythematosus"), ordered=TRUE)
library(SingleCellExperiment)

reducedDim(sce, "PCA") <- Embeddings(dc_lupus, "pca")
reducedDim(sce, "UMAP") <- Embeddings(dc_lupus, "umap")
library(miloR)

milo <- Milo(sce)

# build KNN graph 
milo <- buildGraph(milo, k = 30, d = 30, reduced.dim = "PCA")
# "We suggest using prop=0.1 for datasets of less than 30k cells. For bigger datasets using prop=0.05 should be sufficient (and makes computation faster)"
milo <- makeNhoods(milo, prop = 0.1, k = 30, d=30, refined = TRUE, 
                   reduced_dims = "PCA", refinement_scheme="graph")
plotNhoodSizeHist(milo)
# Counting cells in neighborhoods
milo <- countCells(milo, meta.data = data.frame(colData(milo)), samples="donor_id") #sample_ID #donor
head(nhoodCounts(milo))

# Defining experimental design
milo_design <- data.frame(colData(milo))[,c("donor_id", "disease")]  #sample_ID #disease__ontology_label "donor", "dis"
milo_design <- distinct(milo_design)
rownames(milo_design) <- milo_design$donor_id

# Computing neighborhood connectivity
milo <- calcNhoodDistance(milo, d=30, reduced.dim="PCA")
# Visualize results:
milo <- buildNhoodGraph(milo)

# run the differential abundance test
da_results <- testNhoods(milo, design = ~ disease, design.df=milo_design, reduced.dim="PCA")

da_results <- testNhoods(milo,
                         design = ~ disease, #disease__ontology_label #dis
                         design.df = milo_design,
                         fdr.weighting = 'graph-overlap', # make sure to use graph-overlap since we used graph-based approach for makeNhoods()
                         reduced.dim = 'PCA')
head(da_results)


da_results %>%
  arrange(SpatialFDR) %>%
  head() 
da_results <- annotateNhoods(milo, da_results, coldata_col = "CellType") #celltype_granular #Ident2
da_results$CellType <- ifelse(da_results$CellType < 0.7, "Mixed", da_results$CellType)
p = plotDAbeeswarm(da_results, group.by = "CellType")
ggplot2::ggsave(plot = p, 'milo.pdf',width = 6, height= 5)

library(scales)
plotNhoodGraphDA(milo, da_results, alpha=0.1) +
  scale_fill_gradient2(low="#070091",
                       mid="lightgrey",
                       high="#910000", 
                       name="log2FC",
                       limits=c(-5,5),
                       oob=squish) 

