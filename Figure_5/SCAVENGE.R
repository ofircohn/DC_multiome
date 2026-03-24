# Figure 5 - SCAVENGE

library(SCAVENGE)
library(chromVAR)
library(gchromVAR)
library(SummarizedExperiment)
library(data.table)
library(BiocParallel)
library(dplyr)
library(igraph)
library(parallel)
library(readxl)
library(dplyr)
library(purrr)
library(ComplexHeatmap)
library(circlize)
library(Matrix)
library(cluster)
library(RColorBrewer)
library(tidyverse)
library(ggpubr)
library(ggplot2)
library(tidyverse)
library(patchwork)
library(BSgenome.Hsapiens.UCSC.hg38)
set.seed(1234)


# export SingleCellExperiment object
DefaultAssay(DC) <- "ATAC"
DC.sce <- as.SingleCellExperiment(DC, assay = "ATAC")
peak.matrix <- GetAssayData(object = DC, slot = "counts")
idx.keep <- rowSums(x = peak.matrix) > 0
peak.matrix <- peak.matrix[idx.keep, , drop = FALSE]
peak.ranges <- granges(x = DC)
peak.ranges <- peak.ranges[idx.keep]
chromvar.obj <- SummarizedExperiment::SummarizedExperiment(
  assays = list(counts = peak.matrix),
  rowRanges = peak.ranges,
  colData = colData(DC.sce)
)

#Add UMAP coordinates to column data in summarized experiment
embeds = as.data.frame(Embeddings(data[["wnn.umap"]]))
colData(chromvar.obj)$UMAP1 <- embeds$wnn.umap
colData(chromvar.obj)$UMAP2 <- embeds$wnn.umap
chromvar.obj <- chromVAR::addGCBias(
  object = chromvar.obj,
  genome = BSgenome.Hsapiens.UCSC.hg38
)
# Remove NA values https://github.com/GreenleafLab/chromVAR/issues/26
row.data <- data.frame(SummarizedExperiment::rowData(x = chromvar.obj))
row.data[is.na(x = row.data)] <- 0
SummarizedExperiment::rowData(x = chromvar.obj) <- row.data
bg <- chromVAR::getBackgroundPeaks(
  object = chromvar.obj,niterations=200)


dir_path <- "path/to/finemapped/data"
for (file_path in list.files(dir_path,full.names=T)) {
  tryCatch({
    file_name <- basename(file_path)  
    # print progress
    print(paste0("Processing: ", file_name))
    trait <- importBedScore(rowRanges(chromvar.obj), file_path, colidx = 5) # finngen 5
    SE_blood_DEV <- computeWeightedDeviations(chromvar.obj, trait, background_peaks = bg)
    ## Reformat results
    z_score_mat <- data.frame(colData(chromvar.obj), z_score=t(assays(SE_blood_DEV)[["z"]]) %>% c)
    z_score_mat$z_score[is.nan(z_score_mat$z_score)]<-0 
    z_score_mat$z_score[which(!is.finite(z_score_mat$z_score))] <- 0
    seed_idx <- seedindex(z_score_mat$z_score, 0.05)
    ## Calculate scale factor
    scale_factor <- cal_scalefactor(z_score=z_score_mat$z_score, 0.01)
    ## Construct m-knn graph
    peak_by_cell_mat <- assay(chromvar.obj)
    tfidf_mat <- tfidf(bmat=peak_by_cell_mat, mat_binary=TRUE, TF=TRUE, log_TF=TRUE)
    lsi_mat <- do_lsi(tfidf_mat, dims=30)
    ## Calculate m-knn graph
    mutualknn30 <- getmutualknn(lsi_mat, 30)
    ### Network propagation  
    np_score <- randomWalk_sparse(intM=mutualknn30, rownames(mutualknn30)[seed_idx], gamma=0.05)
    ## Trait relevant score (TRS) with scaled and normalized 
    omit_idx <- np_score==0
    sum(omit_idx)
    mutualknn30 <- mutualknn30[!omit_idx, !omit_idx]
    np_score <- np_score[!omit_idx]
    TRS <- np_score %>% capOutlierQuantile(., 0.95) %>% max_min_scale
    TRS <- TRS * scale_factor
    zscoreWeighted1 <- data.frame(z_score_mat[!omit_idx, ], seed_idx[!omit_idx], np_score, TRS)
    permut <- get_sigcell_simple(knn_sparse_mat=mutualknn30, seed_idx=zscoreWeighted1$seed_idx, topseed_npscore=zscoreWeighted1$np_score, permutation_times=1000, true_cell_significance=0.05, rda_output=F, mycores=3, rw_gamma=0.05)
    zscoreWeighted2 <- data.frame(zscoreWeighted1, permut)
    openxlsx::write.xlsx(as.data.frame(zscoreWeighted2), paste0("path/to/out/", file_name, ".xlsx"), rowNames=F)
  }, error = function(e) {
    # Skip to the next file if there is an error
    message("Skipping file: ", file_name)
  })
  
}


# plots
# set a working dir
dir = "/TRS/outputs"

# merge all excel files and delete strings
file.list <- list.files(path = dir, pattern='*.xlsx',full.names=T)
df.list <- lapply(file.list, read_excel)
attr(df.list, "names") <- file.list
df2 <- bind_rows(df.list, .id = "id")
df2$id <-gsub('.{9}$', '', df2$id)
df2$id <-gsub('.FINEMAP.snp', '', df2$id)
df2$id <-gsub(dir, '', df2$id)
df2$id <-gsub("/", '', df2$id)
df2$id <-gsub("TRS", '', df2$id)
df2$id <-gsub(".bed.xlsx", '', df2$id)
df2$id = gsub(" ", "_", df2$id)

colnames(df2)[which(names(df2) == "CellTypeAll")] <- "cell_cluster"

# remove duplicates and non-defined diseases
remove = unique(c("MR05751","MR06525","MR11759","MR12973","PE06525",
                  "MR14191","PE07093","PE06880","PE06403","PE06293","PE06125","PE06124","PE04994","PE05693",
                  "PE05621","PE05476","PE05474","PE05056","PE04953","PE04950", "PE00356","PE06126","PE05501","GCST003738",
                  "PE04876","PE04658","PE04552","PE04548","PE01720","PE00315","GCST90029016","MR14174","PE04316",
                  "PE00228","PE00215","PA01086","PA01060","PA01059","PA00763","PE04949","PE06275","PE04987",
                  "PA00646","PA00523","MR35913","MR35644","MR34447","MR17165","PE05623","PE05626","PE06518","PE06506","PE04994",
                  "PE06507","GD09446","GD09037","PE06752","PA01075","PE06512","GD09507","GD00840","GD00851"))

# create a matrix
df_summarize = df2 %>% dplyr::select(id, cell_cluster, TRS, z_score) %>% 
  dplyr::group_by(id,cell_cluster) %>% dplyr::summarize(TRS_mean = mean(TRS), TRS_median = median(TRS), z_score_mean = mean(z_score), z_score_median = median(z_score) )
df_summarize = df_summarize %>% filter(.,!is.na(cell_cluster))
dim_x <- unique(df_summarize$id)
dim_y <- unique(df_summarize$cell_cluster)
map_x <- setNames(seq_along(dim_x), dim_x)
map_y <- setNames(seq_along(dim_y), dim_y)
mat <- as.matrix(sparseMatrix(
  i=map_x[df_summarize$id], j=map_y[df_summarize$cell_cluster], x=df_summarize$TRS_mean,   #or TRS_median TRS_mean
  dims=c(length(dim_x), length(dim_y)), dimnames=list(dim_x, dim_y)
))

# plot heatmap
column_order= c("pDC_1","pDC_2","pDC_3","pDC_4","pDC_5","pDC_6",
                "MEP / Ery / MAST_1","MEP / Ery / MAST_2","MEP / Ery / MAST_3","T /NK /B",
                "mregDC","cDC1_1","cDC1_2",
                "cDC2_1","cDC2_2","cDC2_3","cDC2_4","cDC2_5","cDC2_6","cDC2_7","cDC2_8 / DC3","cDC2_9",
                "cDC2_10","cDC2_11")

mat_scaled = mat[ , column_order]
mat_scaled = t(scale(t(as.matrix(mat_scaled))))
fa = rep(c("group A", "group B", "group C", "group D", "group E"), times = c(6,4,1,2,11))
### new plot
order1 <- seriation::seriate(dist(mat_scaled),method="OLO_ward")
###
mat_scaled2 <- mat_scaled[, order(fa)]

p1 = ComplexHeatmap::Heatmap(
  mat_scaled2,   
  width  = ncol(mat_scaled2)*unit(3, "mm"),
  height = nrow(mat_scaled2)*unit(3, "mm"),
  name="TRS",
  col = circlize::colorRamp2(c(0, max(mat_scaled2)), c("white", "#7d191f")),
  column_names_side = "top",
  cluster_columns = FALSE, 
  clustering_distance_rows = "euclidean",
  clustering_method_rows ="complete",
  column_split = fa,
  top_annotation = HeatmapAnnotation(foo = anno_block(gp = gpar(fill = c("#d9cb84", "#6f9ed6","#a9cc8f","#e0799f","#e3b586")),
                                                      labels_gp = gpar(col = "white", fontsize = 10))),

  row_names_side = "left",
  show_column_dend = FALSE,
  border = "black",
  row_names_gp = gpar(fontsize=8),
  column_names_gp = gpar(fontsize=8),
  heatmap_legend_param = list(
    title_gp   = gpar(fontsize = 8, fontface = "bold"),
    labels_gp  = gpar(fontsize = 8),
    border     = "black"
  )
)

draw(p1)


# UMAP plot
disease.list = c('BE109_Keratosis_Seborrheic','GD09159_Dermatitis_Atopic',"GD03661_Soft_Tissue_Infections")
library(scattermore)
plot.list <- list()
for (i in 1:length(disease.list)) {
  plot.list[[i]] <- ggplot(data=subset(df2,id == disease.list[i]), aes(UMAP1, UMAP2, color=TRS))+
    geom_scattermore(alpha = 0.95, pointsize = 9, pixels    = c(3000, 3000)) +
    xlab("UMAP 1") +  
    ylab("UMAP 2") +
    labs(colour = "TRS") +
    ggtitle(disease.list[i])+
    scale_color_gradientn(
      colors = c("gray95", "#bfa1a3", "#a61174"),
      oob = scales::squish,
      guide = guide_colorbar(
        ticks.colour = "black",
        frame.colour = "black",
        barwidth = 1.3
      ))+
    theme_bw() +
    theme(
      legend.background = element_rect(color = NA),
      axis.title = element_blank(),           
      axis.line = element_blank(),          
      strip.placement = "outside",
      strip.background = element_rect(fill = "gray95", colour = "white"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      plot.background =  element_rect(fill = "white", color = "white"),
      panel.background = ggplot2::element_rect(fill = "white", color = "white"),
      panel.border = ggplot2::element_rect(fill = NA, color = "black", linewidth = 1),
      plot.margin = ggplot2::margin(t = 0, r = 0, b = 0, l = 0),
      panel.grid.minor.y = ggplot2::element_line(color = "white", linewidth = 1))}
final_plot <- ggpubr::ggarrange(plotlist = plot.list, ncol=1, nrow=3, common.legend = TRUE, legend="right")
ggplot2::ggsave(
  filename = "umaps.pdf",
  plot = final_plot,
  width = 3,
  height = 6,
  device = cairo_pdf,
  bg = "white"
)


# enrichment by major DC populations
####################### ########### ########### ########### ########### 
Enriched = df2 %>% dplyr::group_by(cell_cluster,id) %>% 
  dplyr::summarise(enriched_cell=sum(true_cell_top_idx), TRS.mean= mean(TRS)) 
Enriched = Enriched[grepl("DC",Enriched$cell_cluster),]
Enriched$category = substr(Enriched$cell_cluster, 1, 4)
Enriched <- Enriched   %>% group_by(category,id) %>% dplyr::summarise(enriched_cell=sum(enriched_cell))  %>% 
  arrange(desc(enriched_cell))  %>% dplyr::slice(1:10)  %>%  ungroup()
Enriched$category <- as.factor(Enriched$category)
genes.list = c("cDC1","cDC2" ,"pDC_")
plot.list <- list()
for (gene in genes.list) {
  subset_data <- subset(Enriched, category == gene)
  plot <- ggplot(subset_data,aes(y =reorder(id,enriched_cell),fill =enriched_cell, x = enriched_cell))+
    geom_dotplot(binaxis = "y", stackdir = "center", position = "dodge",
                 stackratio = 1, dotsize = 1.5, binwidth=0.75)+
    ylab("") + xlab("Cell enrichment (n)") +
    scale_fill_gradient(low="#e3b7ac",high='#b53718',na.value="white") +
    theme_bw() +
    theme(axis.text.x = element_text(size=8, hjust=1, color="black"),
          axis.text.y = element_text(size=8, color="black"),
          axis.title.y = element_blank(),
          legend.position = "none",
          panel.spacing.x = unit(0.25, "cm"))+
    ggtitle(gene)
  plot.list[[gene]] <- plot
}

height_cDC1 <- length(unique(subset(Enriched, category == "cDC1")$id))
height_cDC2 <- length(unique(subset(Enriched, category == "cDC2")$id))
height_pDC_ <- length(unique(subset(Enriched, category == "pDC_")$id))
relative_heights <- c(height_cDC1, height_cDC2, height_pDC_)
combined_plot <- (plot.list[["cDC1"]] /
                    plot.list[["cDC2"]] /
                    plot.list[["pDC_"]]) +
  patchwork::plot_layout(ncol = 1, heights = relative_heights)


ggsave(
  filename = "enriched.splited.pdf",
  plot = combined_plot,
  width = 4.5,
  height = 6 # The height here now acts as a total height constraint
)

ggpubr::ggarrange(plotlist = plot.list, ncol=1, nrow=3) %>% 
  ggpubr::ggexport(filename = "enriched.splited.integrated.pdf", width = 2, height = 6)
