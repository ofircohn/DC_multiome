# Figure 3 and 4
library(CytoTRACE2) 
library(readr)
library(ggplot2)
library(ggrepel)
library(scales)
library(scattermore)
library(stringr)
library(ccaPP)
DC <- qs::qread("integrated.multiomics.chromvar.clustered.qsave")
DC <- FindVariableFeatures(DC)
cytotrace2_result <- cytotrace2(DC,   
                                species = "human",
                                is_seurat = TRUE,
                                slot_type = "counts",
                                full_model = TRUE,
                                batch_size = 100000,
                                smooth_batch_size = 10000,
                                parallelize_models = F,
                                parallelize_smoothing = F,
                                ncores = 1,
                                max_pcs = 200,
                                seed = 14)  

annotation <- data.frame(phenotype = cytotrace2_result@meta.data$CellTypeAll) %>% set_rownames(., colnames(cytotrace2_result))

meta <- as.data.frame(cytotrace2_result@meta.data)
df <- data.table(Embeddings(cytotrace2_result, "wnn.umap"),cytotrace2_result@meta.data)
readr::write_tsv(df, "meta.tsv")

# plot umap

lo <- quantile(meta$CytoTRACE2_Relative, 0.02, na.rm = TRUE)
hi <- quantile(meta$CytoTRACE2_Relative, 0.98, na.rm = TRUE)
centroids <- meta %>%
  group_by(CellTypeAll) %>%
  summarize(wnnUMAP_1 = mean(wnnUMAP_1), wnnUMAP_2 = mean(wnnUMAP_2))

p= ggplot(meta) +
  scattermore::geom_scattermore(
    aes(wnnUMAP_1, wnnUMAP_2, color = CytoTRACE2_Relative),
    alpha = 0.95, pointsize = 9, pixels = c(3000, 3000)
  ) +
  geom_text(
    data = centroids,
    aes(x = wnnUMAP_1, y = wnnUMAP_2, label = CellTypeAll),
    size = 3
  ) +
  scale_color_gradientn(
    name = "Differentiation\norder",
    colours = c("#2c7bb6", "#ffffbf", "#d7191c"),
    limits = c(lo, hi),
    oob = scales::squish,
    guide = guide_colorbar(ticks.colour = "black",frame.colour = "black",barwidth = 1.5)
  ) +
  theme_void()

ggplot2::ggsave(
  filename = "umap_cytrotrace.pdf",
  plot = p,
  width = 12,
  height = 8,
  units = "in",
  device = cairo_pdf,
  bg = "white",
  dpi = 600
)

# cytotrace2-gene expression correaltion
cytotrace2_score <- FetchData(cytotrace2_result, vars = "CytoTRACE2_Score")[,1]
gene_expression_data <- GetAssayData(
  cytotrace2_result,
  assay = "RNA",
  slot = "data"   # or "counts" if you want raw counts
)
all(colnames(gene_expression_data) == names(cytotrace2_score))
ct2genes <- sapply(
  1:nrow(gene_expression_data),
  function(x) ccaPP::corPearson(gene_expression_data[x, ], cytotrace2_score)
)
names(ct2genes) <- rownames(gene_expression_data)

ct2genes_df <- data.frame(
  gene = names(ct2genes),
  cor = as.numeric(ct2genes)
)

ct2genes_df <- ct2genes_df[order(ct2genes_df$cor, decreasing = TRUE), ]

write_tsv(ct2genes_df, "ct2genes_cor.tsv")

RNA = qs::qread("integrated_clustered_RNA_diet.qsave")
RNA <- FindVariableFeatures(RNA, assay = "RNA", selection.method = "vst")
var_genes <- VariableFeatures(RNA)
ct2genes_df_var <- ct2genes_df[ct2genes_df$gene %in% var_genes, ]
ct2genes_df_var <- ct2genes_df_var[order(ct2genes_df_var$cor, decreasing = TRUE), ]
DefaultAssay(RNA) = "RNA"
ct2genes_df$rank <- seq_len(nrow(ct2genes_df))
highlight_genes <- c("KIT","CD34","MECOM","NFIL3","ZEB2","ID2","BATF3","ZBTB18","KLF4","MYCL","BCL11A","TCF4","SPI1")
highlight_df <- ct2genes_df[ct2genes_df$gene %in% highlight_genes, ]


p = ggplot(ct2genes_df, aes(x = rank, y = cor)) +
  geom_line(color = "#E4572E", linewidth = 0.7, alpha = 0.8) +
  geom_point(data = highlight_df, color = "#1B9E77", size = 2.6) +
  geom_text_repel(
    data = highlight_df,
    aes(label = gene),
    color = "#1B9E77",
    size = 3.5,
    fontface = "bold",
    max.overlaps = Inf
  ) +
  labs(x = "Progeny rank", y = "Pearson correlation") +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6),
    axis.line = element_line(color = "black", linewidth = 0.5),
    axis.ticks = element_line(color = "black", linewidth = 0.5),
    axis.ticks.length = unit(3, "pt"),
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black")
  )
ggplot2::ggsave(
  filename = "progeny_rank.pdf",
  plot = p,
  width = 6,
  height = 4,
  device = cairo_pdf,
  bg = "white"
)


## box plot 
# adapted from https://github.com/digitalcytometry/cytotrace2
#mtd <- seurat@meta.data[c("Phenotype_CT2", "CytoTRACE2_Score")]
mtd = readr::read_tsv('meta.tsv') %>% select(CytoTRACE2_Score,CellTypeAll)

medians <- mtd %>%
  group_by(CellTypeAll) %>%
  dplyr::summarise(CytoTRACE2_median_per_pheno = median(CytoTRACE2_Score, na.rm = TRUE)) %>%
  dplyr::arrange(dplyr::desc(CytoTRACE2_median_per_pheno))

mtd <- mtd %>%
  inner_join(medians, by = "CellTypeAll")

mtd$CellTypeAll <- factor(mtd$CellTypeAll, levels = medians$CellTypeAll)
labels <- c('Differentiated', 'Unipotent', 'Oligopotent', 'Multipotent', 'Pluripotent', 'Totipotent')
colors <- c("#9E0142", "#F46D43", "#FEE08B", "#E6F598", "#66C2A5", "#5E4FA2")
potencyBoxplot_byPheno <- ggplot(mtd[!is.na(mtd$CellTypeAll), ], aes(x = CellTypeAll, y = CytoTRACE2_Score)) +
  geom_boxplot(aes(fill = CytoTRACE2_median_per_pheno), width = 0.8, alpha = 0.5, outlier.shape = NA) +
  #geom_jitter(aes(fill = CytoTRACE2_median_per_pheno), width = 0.05, height = 0, alpha = 0.5, shape = 16, stroke = 0.1, size = 1) +
  theme_classic() +
  scale_y_continuous(breaks = seq(0, 1, by = 0.2), limits = c(0, 1),
                     sec.axis = sec_axis(trans = ~., breaks = seq(0, 1, by = 1/12),
                                         labels = c("", 'Differentiated', "", 'Unipotent', "", 'Oligopotent', "", 'Multipotent', "", 'Pluripotent', "", 'Totipotent', ""))) +
  scale_fill_gradientn(colors = rev(colors), breaks = c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0), limits = c(0, 1), labels = labels) +
  scale_color_gradientn(colors = rev(colors), breaks = c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0), limits = c(0, 1), labels = labels) +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 10)) +
  labs(x = "", y = 'Potency score') +
  theme(legend.position = "None", axis.text = element_text(size = 8), axis.text.x = element_text(angle = 45, vjust = 0.5),
        axis.title = element_text(size = 12), legend.text = element_text(size = 12),
        plot.title = element_text(size = 12, face = "bold", hjust = 0.5, margin = margin(b = 20)),
        axis.ticks.y.right = element_line(color = c("black", NA, "black", NA, "black", NA, "black", NA, "black", NA, "black", NA, "black")),
        axis.ticks.length.y.right = unit(0.3, "cm"))
ggsave("potency.pdf", potencyBoxplot_byPheno, width = 8, height = 4, units = "in",useDingbats = FALSE)
