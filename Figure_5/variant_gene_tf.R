# Figure 5 - variant-to-gene
library(readxl)
library(data.table)
library(GenomicRanges)
library(purrr)
library(dplyr)
library(ggplot2)
library(Seurat)
library(scCustomize)

# load finemapped files
files <- list.files(path = "path/to/finemapped/files", full.names = T,pattern = "\\.bed$")
ids_new <- file.path(
  "/Users/cohn/DendriticProject/FileMappedData/causaldb.V2/v2.0/immune",
  basename(ids)

# load scE2G predictions
scE2G_mDC = fread('encode_e2g_predictions_threshold0.177.tsv.gz')
scE2G_mDC.gr = scE2G_mDC %>% dplyr::select(chr, start, end, TargetGene, ABC.Score, CellType, Kendall, E2G.Score.qnorm, class) %>% makeGRangesFromDataFrame(., keep.extra.columns = TRUE)
scE2G_cDC1 = fread('encode_e2g_predictions_threshold0.177.tsv.gz')
scE2G_cDC1.gr = scE2G_cDC1 %>% dplyr::select(chr, start, end, TargetGene, ABC.Score, CellType, Kendall, E2G.Score.qnorm,class) %>% makeGRangesFromDataFrame(., keep.extra.columns = TRUE)
scE2G_cDC2 = fread('encode_e2g_predictions_threshold0.177.tsv.gz')
scE2G_cDC2.gr = scE2G_cDC2 %>% dplyr::select(chr, start, end, TargetGene, ABC.Score, CellType, Kendall, E2G.Score.qnorm, class) %>% makeGRangesFromDataFrame(., keep.extra.columns = TRUE)
scE2G_pDC = fread('encode_e2g_predictions_threshold0.177.tsv.gz')
scE2G_pDC.gr = scE2G_pDC %>% dplyr::select(chr, start, end, TargetGene, ABC.Score, CellType, Kendall, E2G.Score.qnorm,class) %>% makeGRangesFromDataFrame(., keep.extra.columns = TRUE)

merged_scE2G  =  rbind(as_tibble(scE2G_mDC.gr),as_tibble(scE2G_cDC1.gr),as_tibble(scE2G_cDC2.gr),as_tibble(scE2G_pDC.gr)) %>% 
  dplyr::select(seqnames,start,end, CellType, TargetGene,ABC.Score,E2G.Score.qnorm,Kendall,class)  %>% makeGRangesFromDataFrame(., keep.extra.columns = TRUE)

safe_fun <- safely(function(f) {
  dat <- read.table(f, header = FALSE)
  dat$sample <- gsub("\\.bed$", "", basename(f))
  names(dat) <- c("chrom", "start", "end", "region", "PP", "rs", "traits")
  
  dat <- makeGRangesFromDataFrame(dat, keep.extra.columns = TRUE)
  olap.peaks <- plyranges::join_overlap_left(dat, merged_ABC)
  
  olap.peaks.df <- olap.peaks %>% 
    as_tibble() %>% 
    na.omit() %>%
    #group_by(TargetGene, rs) %>% 
    #summarise(E2G.Score.qnorm = max(E2G.Score.qnorm), .groups = "drop") %>%
    mutate(traits = gsub("\\.bed$", "", basename(f)))
  
  as.data.frame(olap.peaks.df)
})

seg_datlist <- map_df(files, function(f) {
  res <- safe_fun(f)
  if (is.null(res$result)) return(NULL)   
  res$result
}, .progress = TRUE)


trait_scE2G = seg_datlist %>% select(rs, E2G.Score.qnorm, TargetGene, CellType, traits, class)

# add DEG
RNA <- qs::qread("integrated_clustered_RNA_diet.qsave")
Idents(RNA) <- RNA$group
DefaultAssay(RNA) <- "RNA"
all_markers_pct <- FindAllMarkers(object = RNA, verbose = TRUE) %>%
  Add_Pct_Diff() %>%
  select(cluster, gene, avg_log2FC, p_val_adj, pct.1, pct.2, pct_diff)
plot_df <- trait_scE2G %>%
  group_by(CellType, traits, TargetGene,rs, class) %>%
  summarise(
    E2G.Score.qnorm = max(E2G.Score.qnorm, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(all_markers_pct, by = c("CellType" = "cluster", "TargetGene" = "gene")) %>%
  filter(!is.na(E2G.Score.qnorm) & !is.na(avg_log2FC) & !is.na(pct_diff)) %>%
  mutate(
    score = scale(E2G.Score.qnorm)[,1] +
      scale(avg_log2FC)[,1] +
      scale(pct_diff)[,1]
  )

# plot top 20 target genes 
top_df <- plot_df %>%
  group_by(CellType, TargetGene, traits) %>%
  summarise(score = max(score, na.rm = TRUE), .groups = "drop") %>%
  filter(!is.na(score), !is.na(traits)) %>%          # drop missing
  group_by(CellType) %>%
  slice_max(score, n = 20, with_ties = FALSE) %>%
  ungroup()

trait_cols <- c(
  "CA022_Lupus_Erythematosus_Systemic"           = "#1f77b4",
  "CA075_Celiac_Disease"                         = "#ff7f0e",
  "CA152_Vitiligo"                               = "#2ca02c",
  "CA210_Crohn_Disease"                          = "#d62728",
  "CA211_Inflammatory_Bowel_Diseases"            = "#9467bd",
  "CA212_Colitis_Ulcerative"                     = "#8c564b",
  "CA265_Diabetes_Mellitus_Type_1"               = "#e377c2",
  "CA391_Arthritis_Rheumatoid"                   = "#7f7f7f",
  "GCST003740_Barrett Esophagus"                 = "#bcbd22",
  "GCST90019016_Psoriasis"                       = "#17becf",
  "GCST90086158_Brugada Syndrome"                = "#4e79a7",
  "GD00307_Latent_Autoimmune_Diabetes_in_Adults" = "#f28e2b",
  "GD00398_Diverticular_Diseases"                = "#59a14f",
  "GD09058_Asthma"                               = "#e15759",
  "GD09159_Dermatitis_Atopic"                    = "#b07aa1",
  "GD09414_Vitiligo"                             = "#9c755f",
  "GD09476_Asthma"                               = "#76b7b2",
  "GD09493_Asthma"                               = "#edc948",
  "GD09629_Arthritis_Juvenile"                   = "#bab0ab",
  "GD09639_Lupus_Erythematosus_Systemic"         = "#af7aa1",
  "GD09657_Diabetes_Mellitus_Type_1"             = "#ff9da7",
  "MR04390_Eczema"                               = "#86bcb6",
  "PH378_Arthritis_Rheumatoid"                   = "#c5b0d5"
)


top_traits <- names(trait_cols)

plots <- lapply(unique(top_df$CellType), function(ct) {
  df_ct <- top_df %>% filter(CellType == ct)
  ggplot(df_ct, aes(x = score, y = reorder(TargetGene, score), color = traits)) +
    geom_point(size = 4, alpha = 0.9) +
    scale_color_manual(values = trait_cols, limits = top_traits, drop = FALSE) +
    labs(title = ct, x = "Score", y = "Gene") +
    theme_minimal(base_size = 12) +
    theme(
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6),
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_blank(),
      axis.ticks = element_line(color = "black"),
      axis.ticks.length = unit(2, "pt"),
      axis.text = element_text(color = "black"),
      axis.title.y = element_blank(),
      plot.title = element_text(face = "bold", hjust = 0.5)
    )
})

final_plot <- ggpubr::ggarrange(
  plotlist = plots,
  ncol = 2,
  nrow = ceiling(length(plots) / 2),
  common.legend = TRUE,
  legend = "right"
)
ggsave(
  "top_genes_by_celltype.pdf",
  final_plot,
  width = 12,
  height = 8
)


# add TF information
SECNIC <- read_tsv('eRegulon_direct.tsv') 
SECNIC_variants <- top_df %>%
  inner_join(SECNIC, by = c("TargetGene" = "Gene"), relationship = "many-to-many")
tf_summary <- SECNIC_variants %>%
  group_by(TF, TargetGene) %>%
  summarise(
    importance = mean(importance_R2G, na.rm = TRUE),
    rho = mean(rho_R2G, na.rm = TRUE),
    TF2G = mean(importance_TF2G, na.rm = TRUE),
    triplet = mean(triplet_rank, na.rm = T),
    rho_TF2G = mean(rho_TF2G, na.rm = T),
    .groups = "drop"
  ) %>%
  arrange(desc(importance))
tf_summary <- tf_summary %>%
  mutate(
    z_R2G = scale(importance)[,1],
    z_TF2G = scale(TF2G)[,1],
    z_rho = scale(rho)[,1],
    z_triplet = -scale(triplet)[,1],  
    combined_score = z_R2G + z_TF2G + z_rho + z_triplet
  )
top_tf <- tf_summary %>%
  group_by(TargetGene) %>%
  slice_max(order_by = z_triplet, n = 1, with_ties = FALSE) %>%
  ungroup()

library(circlize)
library(dplyr)

plot_df <- top_tf %>%
  mutate(weight = z_triplet) %>%
  select(TF, TargetGene, weight)
tf_names <- unique(plot_df$TF)
tf_palette <- colorRampPalette(c("#2b6cb0", "#c53030", "#2f855a", "#805ad5", "#d69e2e", "#dd6b20"))(length(tf_names))
tf_cols <- setNames(tf_palette, tf_names)
link_cols <- tf_cols[plot_df$TF]
gene_names <- unique(plot_df$TargetGene)
grid_cols <- c(
  setNames(tf_cols, tf_names),
  setNames(rep("#d1d5db", length(gene_names)), gene_names)
)

dev.new(width = 9, height = 9)
pdf("tf_gene_chord.pdf", width = 9, height = 9)

circos.clear()
circos.par(start.degree = 90, gap.degree = 2, track.margin = c(0.01, 0.01))

chordDiagram(
  x = plot_df,
  col = link_cols,
  grid.col = grid_cols,
  transparency = 0.2,
  annotationTrack = "grid",
  preAllocateTracks = list(track.height = 0.15)
)

circos.trackPlotRegion(
  track.index = 1,
  panel.fun = function(x, y) {
    sector.name = get.cell.meta.data("sector.index")
    xlim = get.cell.meta.data("xlim")
    ylim = get.cell.meta.data("ylim")
    circos.text(mean(xlim), ylim[1] + 0.1, sector.name,
                cex = 0.7, facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5))
  },
  bg.border = NA
)

dev.off()

