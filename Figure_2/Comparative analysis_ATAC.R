# Comparative analysis (ATAC)
library(pheatmap)
library(Seurat)
library(Signac)
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)
library(ComplexHeatmap)
library(circlize)
library(grid)
library(Matrix)
# inputs
cb = qs::qread("dc_culture.qsave")
adult = qs2::qs_read('panDC.qs2') # similar for BM

cb_counts <- FeatureMatrix(
    fragments = Fragments(cb[["ATAC"]]),
    features = combined_peaks,
    cells = colnames(cb)
)
  
adult_counts <- FeatureMatrix(
    fragments = Fragments(adult[["peaks"]]),
    features = combined_peaks,
    cells = colnames(adult)
)

pseudobulk <- function(counts, groups) {
  sapply(split(colnames(counts), groups), function(cells) Matrix::rowSums(counts[, cells, drop = FALSE]))
}

cb_pb    <- pseudobulk(cb_counts, DC$group)
adult_pb <- pseudobulk(adult_counts, adult$celltype)

# CPM-normalize
cb_cpm    <- sweep(cb_pb, 2, colSums(cb_pb), "/") * 1e6
adult_cpm <- sweep(adult_pb, 2, colSums(adult_pb), "/") * 1e6
cor_mat <- cor(log1p(cb_cpm), log1p(adult_cpm), method = "spearman")

cor_z <- scale(cor_mat)  # z-scores each column
p1 <- ComplexHeatmap::Heatmap(
  cor_z,
  width  = ncol(cor_z) * unit(8, "mm"),
  height = nrow(cor_z) * unit(6, "mm"),
  name = "z-score",
  col = circlize::colorRamp2(c(min(cor_z), 0, max(cor_z)), c("#2166ac", "white", "#7d191f")),
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  column_names_side = "top",
  row_names_side = "left",
  show_column_dend = FALSE,
  show_row_dend = FALSE,
  border = "black",
  row_names_gp = gpar(fontsize = 9),
  column_names_gp = gpar(fontsize = 9),
  heatmap_legend_param = list(
    title_gp  = gpar(fontsize = 8, fontface = "bold"),
    labels_gp = gpar(fontsize = 8),
    border    = "black"
  )
)


# DAR
### DAR
cb = qs::qread("dc_culture.qsave")
adult = qs2::qs_read('panDC.qs2')

DefaultAssay(cb) = "peaks"
DefaultAssay(adult) = "peaks"

da_peaks <- FindMarkers(
  object = adult,
  ident.1 = 'pDC', #cDC1/pDC/cDC2
  ident.2 = NULL,
  test.use = 'LR',
  latent.vars = 'nCount_peaks'
)
write.csv(da_peaks, 'DARs_panDC.csv')

DAR_panDC <- read.csv('panDC.dar.csv')
adult_dar_list <- setNames(lapply(celltypes, function(ct) {
  DAR_panDC %>% filter(p_val_adj < 0.05, avg_log2FC > 0.25, celltype == ct)
}), celltypes)
culture_paths <- c(
  pDC  = "pDC.csv",
  cDC1 = "cDC1.csv",
  cDC2 = "cDC2.csv"
)
celltypes <- c("pDC", "cDC1", "cDC2")
culture_da_list <- setNames(lapply(celltypes, function(ct) {
  read.csv(culture_paths[[ct]], row.names = 1) %>%
    filter(p_val_adj < 0.05, avg_log2FC > 0)
}), celltypes)
adult_gr_list   <- setNames(lapply(celltypes, function(ct) StringToGRanges(adult_dar_list[[ct]]$peak)), celltypes)
culture_gr_list <- setNames(lapply(celltypes, function(ct) StringToGRanges(rownames(culture_da_list[[ct]]))), celltypes)
results <- lapply(celltypes, function(ct) {
  adult_gr   <- adult_gr_list[[ct]]
  culture_gr <- culture_gr_list[[ct]]
  
  hits_a_in_c <- findOverlaps(adult_gr, culture_gr)
  pct_adult_in_culture <- length(unique(queryHits(hits_a_in_c))) / length(adult_gr) * 100
  
  hits_c_in_a <- findOverlaps(culture_gr, adult_gr)
  pct_culture_in_adult <- length(unique(queryHits(hits_c_in_a))) / length(culture_gr) * 100
  
  ov_df <- as.data.frame(hits_a_in_c)
  ov_df$adult_lfc   <- adult_dar_list[[ct]]$avg_log2FC[ov_df$queryHits]
  ov_df$culture_lfc <- culture_da_list[[ct]]$avg_log2FC[ov_df$subjectHits]
  ov_df$celltype <- ct
  effect_cor <- suppressWarnings(cor(ov_df$adult_lfc, ov_df$culture_lfc, method = "spearman"))
  
  list(
    summary = data.frame(
      celltype = ct,
      n_adult_peaks = length(adult_gr),
      n_culture_peaks = length(culture_gr),
      pct_adult_in_culture = round(pct_adult_in_culture, 1),
      pct_culture_in_adult = round(pct_culture_in_adult, 1),
      effect_size_cor = round(effect_cor, 3)
    ),
    overlap_df = ov_df
  )
})
names(results) <- celltypes 

plot_df <- do.call(rbind, lapply(results, `[[`, "overlap_df"))

ggplot(plot_df, aes(x = adult_lfc, y = culture_lfc)) +
  geom_point(alpha = 0.4, size = 1, color = "#7d191f") +
  geom_smooth(method = "lm", color = "black", se = FALSE, linewidth = 0.5) +
  stat_cor(method = "spearman", size = 4, label.x.npc = "left", label.y.npc = "top") +
  facet_wrap(~celltype, scales = "free") +
  labs(x = "Adult log2FC (in vivo)", y = "Cultured log2FC (cord blood)") +
  theme_minimal()+
  theme_minimal() +
  theme(
    legend.position = "right",
    panel.grid = element_blank(), # Remove grids for the clean "line-to-axis" look
    axis.line.y = element_line(color = "black"),
    axis.line.x = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    axis.text = element_text(color = "black"),
    axis.title = element_text(face = "bold"),
    plot.title = element_text(face = "bold", size = 14)
  )
dev.off()


df_stacked_list <- lapply(celltypes, function(ct) {
  adult_gr   <- adult_gr_list[[ct]]
  culture_gr <- culture_gr_list[[ct]]
  
  shared_adult <- length(unique(queryHits(findOverlaps(adult_gr, culture_gr))))
  unique_adult <- length(adult_gr) - shared_adult
  
  data.frame(
    CellType = ct,
    Shared   = shared_adult,
    Unique   = unique_adult
  )
})

df_stacked <- bind_rows(df_stacked_list) %>%
  pivot_longer(cols = c(Shared, Unique), names_to = "group", values_to = "total") %>%
  group_by(CellType) %>%
  mutate(pct = total / sum(total) * 100) %>%   # correct, per-cell-type percentage
  ungroup()

df_stacked$group <- factor(df_stacked$group, levels = c("Unique", "Shared"))
df_stacked$CellType <- factor(df_stacked$CellType, levels = celltypes)

p <- ggplot(df_stacked, aes(x = CellType, y = total, fill = group)) +
  geom_bar(stat = "identity", position = "fill", color = "black", width = 0.6) +
  labs(x = "", y = "Proportion of in vivo DA Peaks") +
  theme_minimal() +
  scale_fill_manual("", values = c("Unique" = "#E87D72", "Shared" = "grey60")) +
  geom_text(
    aes(label = ifelse(group == "Shared", percent(pct / 100, accuracy = 0.1), "")),
    position = position_fill(vjust = 0.5),
    size = 5, color = "white"
  ) +
  theme(
    axis.text.x = element_text(size = 8, color = "black"),
    axis.line.x = element_line(linewidth = 0.4, color = "black"),
    axis.line.y = element_line(linewidth = 0.4, color = "black"),
    axis.text.y = element_text(size = 8, color = "black"),
    axis.ticks = element_line(linewidth = 0.25, color = 'black'),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 8, color = "black"),
    legend.key.size = unit(0.8, 'cm'),
    legend.text = element_text(size = 8),
    legend.key = element_rect(fill = alpha("white", .5), color = "black"),
    panel.background = element_blank(),
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.background = element_blank()
  )+ scale_y_continuous(
    expand = c(0, 0))
p


# bulk ATAC
bam_dir   <- "pDC/GSE279911/"
macs2_dir <- "pDC/GSE279911/PEAKS/"

adult_bam_files  <- list.files(bam_dir, pattern = "\\.bam$", full.names = TRUE)
adult_peak_files <- list.files(macs2_dir, pattern = "\\.narrowPeak$", full.names = TRUE)

# 500-BP CONSENSUS PEAK SET
adult_peaks_list <- lapply(adult_peak_files, import)
adult_peaks      <- unlist(GRangesList(adult_peaks_list))
cb_peaks <- granges(cb[["peaks"]])
all_peaks       <- c(cb_peaks, adult_peaks)
peak_centers    <- resize(all_peaks, width = 1, fix = "center")
consensus_peaks <- reduce(resize(peak_centers, width = 500, fix = "center"))
consensus_peaks <- keepStandardChromosomes(consensus_peaks, pruning.mode = "coarse")

# QC
library(Rsamtools)
bam_chr_names   <- names(scanBamHeader(adult_bam_files[1])[[1]]$targets)
uses_chr_prefix <- any(grepl("^chr", bam_chr_names))
chr_vec <- as.character(seqnames(consensus_peaks))
if (!uses_chr_prefix) chr_vec <- gsub("^chr", "", chr_vec)
saf <- data.frame(
  GeneID = GRangesToString(consensus_peaks),
  Chr    = chr_vec,
  Start  = start(consensus_peaks),
  End    = end(consensus_peaks),
  Strand = "*"
)
fc_adult <- featureCounts(
  files = adult_bam_files,
  annot.ext = saf,
  isPairedEnd = TRUE,
  nthreads = 1
)
adult_counts <- fc_adult$counts
rownames(adult_counts) <- GRangesToString(consensus_peaks)
colnames(adult_counts) <- gsub("\\.bam$", "", basename(colnames(adult_counts)))
assigned    <- fc_adult$stat[fc_adult$stat$Status == "Assigned", -1]
no_features <- fc_adult$stat[fc_adult$stat$Status == "Unassigned_NoFeatures", -1]
frip_scores <- as.numeric(assigned / (assigned + no_features))
names(frip_scores) <- colnames(adult_counts)
pass_qc_adult <- names(frip_scores[frip_scores >= 0.20])

# map
cb_pseudo_bulk <- AggregateExpression(
  cb, assays = "peaks", group.by = "group", slot = "counts"
)$peaks
rownames(cb_pseudo_bulk) <- GRangesToString(cb_peaks)
overlaps <- findOverlaps(consensus_peaks, cb_peaks)

mapping_matrix <- sparseMatrix(
  i = queryHits(overlaps),
  j = subjectHits(overlaps),
  x = 1,
  dims = c(length(consensus_peaks), nrow(cb_pseudo_bulk))
)
cb_consensus_counts <- as.matrix(mapping_matrix %*% cb_pseudo_bulk)
rownames(cb_consensus_counts) <- GRangesToString(consensus_peaks)
colnames(cb_consensus_counts) <- colnames(cb_pseudo_bulk)
adult_qc_counts     <- adult_counts[, pass_qc_adult, drop = FALSE]
adult_pooled_counts <- rowSums(adult_qc_counts)
adult_label <- paste0("Adult Target (N=", length(pass_qc_adult), ")")

combined_counts <- cbind(
  adult_pooled_counts,
  cb_consensus_counts
)
colnames(combined_counts)[1] <- adult_label

# Calculate CPM, log2 transform, and apply Quantile Normalization
cpm_matrix   <- cpm(combined_counts)
log_cpm      <- log2(cpm_matrix + 1)
log_cpm_norm <- normalizeBetweenArrays(log_cpm, method = "quantile")

peak_vars   <- apply(log_cpm_norm, 1, var)
valid_peaks <- peak_vars > 0 & !is.na(peak_vars)
top_20k_idx <- head(order(peak_vars[valid_peaks], decreasing = TRUE), 20000)
global_matrix <- log_cpm_norm[valid_peaks, ][top_20k_idx, ]

# Compute global Spearman correlation matrix
global_cor_matrix <- cor(global_matrix, method = "spearman")
qs2::qs_save(
  global_cor_matrix, 
  'cor_pDC.qs2'
)


###
cor_pDC  <- qs2::qs_read('cor_pDC.qs2')
cor_cDC1 <- qs2::qs_read('cor_cDC1.qs2')
cor_cDC2 <- qs2::qs_read('cor_cDC2.qs2')
adult_row <- function(mat) {
  adult_col_name <- grep("^Adult Target", colnames(mat), value = TRUE)
  mat[adult_col_name, setdiff(colnames(mat), adult_col_name)]
}
cor_matrix_combined <- rbind(
  pDC  = adult_row(cor_pDC),
  cDC1 = adult_row(cor_cDC1),
  cDC2 = adult_row(cor_cDC2)
)
celltypes <- c("pDC", "cDC1", "cDC2")
for (ct in celltypes) {
  matched <- cor_matrix_combined[ct, ct]
  other <- mean(cor_matrix_combined[ct, setdiff(colnames(cor_matrix_combined), ct)])
  cat(ct, ": matched =", round(matched, 3), " vs other =", round(other, 3), "\n")
}
col_order <- c("pDC", "cDC1", "cDC2", "mDC", "HSPC", "Lineage")
row_order <- c("pDC", "cDC1", "cDC2")

mat <- cor_matrix_combined[row_order, col_order]
max_abs <- max(abs(mat))
col_fun <- colorRamp2(c(-max_abs, 0, max_abs), c("#2166ac", "white", "#b2182b"))
matched_pairs <- data.frame(row = row_order, col = row_order)  # same name = matched

ht <- Heatmap(
  mat,
  name = "Spearman r",
  col = col_fun,
  
  rect_gp = gpar(col = "white", lwd = 1.2),
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  column_names_side = "top",
  row_names_side = "left",
  column_title = "Cord-blood-derived cell type",
  row_title = "Adult bulk ATAC",
  
  cell_fun = function(j, i, x, y, width, height, fill) {
    r_val <- mat[i, j]
    is_matched <- colnames(mat)[j] == rownames(mat)[i]
    
    grid.text(
      sprintf("%.2f", r_val), x, y,
      gp = gpar(fontsize = 9,
                fontface = ifelse(is_matched, "bold", "plain"),
                col = ifelse(abs(r_val) > max_abs * 0.6, "white", "black"))
    )
    
    if (is_matched) {
      grid.rect(x, y, width = width, height = height,
                gp = gpar(col = "black", lwd = 2.5, fill = NA))
    }
  },
  
  row_names_gp = gpar(fontsize = 11, fontface = "bold"),
  column_names_gp = gpar(fontsize = 11, fontface = "bold"),
  heatmap_legend_param = list(
    title_gp  = gpar(fontsize = 9, fontface = "bold"),
    labels_gp = gpar(fontsize = 9),
    border    = "black"
  ),
  border = "black"
)

draw(ht)
