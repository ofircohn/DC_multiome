# Figure 6 - TRS and scE2G
### box plot
library(ggthemes)
library(dplyr)
library(tidyr)
library(tibble)
library(GenomicRanges)
library(data.table)
library(readxl)
library(scales)
library(ggplot2)
set.seed(1234)

df2 = read_excel('path/to/TRS/results')
df2_filter = df2 %>% filter(df2$id == "GD09639_Lupus_Erythematosus_Systemic")
df2_filter$category = substr(df2_filter$cell_cluster, 1, 4)
df2_filter = filter(df2_filter, category %in% c("cDC1","cDC2","pDC_")) 
df2_filter <- df2_filter %>%
  mutate(category = fct_reorder(category, TRS, .fun = median, .desc = TRUE))
p = ggplot(df2_filter, aes(x=category, y=TRS)) + 
  geom_boxplot(aes(fill = category, alpha = 0.2), colour = "grey50", show.legend = FALSE, 
               outlier.shape = NA, alpha = 0.2) +
  xlab("")+
  scale_fill_brewer(palette = "Paired", direction = 1, guide = FALSE) +
  ggthemes::theme_tufte(base_family = "Helvetica", ticks = TRUE) +
  theme(panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = ggplot2::element_blank(),
        axis.line = ggplot2::element_line(colour = "black"),
        axis.text = element_text(size = 14, colour = "black"),
        axis.title.x = element_text(margin = margin(t = 10), hjust = 0.5))
p
p + scale_fill_viridis_d(option = "plasma", direction = -1) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
p 
ggsave(plot= p, "boxplot_TRS.pdf", width=3, height=3, bg = "transparent")

# scE2G: variant to gene
snps = read.table("/Users/cohn/DendriticProject/FileMappedData/causaldb/BED_HLA_hg38/GD09639_Lupus_Erythematosus_Systemic.bed", 
                  col.names = c("chrom", "start", "end","region1", "PP","rs")) 
snps$end =  snps$start 
snps = snps%>% GenomicRanges::makeGRangesFromDataFrame(.,keep.extra.columns = T)

scE2G_pDC = fread('/Users/cohn/DendriticProject/FileMappedData/caQTL/ABC_results/scE2G/pDC/encode_e2g_predictions_threshold0.177.tsv.gz')
scE2G_pDC.gr = scE2G_pDC %>% dplyr::select(chr, start, end, TargetGene, ABC.Score, CellType, Kendall, E2G.Score.qnorm) %>%
  makeGRangesFromDataFrame(., keep.extra.columns = TRUE) 
scE2G_pDC.gr <- plyranges::join_overlap_left(snps, scE2G_pDC.gr) 

plot_df <- scE2G_pDC.gr %>% as_tibble() %>%
  dplyr::filter(!is.na(TargetGene), !is.na(rs), !is.na(E2G.Score.qnorm)) %>%
  dplyr::distinct(rs, TargetGene, E2G.Score.qnorm)
plot_df <- plot_df %>%
  arrange(TargetGene, rs) %>%
  mutate(rs = factor(rs, levels = unique(rs)))




p1 <- ggplot(plot_df, aes(x = rs, y = TargetGene, fill = E2G.Score.qnorm)) +
  geom_tile(color = "black", size = 0.3) +   
  scale_fill_gradientn(
    colors = c("gray95", "white", "#7d1646"),  
    guide = guide_colorbar(
      ticks.colour = "black",
      frame.colour = "black",
      barwidth = 1.3
    ),
    oob = squish
  ) +
  theme_minimal(base_size = 10) +
  labs(
    x = "SNPs",
    y = ""  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    legend.key = element_rect(fill = "white", color = "black"),  # black border around color bar
    legend.background = element_rect(color = NA, fill = NA),
    legend.title = element_text(color = "black"),
    legend.text = element_text(color = "black"),
    axis.text.y = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    panel.grid = element_blank(),
    axis.ticks = element_line(color = "black")  # add ticks
  )

p1