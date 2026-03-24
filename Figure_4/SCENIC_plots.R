# Figure 4 - SECNIC+ plots

library(data.table)
library(dplyr)
library(tibble)
library(scattermore)
library(tidyverse)

# load eRegulon specificity score (eRSS)
rss = fread('rss.csv') %>% as_tibble() %>% column_to_rownames('V1')
rss_values <- sweep(rss,2,colSums(rss),`/`)*100
rss_values <- rss_values[,sort(colnames(rss_values))]
rss_values <- t(apply(rss_values, 1, function(x)(x-min(x))/(max(x)-min(x))))

RSS_long <- rss_values %>% pivot_longer(cols = MEP:Lineage,names_to = "cell_type",values_to = "RSS")
exp = fread('TF_expression_gene_based.csv') %>% as_tibble() %>% column_to_rownames('group') %>% t() %>% as.data.frame()
exp = t(scale(t((exp))))  %>% as.data.frame() %>% rownames_to_column('gene') %>% filter(gene %in% rss_values$gene)
expression_long <- exp %>% pivot_longer(cols = pDC:Lineage,names_to = "cell_type",values_to = "expression")
scale_min_max <- function(x) {
  (x - min(x)) / (max(x) - min(x))
}
expression_long$expression_scaled = scale_min_max(expression_long$expression)
merged = left_join(RSS_long, expression_long, by = c("gene","cell_type"))
merged$regulon_gene = paste0(merged$gene, "_", merged$regulon)


# cDC1 negative regulators
cDC1 <- merged %>%
  filter(cell_type == "cDC1" & str_detect(regulon_gene, "direct_-"))
corResults = read.csv("corResults.csv",header = T, row.names = 1) # pseudotime
pseudotime_regulators_cDC1 = cDC1 %>% left_join(., corResults, by= "gene" ) %>% na.omit()

# pDC positive regulators
pDC <- merged %>%
  filter(cell_type == "pDC" & str_detect(regulon_gene, "direct_\\+/\\+_"))
corResults = read.csv("/Users/cohn/Fellowships/Lupus/application/plots/corResults.csv",header = T, row.names = 1) #pseudotime
pseudotime_regulators_pDC = pDC %>% left_join(., corResults, by= "gene" ) %>% na.omit()

pseudotime_regulators_cDC1 <- pseudotime_regulators_cDC1 %>%
  mutate(regulon_gene_ordered = fct_reorder(regulon_gene, r_expression, .na_rm = TRUE))
pal = c("#5fbaad","#8cded2","#b5ebe3","#dff5f2","#eb9bbb","#e36b98","#c7386f","#ba205b")
p=ggplot(pseudotime_regulators_cDC1,aes(fct_reorder(regulon_gene, r_expression), r_expression, fill = r_expression)) +
  scale_y_continuous(expand = expansion(mult = c(0.01, 0.05))) +  # Add a small padding
  ylab("cDC1\n Negative regulators") +  xlab("") +
  scale_fill_gradientn(colours =pal , limits = c(-1, 1),name = "Pseudotime\nCorrelation",
                       oob = scales::squish, breaks = c(-1, 0, 1),
                       guide = guide_colorbar(barheight = 4, barwidth = 1,
                                              ticks.colour = "black",
                                              title.position = 'top',
                                              frame.colour = "black",
                                              title.vjust = 0)) +
  geom_point(data = pseudotime_regulators_cDC1,pch = 21, aes(size = RSS))+
  geom_segment(aes(x = regulon_gene, xend = regulon_gene, y = 0, yend = r_expression)) +
  #scale_size(range = c(0, 10))+
  guides(fill = guide_colorbar(frame.colour = "black", ticks.colour = "black"))+
  ggplot2::theme(legend.position="bottom",
                 panel.grid = element_blank(),
                 strip.background = element_blank(),
                 panel.grid.major = element_blank(),
                 axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 12),
                 axis.title.y = element_text(size = 14),
                 panel.spacing = unit(0, "mm"),
                 panel.border = element_rect(colour = "black", fill = NA),
                 axis.ticks.length = unit(0.25, "cm"),
                 plot.title = ggplot2::element_text(size=14),
                 axis.line.y.left = ggplot2::element_line(),
                 plot.background = element_rect(fill = "transparent", color = NA),
                 panel.background = element_rect(fill = "transparent", color = NA),
                 axis.line.x.bottom = ggplot2::element_line())
p

## R version of scenicplus heatmap_dotplot

subset <- merged %>%
  filter(regulon_gene %in% c('ZEB2_direct_-/+_(17g)','TCF4_direct_+/+_(613g)'))
library(subset)
p=ggplot(subset2, aes(x = cell_type, y = regulon_gene, fill = expression))+
  ggplot2::geom_tile(color = "white", linewidth = 0.5) +
  ggplot2::scale_y_discrete(expand = c(0, 0)) +
  ggplot2::scale_x_discrete(expand = c(0, 0),position = "top",guide = guide_axis(angle = 90))+
  ggplot2::coord_equal() +
  geom_point(mapping = aes(size = RSS))+
  ylab('')+ xlab("")+
  scale_fill_gradient2(low ="#6baed6", high =  "#d6936b", midpoint = 0.001, na.value = "grey50",
                       guide = guide_colorbar(ticks.colour = "black",
                                              frame.colour = "black",
                                              barwidth=1.3))+
  theme(
    plot.background =  element_rect(fill = "white", color = "white"),
    panel.background = ggplot2::element_rect(fill = "white", color = "white"),
    panel.grid.major = ggplot2::element_blank(),
    panel.border = ggplot2::element_rect(fill = NA, color = "black", linewidth = 1),
    plot.margin = ggplot2::margin(t = 0, r = 0, b = 0, l = 0),
    panel.grid.minor.y = ggplot2::element_line(color = "white", linewidth = 1),
    legend.background = ggplot2::element_rect(fill = "white", color = "white"))
p



## R version of UMAP eRegulon enrichment scores
umap = read.csv('umap.csv', header = T, sep="\t") %>% as_tibble()
umap$group = as.factor(umap$group)
my_cols <- c('pDC'='#97bee8','cDC2'='#e0a09d','cDC1'='#F2CC8F', 'mDC'='#e3bcd4','MEP' = "#b3d1b0", 'DC3' ='#e9f0ad', 'Lineage' = '#a5cdcf' )
arr <- list(x = -3, y = -3, x_len = 2, y_len = 2)
ggplot(umap) + aes(x = UMAP_1, y = UMAP_2, color = group)+
  geom_scattermore(alpha = 0.95, pointsize = 3.5, pixels    = c(1000, 1000)) +
  labs(colour = "Celltype") +
  xlab("") +  
  ylab("") +
  scale_colour_manual(values = my_cols)+
  annotate("segment", 
           x = arr$x, xend = arr$x + c(arr$x_len, 0), 
           y = arr$y, yend = arr$y + c(0, arr$y_len), 
           arrow = arrow(type = "closed", length = unit(10, 'pt')))+
  theme_bw() +
  theme(
    legend.background = element_rect(color = NA),
    axis.title = element_blank(),           
    axis.text = element_blank(),           
    axis.ticks = element_blank(),          
    
    axis.line = element_blank(),          
    strip.placement = "outside",
    strip.background = element_rect(fill = "gray95", colour = "white"),
    panel.border = element_blank(),        
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.background = element_blank())
