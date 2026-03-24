# Figure 3 - Pando GRN
library(Pando)
library(Seurat)
library(BSgenome.Hsapiens.UCSC.hg38)
library(Signac)
library(dplyr)
library(tidyr)
library(purrr)
library(doParallel)
library(plyr)
library(foreach)
library(stringr)
library(readr)
library(ggplot2)
library(tidygraph)
library(ggraph)
library(RColorBrewer)
set.seed(1234)

# load data
DC <- qs::qread("integrated.multiomics.chromvar.clustered.qsave")
DC <- FindVariableFeatures(DC, assay='RNA', nfeatures=2000)
DC <- RegionStats(DC, genome = BSgenome.Hsapiens.UCSC.hg38,assay = 'ATAC')

# link peaks to genes
DC <- LinkPeaks(
  object = DC,
  peak.assay = "peaks",
  expression.assay = "RNA",
  genes.use = VariableFeatures(DC, assay='RNA')
)
DefaultAssay(DC) <- "ATAC"

# Initiate GRN object and select candidate regions
DC = initiate_grn(
  DC, 
  peak_assay = 'ATAC', 
  rna_assay = 'RNA', 
  exclude_exons = F
)

DefaultAssay(DC) <- "ATAC"

DC = infer_grn(
  DC,
  peak_to_gene_method = "Signac",
  only_tss = T,
  parallel = F,
  tf_cor = 0.1,
  aggregate_peaks_col = "group",
  method = "glm", 
  verbose = 2)

DC <- find_modules(DC)
modules <- NetworkModules(DC)
modules@meta

glm_coefs <- coef(DC2, network='glm_network')

GetGRN(DC)
coef(DC)

# extract glm network
grn_net <- glm_coefs %>% 
  filter(padj<0.05) %>% 
  filter(!str_detect(target, '^(AMEX|LOC)')) %>% 
  group_by(tf, target) %>% 
  filter(pval==min(pval)) %>% 
  summarize(
    estimate=mean(estimate),
    pval=min(pval),
    padj=min(padj),
    regions=region[1]
  )
grn_net %>% write_tsv('glm_modules.tsv')


# DAR
da_peaks <- FindMarkers(
  object = DC2,
  ident.1 = 'mDC',
  ident.2 = NULL,
  test.use = 'LR',
  latent.vars = 'nCount_peaks'
)
#write.csv(da_peaks, file = "mDC.csv", row.names = T, col.names = T, quote = F)
cDC1 = read.csv(file = "cDC1.csv",row.names =1) %>% mutate(cluster = "cDC1")  %>% filter(., p_val_adj < 0.05 & avg_log2FC > 0) 
pDC = read.csv(file = "pDC.csv",row.names =1) %>% mutate(cluster = "pDC") %>% filter(., p_val_adj < 0.05 & avg_log2FC > 0) 
mDC = read.csv(file = "mDC.csv",row.names =1) %>% mutate(cluster = "mregDC")  %>% filter(., p_val_adj < 0.05 & avg_log2FC > 0) 
cDC2 = read.csv(file = "cDC2.csv",row.names =1) %>% mutate(cluster = "cDC2")   %>% filter(., p_val_adj < 0.05 & avg_log2FC > 0) 
l = list(A=cDC1, B=pDC, C=mDC, D= cDC2)
df <- data.table::rbindlist(l)
df[,Col := unlist(lapply(l, rownames))]
df <- df %>% dplyr::select(Col, everything())

cluster_stats <- unique(df$cluster) %>%
  map_dfr(~ {
    grn_net %>%
      dplyr::filter(regions %in% df$Col[df$cluster == .x]) %>%
      dplyr::group_by(tf) %>%
      dplyr::summarize(ngenes = n()) %>%
      mutate(tf = str_to_title(tf),
             tf = factor(tf, levels = unique(tf)),
             cluster = .x)
  }) %>%
  dplyr::select(cluster, tf, ngenes) %>%
  pivot_wider(names_from = cluster, values_from = ngenes, values_fill = 0)

expression_long <- cluster_stats %>% pivot_longer(cols = cDC1:cDC2,names_to = "cell_type",values_to = "connections")
expression_long$tf = toupper(expression_long$tf)


# add AUC chromvar pvalues

DA_motifs_ct <- presto::wilcoxauc(DC, group_by = "group", seurat_assay = "chromvar")
DA_motifs_ct$gene <- ConvertMotifID(DC, id = DA_motifs_ct$feature, assay = "peaks")
DA_motifs_ct$gene = toupper(DA_motifs_ct$gene)
DA_motifs_ct <- DA_motifs_ct %>% rename("cell_type"=group, "tf" = gene)
DA_motifs_ct <- DA_motifs_ct %>% filter(padj < 0.01) %>% group_by(cell_type)

DAR = subset(da_peaks, p_val_adj < 0.05 & avg_log2FC > 0)


# merge chromvar and pando
combined_long <- expression_long %>%
  left_join(DA_motifs_ct, by = c("cell_type", "tf"))
combined_long <- combined_long[!(is.na(combined_long$feature)), ]

top_tfs <- combined_long %>%
  group_by(cell_type) %>%                           
  slice_max(order_by = connections,                   
            n = 10,                                  
            with_ties = TRUE) %>% distinct(tf)    

genes.list = c("cDC1","mDC" ,"pDC","cDC2")

plot.list <- list()
for (gene in genes.list) {
  top_tfs <- combined_long %>%
    filter(cell_type == gene) %>%
    arrange(desc(connections)) %>%
    distinct(tf, .keep_all = TRUE) %>%
    slice_head(n = 10) %>%
    select(tf, cell_type, connections, auc)
  plot <- ggplot(top_tfs, aes(x = reorder(tf, connections), y = connections)) +
    geom_segment(aes(xend = tf, y = 0, yend = connections), color = "grey50") +
    geom_point(aes(size = auc, color = cell_type)) +
    coord_flip(clip = "off") +
    scale_colour_manual(values = c("cDC1" = "#96bd95", "mregDC" = "#c288ab", "pDC" = "#8dc0e3", "cDC2" = "#e3d58d")) +
    theme_bw() +
    theme(
      axis.title = element_blank(),
      panel.grid.major.y = element_blank(),
      legend.position = "bottom",
      axis.ticks.y = element_blank())
  plot.list[[gene]] <- plot
}
ggpubr::ggarrange(plotlist = plot.list, ncol=1, nrow=4,common.legend = T)  %>% 
  ggpubr::ggexport(filename = "net_stats_clusters.pdf", width = 3, height = 8.5)

#### Plot network ####

# add latent time
DC <- NormalizeData(DC)
RNAdata <-FetchData(DC,vars = c("latent_time","barcode","velocity_pseudotime","group",tf), slot = 'data')
RNA_melt <- data.table::melt(setDT(RNAdata), id.vars = c("latent_time","velocity_pseudotime","group","barcode"),
                             variable.name = "gene", value.name = "expression") %>% tibble::as_tibble()
process_gene <- function(subsetData) {
  fit_expression <- mgcv::gam(expression ~ s(latent_time, k = 6, bs = "cr", sp = 10),family = gaussian, data = subsetData)
  subsetData$expression_fit <- fit_expression$fitted.values
  return(subsetData)
}
data.plot <- RNA_melt %>% group_by(gene) %>% do(process_gene(.)) %>% as_tibble()
corResults <- setDT(data.plot)[,
                               .(
                                 r_expression = cor(latent_time, expression_fit, method = "pearson"),
                                 p_expression = cor.test(latent_time, expression_fit)$p.value
                               ),
                               by = gene
]

reg_mat <- grn_net %>% 
  filter(!str_detect(target, '^(AMEX|LOC)')) %>% 
  distinct(target, tf, estimate) %>%
  pivot_wider(names_from=tf, values_from=estimate, values_fill=0) %>% 
  column_to_rownames('target') %>% as.matrix() %>% Matrix::Matrix(sparse=T)
reg_factor_mat <- abs(reg_mat) + 1
gene_cor_subset <- as.matrix(gene_cor[rownames(reg_factor_mat), colnames(reg_factor_mat)])
reg_factor_mat_dense <- as.matrix(sqrt(reg_factor_mat))
weighted_coex_mat <- gene_cor_subset * reg_factor_mat_dense
weight_coex_umap <- uwot::umap(weighted_coex_mat, n_neighbors=5)
rownames(weight_coex_umap) <- rownames(weighted_coex_mat)
colnames(weight_coex_umap) <- c('UMAP1', 'UMAP2')
weight_coex_meta <- weight_coex_umap %>% 
  as_tibble(rownames='gene')  %>%
  left_join(corResults)


glm_graph <- as_tbl_graph(grn_net) %>% 
  activate(edges) %>% 
  mutate(from_node=.N()$name[from], to_node=.N()$name[to]) %>% 
  activate(nodes) %>% 
  mutate(
    central_betw=centrality_betweenness(),
    central_eig=centrality_eigen(),
    central_deg=centrality_degree(),
    outdegree=centrality_degree(mode='out'),
    indegree=centrality_degree(mode='in')
  ) %>% 
  inner_join(weight_coex_meta, by=c('name'='gene')) %>% 
  activate(edges) %>%
  filter(padj<1e-2)

grn_graph <- glm_graph %>%
  activate(nodes) %>%  
  filter(name %in% grn_net$tf)
pal = c("#5fbaad","#8cded2","#b5ebe3","#dff5f2","#eb9bbb","#e36b98","#c7386f","#ba205b")

library(rcartocolor)
p=ggraph(grn_graph, x=UMAP1, y=UMAP2) + 
  geom_edge_diagonal(aes(alpha=-log10(padj), 
                         color=factor(sign(estimate))), width=0.5,
                     #arrow = arrow(length = unit(4, 'mm'), type = 'closed'))+
                     geom_node_point(aes(fill=r_expression,size = outdegree), shape=21, color='black',position = position_jitter(width = 0.1, height = 0.2)) +
                       geom_node_text(aes(label=name), size=8/ggplot2::.pt, repel=T, colour = "grey30") +
                       scale_edge_color_manual(name = "Interactions", values=c('#a3c7d1', '#ebdece'),
                                               labels = c("Negative","Positive"),
                                               guide = ggplot2::guide_legend(
                                                 title.position = "top")) +
                       scale_edge_alpha_continuous(name = "Strength\n(p-value)",range=c(0.01,0.8), limits=c(2,20),guide = ggplot2::guide_legend(
                         title.position = "top")) +
                       scale_fill_gradientn(colours = pal, limits = c(-1, 1),name = "Pseudotime\nCorrelation",
                                            oob = scales::squish, breaks = c(-1, 0, 1),
                                            guide = guide_colorbar(barheight = 4, barwidth = 1,
                                                                   ticks.colour = "black",
                                                                   title.position = 'top',
                                                                   frame.colour = "black",
                                                                   title.vjust = 0)) +
                       scale_size_binned(name = "Centrality")+
                       ggplot2::theme(
                         panel.background = ggplot2::element_rect(fill = NA, colour = NA),
                         plot.background = ggplot2::element_rect(fill = NA, colour = NA),
                         plot.title = ggplot2::element_blank(),
                         plot.subtitle = ggplot2::element_blank(),
                         legend.position="top",
                         legend.direction = "vertical",
                         legend.justification = "left",
                         plot.margin = ggplot2::unit(c(0, 0, 0, 0), "cm"),
                         axis.title.x = ggplot2::element_blank(),
                         axis.title.y = ggplot2::element_blank(),
                         axis.text.x = ggplot2::element_blank(),
                         axis.text.y = ggplot2::element_blank(),
                         axis.ticks.x = ggplot2::element_blank(),
                         axis.ticks.y = ggplot2::element_blank(),
                         panel.grid.major = ggplot2::element_blank(),
                         panel.grid.minor = ggplot2::element_blank())
p





