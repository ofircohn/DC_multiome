# Figure 4 - TF expression/peaks-pseudotime correlations

library(plyr)
library(Signac)
library(GenomicRanges)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(ggplot2)
library(circlize)
library(Seurat)
library(dplyr)
library(tibble)
library(tidyr)
library(data.table)
library(purrr)
library(mgcv)
library(readxl)
library(Hmisc)
library(reticulate)
library(Pando)
library(BSgenome.Hsapiens.UCSC.hg38)
library(Matrix.utils)
library("scales")
library(scattermore)

setDTthreads(20)
doMC::registerDoMC(cores=30)

# load Seurat object
DC = qs::qread('MEP_cDC1.qsave')

# load fate/mata/pseudotimes
# cDC1
meta = read.csv("cDC1_meta.csv",header = T, row.names = 1)[,c('velocity_pseudotime','latent_time',"CellTypeAll")]
meta$barcode = rownames(meta)
# pDC
meta = read.csv("pDC_meta.csv",header = T, row.names = 1)[,c('velocity_pseudotime','latent_time',"CellTypeAll")]
meta$barcode = rownames(meta)

# add to metadata to object
meta <- meta[order(match(rownames(meta), rownames(DC@meta.data))), ]
DC <- AddMetaData(DC, meta[,c("barcode","latent_time","velocity_pseudotime")])
pseudotime = DC@meta.data["latent_time"]


# add pando motif collection
DefaultAssay(DC) <- "peaks"
pfm = readRDS("motifListPando.rds")
DC <- AddMotifs(
  object = DC,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  pfm = pfm
)

# Add chromvar z-scores
library(BiocParallel)
register(SerialParam())
register(MulticoreParam(40, progressbar = TRUE))
DC <- RunChromVAR(DC, genome = BSgenome.Hsapiens.UCSC.hg38)

# subset expressed TFs

# peaks-to-gene (Signac)
Links = readRDS("DC.rds")
lnk <- Links(object = Links[["ATAC"]])
lnk = readRDS("/broad/sankaranlab/cohn/data/GRN/velocyto/scVelo/peakToGene_links.rds")

# motifs data
motifsNames = data.frame(fread("TF_Information_all_motifs_plus.txt",header=T))
motifxTF <- unlist(DC@assays$peaks@motifs@motif.names)
motifxTF <- cbind(names(motifxTF), motifxTF)
motifxTF <- data.frame(motif = names(DC@assays$peaks@motifs@motif.names), gene = unlist(DC@assays$peaks@motifs@motif.names))
motifxTF$gene = toupper(motifxTF$gene)
TF_names <- motifxTF$gene 
TF_names2 <- unique(motifsNames$TF_Name)
TF_total = c(TF_names,TF_names2) %>% unique()
lnk.keep <- lnk[(abs(x = lnk$score) > 0) & lnk$gene %in% TF_total]
TF_names <- intersect(TF_total , DC@assays$RNA %>% rownames())
genePeaks <- lnk.keep %>% as.data.frame() %>% select(peak, gene) 
TF_total <- intersect(TF_total , lnk.keep$gene)

### extract RNA/peaks data
DefaultAssay(DC) <- "RNA" 
DC <- NormalizeData(DC)
RNAdata <-FetchData(DC,vars = c("latent_time","barcode","velocity_pseudotime","group",TF_total), slot = 'data')
RNA_melt <- data.table::melt(setDT(RNAdata), id.vars = c("latent_time","velocity_pseudotime","group","barcode"),
                             variable.name = "gene", value.name = "expression") %>% tibble::as_tibble()
# load TFs peaks links
peak_count <- t(GetAssayData(object = DC, slot = "data", assay = "ATAC"))
results_list <- lapply(unique(lnk.keep$gene), function(Genepeaks) {
  print(paste("Processing gene:", Genepeaks))
  genePeaksFiltered <- subset(genePeaks, gene == Genepeaks)
  print(head(genePeaksFiltered)) 
  peaksFiltered <- peak_count[, colnames(peak_count) %in% genePeaksFiltered$peak]
  peaksFiltered <- t(peaksFiltered)
  depth.vec <- Matrix.utils::aggregate.Matrix(x = peaksFiltered, fun = "sum", 
                                              groupings = rep(1, nrow(peaksFiltered)))[1, ] %>% data.frame() %>%
    rownames_to_column(var = "barcode") %>% as_tibble() %>%
    rename(!!Genepeaks := ".")  
  print(head(depth.vec))
  ATAC_data <- cbind(pseudotime, depth.vec) %>% as_tibble()
  return(ATAC_data)
})
combined_results <- Reduce(function(x, y) merge(x, y, by=c("barcode","latent_time")), results_list)
plot_data <- data.table::melt(setDT(combined_results), id.vars = c("barcode","latent_time"),
                              variable.name = "gene", value.name = "peaks")
ATAC_RNA = merge(RNA_melt,plot_data, by = c("barcode","latent_time","gene"))

# fitting to pseudotime using GAM model
process_gene <- function(subsetData) {
  fit_expression <- mgcv::gam(expression ~ s(latent_time, k = 6, bs = "cr", sp = 10),family = gaussian, data = subsetData)
  subsetData$expression_fit <- fit_expression$fitted.values
  fit_expression <- mgcv::gam(peaks ~ s(latent_time, k = 6, bs = "cr", sp = 10),family = gaussian, data = subsetData)
  subsetData$peaks_fit <- fit_expression$fitted.values
  return(subsetData)
}

data.plot <- ATAC_RNA %>% group_by(gene) %>% do(process_gene(.)) %>% as_tibble()

# correlation
corResults <- setDT(data.plot)[,
                               .(
                                 r_expression = cor(latent_time, expression_fit, method = "pearson"),
                                 r_peak = cor(latent_time, peaks_fit, method = "pearson"),
                                 p_expression = cor.test(latent_time, expression_fit)$p.value,
                                 p_peak = cor.test(latent_time, peaks_fit)$p.value
                               ),
                               by = gene
]

# plot trands

id = c("ZEB2","BATF3","RBPJ","GATA2") #cDC1
id = c("GATA2","IRF8","ZBTB16","TCF4") #pDC


plot = filter(data.plot, gene %in% id)
p = ggplot(plot,aes(group = gene, color = gene)) +
  geom_smooth(method="gam", se=F,lwd=2, mapping = aes(x = latent_time, y = expression_fit,fill = gene)) +
  #facet_wrap(~gene, scales = "free") +
  scale_color_manual(values=c("#86d1a6", "#e37868",'#7570b3',"#f0e24f"))+
  scale_fill_manual(values=c("#86d1a6", "#e37868",'#7570b3',"#f0e24f"))+
  xlab("Pseudotime") + ylab("Expression") +
  scale_x_continuous(expand = c(0,0))+
  theme(strip.background = element_rect(colour = "white", fill = "white")) + 
  theme(panel.border = element_blank(), axis.line = element_line()) + 
  theme(legend.position="top", legend.box = "horizontal")+
  theme(panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank()) + 
  theme(panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank()) + 
  theme(axis.ticks.length = unit(0.25, "cm")) +
  theme(panel.background = element_rect(fill = "white"))
ggsave(plot = p, "corr.peak.genes.pdf",height = 6, width = 8)

# plot umap
df <- data.table(Embeddings(DC, "wnn.umap"),pseudotime= DC@meta.data$latent_time)
color_palette <- colorRampPalette(c("white", '#e6aaaa', '#92abd4'))(n)

scvelo_pt <- ggplot(df)+
  geom_scattermore(aes(x=wnnUMAP_1, y=wnnUMAP_2, color=pseudotime), alpha=0.95,pointsize = 9,pixels    = c(3000, 3000))+
  labs(colour = "Pseudotime") +
  scale_color_gradientn(colours = color_palette,
                        guide = guide_colorbar(ticks.colour = "black",
                                               frame.colour = "black",
                                               barwidth=1.3)) +
  theme_bw() +
  theme(legend.title = element_text(color = "black", size = 6),
        legend.text = element_text(color = "black", size = 6),
        legend.key = element_rect(color = "black", size = 0.5),
        axis.line.x=element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=0.75),
        axis.line.y= element_blank(),
        strip.placement = "outside",
        axis.ticks.length = unit(0.25, "cm"),
        strip.background = element_rect(fill = "gray95",colour="white"),
        panel.background=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank())

ggplot2::ggsave(
  filename = "umap.pdf",
  plot = scvelo_pt,
  width = 5,
  height = 4,
  device = cairo_pdf,
  bg = "white"
)

