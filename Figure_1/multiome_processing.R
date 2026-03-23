cat("\f") 
rm(list = ls())
library(Signac)
library(Seurat)
library(JASPAR2024)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(chromVAR)
library(scCustomize)
set.seed(1234)
# load the RNA and ATAC data
counts <- Read10X_h5("filtered_feature_bc_matrix.h5")
fragpath <- "atac_fragments.tsv.gz"

# get gene annotations for hg38
annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotation) <- 'UCSC'
genome(annotation) <- "hg38"

# remove mt genes 
rna_counts <- counts$`Gene Expression`
rna_counts= rna_counts[-which(grepl("^MT-",rna_counts@Dimnames[[1]])),]  

DC <- CreateSeuratObject(
  counts = rna_counts,
  assay = "RNA",
  project = "DC_multiome", min.cells = 3
)
# Use peaks in standard chromosomes
atac_counts <- counts$Peaks
grange.counts <- StringToGRanges(rownames(atac_counts), sep = c(":", "-"))
grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
atac_counts <- atac_counts[as.vector(grange.use), ]

# Add ATAC data
DC[["ATAC"]] <- CreateChromatinAssay(
  counts = atac_counts,
  sep = c(":", "-"),
  fragments = fragpath,
  genome = 'hg38',
  annotation = annotation,
  min.cells = 3)

# QC
DefaultAssay(DC) <- "ATAC"
DC <- NucleosomeSignal(DC)
DC <- TSSEnrichment(DC, fast = F)
DC$high.tss <- ifelse(DC$TSS.enrichment > 2, 'High', 'Low')
TSSPlot(DC, group.by = 'high.tss') + NoLegend()
DC$nucleosome_group <- ifelse(DC$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
FragmentHistogram(object = DC, group.by = 'nucleosome_group')

# RNA analysis
DefaultAssay(DC) <- "RNA"
DC <- SCTransform(DC,vst.flavor = "v2")
DC <- RunPCA(DC)
print(DC[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(DC, dims = 1:5, reduction = "pca")
DimPlot(DC, reduction = "pca")
DC = RunUMAP(DC,dims = 1:50, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')
DimPlot(DC, reduction = "umap.rna")

# ATAC analysis
DefaultAssay(DC) <- "ATAC"
peaks <- CallPeaks(
  object = DC,
  macs2.path = "bin/macs2"
)

# remove peaks on nonstandard chromosomes and in genomic blacklist regions
peaks <- keepStandardChromosomes(peaks, pruning.mode = "coarse")
peaks <- subsetByOverlaps(x = peaks, ranges = blacklist_hg38_unified, invert = TRUE)

# create a new assay using the MACS2 peak set and add it to the Seurat object
macs2_counts <- FeatureMatrix(
  fragments = Fragments(DC),
  features = peaks,
  cells = colnames(DC)
)
DC[["peaks"]] <- CreateChromatinAssay(
  counts = macs2_counts,
  fragments = fragpath,
  annotation = annotation
)
DefaultAssay(DC) <- "peaks"
DC <- RegionStats(DC, genome = BSgenome.Hsapiens.UCSC.hg38)

## add motifs
pfm <- TFBSTools::getMatrixSet(
  x = JASPAR2022,
  opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE)
)

DC <- AddMotifs(
  object = DC,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  pfm = pfm
)

## add chromvar
DC <- RunChromVAR(DC, genome = BSgenome.Hsapiens.UCSC.hg38)

## run multimodal pipline
DC <- RunTFIDF(DC)
DC <- FindTopFeatures(DC, min.cutoff = 'q0')
DC <- RunSVD(DC)
DepthCor(DC)
DC <- RunUMAP(DC, reduction = 'lsi', dims = 2:50, reduction.name = "umap.atac", reduction.key = "atacUMAP_")
DC <- FindMultiModalNeighbors(
  DC, reduction.list = list("pca", "lsi"), dims.list = list(1:50, 2:50))
DC <- RunUMAP(DC, 
              nn.name = "weighted.nn", 
              reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
DC <- RunTSNE(DC, 
              nn.name = "weighted.nn", 
              reduction.name = "wnn.tsne", reduction.key = "wnnTSNE_")
DC <- FindClusters(DC, graph.name = "wsnn", algorithm = 3, verbose = FALSE)


DefaultAssay(DC) <- "SCT"
DimPlot(DC, reduction = "wnn.umap", group.by = "ident", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("WNN")


# Integrate DC datasets
DC1 = readRDS("DC_multiome.rds")
DC2 = readRDS("DC_multiome_2.rds")
DC3 = readRDS("DC_multiome_3.rds")

## Add metadata to anontate
DC1 <-AddMetaData(DC1, metadata="DC1", col.name="replicate")
DC1 <- RenameCells(DC1, add.cell.id = "DC1")
DC2 <-AddMetaData(DC2, metadata="DC2", col.name="replicate")
DC2 <- RenameCells(DC2, add.cell.id = "DC2")
DC3 <-AddMetaData(DC3, metadata="DC3", col.name="replicate")
DC3 <- RenameCells(DC3, add.cell.id = "DC3")

## Create a object containing all data
obj.list <- list(DC3 = DC3, DC1 = DC1, DC2 = DC2)

## Transform and normalize
obj.list <- lapply(X = obj.list, FUN = SCTransform)
rm(DC)
rm(DC2)
rm(DC3)
## Integrate data: RNA
features <- SelectIntegrationFeatures(object.list = obj.list
                                      obj.list <- lapply(X = obj.list, FUN = function(x) {
                                        x <- ScaleData(x, features = features, verbose = FALSE)
                                        x <- RunPCA(x, features = features, verbose = FALSE)
                                      })
                                      obj.list <- PrepSCTIntegration(object.list = obj.list, anchor.features = features)
)

# Finding integration anchors using the RNA modality
anchors <- FindIntegrationAnchors(object.list = obj.list, normalization.method = "SCT",
                                  anchor.features = features,reduction = "cca",k.anchor = 20)

all_genes <- Reduce(intersect, lapply(obj.list, rownames))
combined.sct.rna <- IntegrateData(anchorset = anchors, normalization.method = "SCT", new.assay.name ='SCTintegrated')
combined.sct.rna <- RunPCA(combined.sct.rna, verbose = FALSE)
combined.sct.rna <- RunUMAP(combined.sct.rna, reduction = "pca", dims = 1:50, verbose = FALSE,reduction.key = 'rnaUMAP_',reduction.name = 'umap.rna')

# Integrate data: ATAC
DefaultAssay(combined.sct) <- "peaks" 
combined.sct <- FindTopFeatures(combined.sct, min.cutoff = 'q0')
combined.sct <- RunTFIDF(combined.sct)
combined.sct <- RunSVD(combined.sct)
combined.sct <- RunUMAP(combined.sct, reduction = "lsi", dims = 2:50,reduction.name = "umap.atac", reduction.key = "atacUMAP_")
# integrate LSI embeddings
combined.sct <- IntegrateEmbeddings(
  anchorset = anchors,
  reductions = combined.sct[["lsi"]],
  new.reduction.name = "integrated_lsi",
  dims.to.integrate = 1:50
)
combined.sct <- RunUMAP(combined.sct, reduction = 'integrated_lsi', dims = 2:50, reduction.name = "umap.atac", reduction.key = "atacUMAP_")


# tranfer embedding from SCTintegrated
combined.sct@assays[["SCTintegrated"]] <- combined.sct.rna@assays[["SCTintegrated"]]
combined.sct@reductions<- append(combined.sct@reductions,combined.sct.rna@reductions)

# WNN
combined.sct <- FindMultiModalNeighbors(
  combined.sct, reduction.list = list("pca", "integrated_lsi"), dims.list = list(1:50, 2:50))
combined.sct <- RunUMAP(combined.sct, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")


combined.sct <- FindClusters(combined.sct, graph.name = "wsnn", algorithm = 3, verbose = FALSE)

# DEGs
markersRNA<- presto::wilcoxauc(DC, 'seurat_clusters', assay = 'data')
TopRNA<- presto::top_markers(markersRNA, n = 5, auc_min = 0.5, pct_in_min = 20, pct_out_max = 20, pval_max = 0.05)

genes  <- unique(c("GZMB", "PTCRA", "SPIB", "CLEC4C", "TCF4", "IL3RA", "IRF7", "TLR7", "PLAC8", #pDC
                   "LEF1","THEMIS","CD79A","NKG7","TRBC2","IGLL1","MS4A1", #NK/B/T
                   "AVP","CD34","SPINK2","PRSS57","CYTL1","EGFL7","SMIM24","MYB","LAPTM4B","CDK6", #HSPC
                   "IGLL1", "CD34", "TAL1","HEMGN","GFI1B","GATA1","SMIM1","NFE2", #progenitor cells
                   "EPOR","CNRIP1","MYL4","TAL1","LAT","HEMGN","NFE2","GFI1B","GMPR","GATA1","KLF1","TFR2","SMIM1", #MEP and erythroblasts
                   "HDC","SLC18A2","CPA3","IL18R1","GRAP2",  #MAST
                   "LAMP3","CD273","IDO1", "CD274", "PDCD1LG2","CCR7","CCL19","CCL12","FSCN1","TNFSF9","CCL22","MARCKSL1","CD80","CD83","CCL21","BIRC3",
                   "CX3CR1", "AXL", "SIGLEC6", #ASDC
                   "CADM1", "CLEC9A", "IDO1", "C1orf54", "BATF3", "SLAMF8", "SNX22", "CPNE3", "GCSAM", "THBD" ,"CLNK","CXCR4","CLEC12A", #cDC1
                   "CD1A", "CTSH", "CD1E","IL1R2","CLEC4A","CLEC10A","CX3CR1","IL13RA1","MRC1","IL22RA2", #cDC2
                   "FCER1A", "ENHO", "GSN","FCN1","S100A9","VCAN","CD14","CSF1R","IL1R1" #DC3/ mono
))
Stacked_VlnPlot(seurat_object = combined.sct, features = DC_cells, x_lab_rotate = TRUE, group.by = "CellTypeAll")

combined.sct <- RenameIdents(combined.sct,
                             "0" = "pDC_1",
                             "1" = "cDC2_1",
                             "2" = "cDC2_2",
                             '3' = "cDC1_1",
                             "4" = "pDC_2",
                             "5" = "cDC2_3",
                             "6" = "pDC_3",
                             "7" = "cDC2_4",
                             "8" = "cDC2_5",
                             "9" = "cDC2_6",
                             "10" = "MEP / Ery / MAST_1",
                             '11' = 'pDC_4',
                             "12" = "mDC",
                             "13" = "pDC_5",
                             "14" = "pDC_6",
                             "15" = "MEP / Ery / MAST_2",
                             "16" = "cDC2_7",
                             "17" = "cDC1_2",
                             "18" = "cDC2_8",
                             "19" = "cDC2_9",
                             "20" = "MEP / Ery / MAST_3",
                             "21" = "cDC2_10",
                             "22" = "cDC2_11",
                             "23" = "T /NK /B"
                             
)

combined.sct$CellTypeAll = Idents(combined.sct)
my_cols <- c('pDC_1'='#97bee8','pDC_2'='#97c9e8','pDC_3'='#9cafdb','pDC_4'='#9cbcdb',pDC_5='#9acfdb',pDC_6='#bfdaf5',
             'cDC2_1'='#e0a09d','cDC2_2'='#e69c85','cDC2_3'='#d9936a','cDC2_4'='#e89974',
             'cDC2_5'='#ea9b65','cDC2_6'='#d9934a', 'cDC2_7'='#e88e51','cDC2_8 / DC3'='#a3614b',
             'cDC2_9'='#e8a679','cDC2_10'='#e89974','cDC2_11'='#e6803c','mDC'='#e683cc','MEP / Ery / MAST_4'='#e5eb8d',
             'MEP / Ery / MAST_1'='#D4D915','MEP / Ery / MAST_2' = "#c6db8a",'MEP / Ery / MAST_3'='#c6d660',
             'cDC1_1'='#78c27c','cDC1_2'='#98e693','T /NK /B'='#d9d2ab')


pdf("umaps.integrated.pdf", height = 8, width = 14)
pdf("WNN.replicates.pdf", height = 8, width = 10) 
p1= DimPlot(combined.sct, reduction = "wnn.umap",cols = my_cols,label = TRUE)  #group.by
p1= DimPlot(combined.sct, reduction = "wnn.umap",group.by ="replicate",cols = c('DC1' = '#ebac98','DC2'= '#bacee3','DC3' = "#c9debd"))  #group.by
plot(p1  &  NoAxes())
dev.off()
combined.sct$CellType = Idents(combined.sct)

saveRDS("combined.sct.rds")