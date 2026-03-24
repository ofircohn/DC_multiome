library(velocyto.R)
library(pagoda2)
library(SeuratWrappers)
library(Seurat)
library(reticulate)

CheckPackage <- function(package, repository, ...) {
  if (!requireNamespace(package = basename(path = package), quietly = TRUE)) {
    if (interactive()) {
      message("Package ", package, " is not yet installed")
      message("Install now?")
      choice <- menu(choices = c('yes', 'no'))
      if (choice == 1) {
        repository <- match.arg(
          arg = tolower(x = repository),
          choices = c('github', 'bioconductor', 'cran')
        )
        switch(
          EXPR = repository,
          'github' = remotes::install_github(repo = package, ...),
          'bioconductor' = BiocManager::install(pkgs = package, ...),
          'cran' = install.packages(pkgs = package, ...),
          stop("Unknown repository ", repository, call. = FALSE)
        )
        return(invisible(x = NULL))
      }
    }
    stop("Unable to find package ", package, ", please install", call. = FALSE)
  }
}
ReadVelocity <- function(file, engine = 'hdf5r', verbose = TRUE) {
  CheckPackage(package = 'velocyto-team/velocyto.R', repository = 'github')
  if (verbose) {
    sink(file = stderr(), type = 'output')
    on.exit(expr = sink())
    ldat <- velocyto.R::read.loom.matrices(file = file, engine = engine)
  } else {
    invisible(x = capture.output(ldat <- velocyto.R::read.loom.matrices(
      file = file,
      engine = engine
    )))
  }
  return(ldat)
}


DC1 <- ReadVelocity(file = "gex_possorted_bam_MZ655.loom")
colnames(DC1$spliced) <- gsub("gex_possorted_bam_MZ655:", "", colnames(DC1$spliced))
colnames(DC1$spliced) <- gsub("x", "", colnames(DC1$spliced))
colnames(DC1$unspliced) <- gsub("gex_possorted_bam_MZ655:", "", colnames(DC1$unspliced))
colnames(DC1$unspliced) <- gsub("x", "", colnames(DC1$unspliced))
colnames(DC1$ambiguous) <- gsub("gex_possorted_bam_MZ655:", "", colnames(DC1$ambiguous))
colnames(DC1$ambiguous) <- gsub("x", "", colnames(DC1$ambiguous))
options(Seurat.object.assay.version = 'v3')
DC1 <- as.Seurat(DC1)

DC2 <- ReadVelocity(file = "gex_possorted_bam_NKHZU.loom")
colnames(DC2$spliced) <- gsub("gex_possorted_bam_NKHZU:", "", colnames(DC2$spliced))
colnames(DC2$spliced) <- gsub("x", "", colnames(DC2$spliced))
colnames(DC2$unspliced) <- gsub("gex_possorted_bam_NKHZU:", "", colnames(DC2$unspliced))
colnames(DC2$unspliced) <- gsub("x", "", colnames(DC2$unspliced))
colnames(DC2$ambiguous) <- gsub("gex_possorted_bam_NKHZU:", "", colnames(DC2$ambiguous))
colnames(DC2$ambiguous) <- gsub("x", "", colnames(DC2$ambiguous))
options(Seurat.object.assay.version = 'v3')
DC2 <- as.Seurat(DC2)

DC3 <- ReadVelocity(file = "gex_possorted_bam_Y0V1L.loom")
colnames(DC3$spliced) <- gsub("gex_possorted_bam_Y0V1L:", "", colnames(DC3$spliced))
colnames(DC3$spliced) <- gsub("x", "", colnames(DC3$spliced))
colnames(DC3$unspliced) <- gsub("gex_possorted_bam_Y0V1L:", "", colnames(DC3$unspliced))
colnames(DC3$unspliced) <- gsub("x", "", colnames(DC3$unspliced))
colnames(DC3$ambiguous) <- gsub("gex_possorted_bam_Y0V1L:", "", colnames(DC3$ambiguous))
colnames(DC3$ambiguous) <- gsub("x", "", colnames(DC3$ambiguous))
options(Seurat.object.assay.version = 'v3')
DC3 <- as.Seurat(DC3)

# Merge Seurat objects
DC <- merge(
  x = DC1,
  y = list(DC2,DC3),
  add.cell.ids = c("DC1","DC2","DC3"),
  merge.data = TRUE
)
colnames(DC) <- paste0(colnames(DC), "-1")

# Import the original seurat object
DC.merged = readRDS("integrated.multiomics.chromvar.clustered.rds")
DC <- subset(DC, cells = colnames(DC.merged))

spliced <- CreateAssayObject(GetAssayData(DC, assay = "spliced"))
unspliced <- CreateAssayObject(GetAssayData(DC, assay = "unspliced"))
ambiguous <- CreateAssayObject(GetAssayData(DC, assay = "ambiguous"))
DC.merged[["spliced"]] <- spliced
DC.merged[["unspliced"]] <- unspliced
DC.merged[["ambiguous"]] <- ambiguous

# export to scVelo

use_condaenv("mambaenv2", required = TRUE)
scv <- import("scvelo")
ad <- import("anndata", convert = FALSE)
scv$logging$print_version()

# extract parameters
RNA <- as.matrix(GetAssayData(DC.merged, assay = "RNA", slot = "counts"))
pc <- Embeddings(DC.merged, reduction = "pca")
umap <- Embeddings(DC.merged, reduction = "wnn.umap")
spliced <- as.matrix(GetAssayData(DC.merged, assay = "spliced", slot = "counts"))
unspliced <- as.matrix(GetAssayData(DC.merged, assay = "unspliced", slot = "counts"))
genes <- rownames(RNA)
spliced <- spliced[genes,]
unspliced <- unspliced[genes,]
dim(spliced)
dim(unspliced)

dfvar <- DC.merged@assays$RNA@meta.features
dfobs = DC.merged@meta.data

adata_DC <- ad$AnnData(
  X=t(RNA),
  obs=dfobs,
  var=dfvar,
  layers=list('spliced'=t(spliced), 'unspliced'=t(unspliced)),
  obsm=list('X_umap'=umap)) 
adata_DC
adata_DC$write('scVelo_annData_DC.h5ad', compression='gzip')

