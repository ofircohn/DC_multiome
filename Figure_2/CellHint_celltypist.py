### Figure 2 -CellHint / CellTypist
import scanpy as sc
import cellhint
import seaborn as sns
import os
import celltypist

# import anndata
DC = sc.read('DC_AnnData.h5ad')
DC.obs_names_make_unique()
liver = sc.read('liver.h5ad')
liver.obs_names_make_unique()
BM = sc.read('BM_new.h5ad')
BM.obs_names_make_unique()
skin_healthy = sc.read('skin_Healthy.h5ad')
skin_atopic = sc.read('skin_Eczema.h5ad')
skin_atopic.obs_names_make_unique()
ileum_Crohn= sc.read('ileum_Crohn.h5ad')
ileum_Crohn.obs_names_make_unique()
ileum= sc.read('ileum_healthy.h5ad')
ileum.obs_names_make_unique()
tonsil= sc.read('tonsil.h5ad')
tonsil.obs_names_make_unique()
colon = sc.read('Colon.h5ad')
colon.obs['CellType'] = colon.obs['cell_type']
print(list(colon.obs['CellType'].unique()))
colon = colon[colon.obs['CellType'].str.contains('cDC2|cDC1|pDC|Lymphoid DC')]
colon.obs_names_make_unique()
gut = sc.read('myeloid_log_counts02_v2.h5ad')
gut.obs['CellType']= gut.obs['annotation']
print(gut.obs.columns)
gut = gut[gut.obs['CellType'].str.contains('cDC2|cDC1|pDC|Lymphoid DC')]
gut.obs_names_make_unique()
lymph = sc.read('Lymph.h5ad')
lymph = lymph[lymph.obs['CellType'].str.contains('cDC1|mDC|pDC')]
lymph.obs_names_make_unique()
spleen = sc.read('spleen.h5ad')
print(spleen.obs['CellType'].unique())
spleen.obs['CellType'] = spleen.obs['CellType'].replace('DC_plasmacytoid', 'pDC')
spleen.obs_names_make_unique()
oesophagus= sc.read('oesophagus.cellxgene.h5ad')
oesophagus.obs['CellType']= oesophagus.obs['Celltypes_updated_July_2020']
oesophagus = oesophagus[oesophagus.obs['CellType'].str.contains('Dendritic_Cells')]
oesophagus.obs['CellType'] = oesophagus.obs['CellType'].replace('Dendritic_Cells', 'cDCs')
oesophagus.obs_names_make_unique()
lungs = sc.read('lungs.h5ad')
lungs.obs_names_make_unique()
IPF = sc.read('IPF.h5ad')
IPF.obs_names_make_unique()
PBMC = sc.read('PBMC_healthy.h5ad')
PBMC.obs_names_make_unique()
PBMC_covid = sc.read('PBMC_covid.h5ad')
PBMC_covid.obs_names_make_unique()

# merge datasets
adatas = {"DC": DC, "PBMC_Healthy": PBMC, "skin": skin_healthy, "BM": BM, "liver": liver,
  "lymph": lymph, "Tonsil": tonsil, "PBMC_covid": PBMC_covid, 
  "IPF": IPF, "Colon": colon, "Atopic": skin_atopic, "lungs": lungs, "gut": gut, 
  "spleen":spleen,"ileum":ileum,"ileum_Crohn":ileum_Crohn,"oesophagus":oesophagus}
adatas = sc.concat(adatas, label="dataset_name",  join="outer", fill_value=0)

# normalize
columns_to_keep = ['CellType', 'dataset_name']
adatas.obs = adatas.obs[columns_to_keep]
adatas.layers["counts"] = adatas.X.copy()
adatas.obs["batch"].value_counts()
adatas.obs.dataset_name.value_counts()
sc.pp.normalize_total(adatas, target_sum = 1e4)
sc.pp.log1p(adatas)
sc.pp.highly_variable_genes(adatas, batch_key = 'dataset_name', subset = True)
sc.pp.scale(adatas, max_value = 10)
sc.tl.pca(adatas)
sc.pp.neighbors(adatas)
sc.tl.umap(adatas)
adata = sc.read('/broad/sankaranlab/cohn/data/TRS/datasets/CellHint/new/adata_before_cellhint.h5ad')

#plot
sc.settings.set_figure_params(dpi=300, fontsize=10, dpi_save=300, figsize=(5,5), format='pdf')
sc.pl.umap(adatas, color = ['dataset_name', 'CellType'], wspace = 0.5, save="umap.pdf")
sc.pl.umap(adatas,color="leiden",frameon=False,title="Cluster",save="umap_leiden.pdf")
sc.pl.umap(adatas, color ='CellType', wspace = 0.5, save="umap.pdf")

# Perform cell type harmonisation


alignment = cellhint.harmonize(adatas, 'dataset_name', 'CellType')
alignment.write('all_DC_alignment.pkl')
dist_mat = alignment.base_distance.to_meta()
dist_mat.iloc[:5, :5]
member_mat = alignment.base_distance.to_meta(turn_binary = True)
member_mat.iloc[:5, :5]
flag = member_mat.index.str.contains('DC')
sns_plot = sns.clustermap(member_mat.loc[flag,flag])
sns_plot.figure.savefig("heatmapDC.pdf")
member_mat.to_csv('member.csv', sep = ',', index = False)

# integrate
toMatch = ["ASDC", "mono", "aDC", "Mono", "Neu", "DC5", "DC4", "MEP", "NK", "FDC", "CD34", "IL7"]
pattern = '|'.join(toMatch)
adatas = adatas[~adatas.obs['CellType'].str.contains(pattern)]
adata.layers["counts"] = adata.X.copy()
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adatas)
adata.raw = adatas
sc.pp.highly_variable_genes(adatas, batch_key='dataset_name', subset=True)
sc.pp.scale(adatas, max_value=10)
sc.tl.pca(adatas)
sc.pp.neighbors(adatas)
sc.tl.umap(adatas)
cellhint.integrate(adatas, 'dataset_name')
sc.tl.umap(adatas)

#UMAP plot
sc.pl.umap(adatas, color='CellType', legend_loc='on data', legend_fontsize=4, save="umap_celltype.pdf")

# Export 
adatas.write_h5ad("adata_cellhint_celltype.h5ad")

# Cell celltypist model
from matplotlib.colors import  LinearSegmentedColormap
import matplotlib

adatas = celltypist.annotate(adatas, model = 'Immune_All_Low.pkl', majority_voting = True).to_adata()
predictions = celltypist.annotate(adata, model = '/broad/sankaranlab/cohn/data/TRS/datasets/CellHint/new/DC_Human.pkl', majority_voting = True)
celltypist.dotplot(predictions, use_as_reference = 'CellType', use_as_prediction = 'majority_voting',
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
dp = celltypist.dotplot(predictions, use_as_reference = 'CellType', use_as_prediction = 'majority_voting',
cmap=LinearSegmentedColormap.from_list('anglemap', ["#d3d3d3", "#a5cde3", "#6dafd6", '#08306b'], N=256, gamma=1.5),
return_fig=True)
axes_dict = dp.get_axes()
dp.make_figure()
axes_dict["mainplot_ax"].set_axisbelow(True)
dp.savefig("CellTypist_model.pdf",dpi=300)


# Format and clean the anotated data
# mregDC = mDC
adatas.obs['CellType_all'] = adatas.obs['CellType_all'].replace({'spleen_pDC': 'pDC', 
'BM_pDCs': 'pDC', 'Colon_pDC': 'pDC', 'DC_pDC': 'pDC','gut_pDC': 'pDC','IPF_pDCs': 'pDC',
'liver_pDCs': 'pDC','lungs_pDCs': 'pDC', 'lymph_pDC': 'pDC', 'PBMC_covid_pDC': 'pDC', 'PBMC_Healthy_pDC': 'pDC',
'spleen_pDC': 'pDC', 'Tonsil_PDC': 'pDC'})
adatas.obs['CellType_all'] = adatas.obs['CellType_all'].replace({
  'spleen_DC_activated': 'mregDC', 'skin_migDC': 'mregDC', 'liver_Mig.cDCs': 'mregDC', 
  'lymph_mDC': 'mregDC','ileum_Mature DCs': 'mregDC','ileum_Crohn_Mature DCs': 'mregDC',
  'gut_Lymphoid DC': 'mregDC', 'Atopic_migDC': 'mregDC','Colon_Lymphoid DC': 'mregDC',
  'DC_mregDC': 'mregDC'})
  
adatas.obs['CellType_all'] = adatas.obs['CellType_all'].replace({
  'spleen_DC_1': 'cDC1', 'skin_DC3': 'cDC1', 'lymph_cDC1': 'cDC1', 'Atopic_DC3': 'cDC1',
  'liver_cDC1s': 'cDC1','ileum_DC1': 'cDC1','ileum_Crohn_DC1': 'cDC1',
  'gut_cDC1': 'cDC1', 'Tonsil_DC1 precursor': 'cDC1','PBMC_covid_DC1': 'cDC1',
  'PBMC_Healthy_DC1': 'cDC1','DC_cDC1': 'cDC1', 'Colon_cDC1': 'cDC1'})  

adatas.obs['CellType_all'] = adatas.obs['CellType_all'].replace({
  'DC_cDC2': 'cDC2', 'PBMC_Healthy_DC2': 'cDC2', 'liver_cDC2s': 'cDC2', 
  'PBMC_covid_DC2': 'cDC2','ileum_DC2 CD1D-': 'cDC2','Atopic_DC1': 'cDC2',
  'Tonsil_DC2': 'cDC2', 'ileum_Crohn_DC2 CD1D-': 'cDC2','skin_DC1': 'cDC2',
  'gut_cDC2': 'cDC2', 'Colon_cDC2': 'cDC2','spleen_DC_2': 'cDC2'})  
adatas.obs['CellType_all'] = adatas.obs['CellType_all'].replace({
  'oesophagus_cDCs': 'cDCs', 'lungs_cDCs': 'cDCs', 'IPF_cDCs': 'cDCs', 
  'BM_cDCs': 'cDCs'})  
  
toMatch = ['PBMC_Healthy_DC3','PBMC_Healthy_ASDC',"PBMC_covid_DC3",'PBMC_Healthy_DC3','Tonsil_DC3',
"PBMC_covid_ASDC", "Tonsil_DC1 mature", "skin_DC2", "Atopic_DC2",'ileum_Crohn_DC2 CD1D','ileum_DC2 CD1D']
adatas.obs['CellType_all'] = adatas.obs['CellType_all'].astype(str)
pattern = '|'.join(toMatch)

# Filter the adatas to exclude the matching rows
adatas = adatas[~adatas.obs['CellType_all'].str.contains(pattern)]
unique_values = adatas.obs['CellType_all'].unique()
sc.pl.dotplot(adatas, marker_protein, groupby="CellType_all",save= "heatmapDC_typist_new.pdf")
sc.tl.umap(adatas)
sc.pl.umap(adatas, color='CellType_all', legend_loc='on data', legend_fontsize=4, save="umap_clean.pdf")
df = pd.DataFrame(adatas.obsm["X_umap"], index=adatas.obs_names, columns=["UMAP_1", "UMAP_2"])
df = pd.DataFrame({
    "original_cellType": adatas.obs["CellType"].values,
    "dataset": adatas.obs["dataset_name"].values,
    "cell_type": adatas.obs["CellType_all"].values,
    "UMAP_1": adatas.obsm["X_umap"][:, 0],
    "UMAP_2": adatas.obsm["X_umap"][:, 1]
}, index=adatas.obs_names)
df.to_csv("umap.csv"), sep="\t",index=False)
adatas.write_h5ad("adata_cellhint_celltypist_cleaned.h5ad")

# check number of cells
cell_counts = adata.obs.groupby(['dataset_name', 'CellType_all']).size().reset_index(name='cell_count')
cell_counts.to_csv("integrated_counts.csv", index=False)
counts = adata.obs[['CellType_all']].value_counts()

# find DEGs
sc.tl.filter_rank_genes_groups(adatas,max_out_group_fraction=0.1,min_in_group_fraction=0.6, min_fold_change = 1 ,key_added="rank_genes_groups_filtered")

# Iterate over each cluster
for CellType_all in groups:
    gene_names = adatas.uns['rank_genes_groups_filtered']['names'][CellType_all]
    valid_genes = [gene for gene in gene_names if str(gene) != 'nan']
    cluster_genes = valid_genes[:topX] if len(valid_genes) >= topX else valid_genes
    cluster_genes_dict[CellType_all] = cluster_genes

for CellType_all in groups:
  diffexp_df = sc.get.rank_genes_groups_df(adatas, key="rank_genes_groups_filtered", group=CellType_all)
  diffexp_df = diffexp_df[~diffexp_df["names"].isnull()]
  diffexp_df['ref_attr_value'] = CellType_all
  result_dataframes.append(diffexp_df)
  diffexp_df = pd.concat(result_dataframes, axis=0).reset_index(drop=True)
combined_df = pd.concat(result_dataframes, axis=0).reset_index(drop=True)
combined_df.to_csv("/broad/sankaranlab/cohn/data/TRS/datasets/CellHint/new/differential_expression_results.csv", index=False)


### plot Cellhint table (in R)
library(grid)
library(ggplot2)

table = read.csv("member.csv")
rownames(table) = colnames(table)
toMatch <- c(".ASDC", "mono", "aDC","Mono","Neu","DC5","DC4","MEP","NK","FDC","CD34","IL7")
DC.mat = table[!grepl(paste(toMatch,collapse="|"),colnames(table)),!grepl(paste(toMatch,collapse="|"),colnames(table))]

q95_pos <- quantile(DC.mat[DC.mat > 0], probs = 0.95)
q95_neg <- quantile(DC.mat[DC.mat < 0], probs = 0.05)
col_fun = circlize::colorRamp2(c(0, (max(DC.mat))), c("#6e0e4c", "#f5e498"))

order1 <- seriation::seriate(dist(DC.mat),method="OLO_ward")
order2 <- seriation::seriate(dist(t(DC.mat)),method="OLO_ward")
col.dend = as.dendrogram(order2[[1]])
dend_k = dendextend::find_k(col.dend)
color = scales::brewer_pal(palette = "Set1")(5)
col.dend = dendextend::color_branches(col.dend, dend_k$k,col = color)
ComplexHeatmap::Heatmap(
  DC.mat,
  name="Cell-type\nSimilarity",
  col =col_fun,
  cluster_columns = col.dend,
  cluster_rows = as.dendrogram(order1[[1]]),
  row_dend_side = "left",
  show_row_dend = FALSE,
  row_names_side = "right",
  border = "black",
  row_title_gp = gpar(fontsize=8),
  row_title_rot = 0,
  row_names_gp = gpar(fontsize=8),
  column_names_gp = gpar(fontsize=8),
  heatmap_legend_param = list(title_gp = gpar(fontsize = 8, fontface = "bold"),
                              labels_gp = gpar(fontsize = 8),border = "black")
)


