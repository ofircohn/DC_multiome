# onek1k analysis

import scanpy as sc
import numpy as np
import pandas as pd
import anndata as ad

# Recluster

outdir = "/datasets/1k1"
OneK1K = ad.read_h5ad(f"{outdir}/OneK1K_updated_14_celltypes_980_donors.h5ad", backed="r")
sc.settings.figdir = outdir
dc = new[(new.obs["cell_label"] == "DC").values].to_memory()
dc.var_names = dc.var["GeneSymbol"].astype(str).values
dc.var_names_make_unique()
dc.layers["counts"] = dc.X.copy()
sc.pp.normalize_total(dc, target_sum=1e4)
sc.pp.log1p(dc)
sc.pp.highly_variable_genes(dc, n_top_genes=2000)
sc.pp.pca(dc, n_comps=30)
sc.pp.neighbors(dc)
sc.tl.leiden(dc, resolution=0.5)
sc.tl.umap(dc)

markers = ["LILRA4", "CLEC4C", "PTCRA", "IL3RA", "TCF4", "SPIB", "GZMB", "IRF7",   # pDC
           "CD1C", "FCER1A", "CLEC10A", "CLEC9A", "XCR1", "CD14", "LYZ"]           # cDC / mono
sc.pl.dotplot(dc, markers, groupby="leiden", standard_scale="var", save="_DC_recluster.pdf")
sc.pl.umap(dc, color=["leiden", "LILRA4", "CD1C"], save="_DC_recluster.pdf")   # remove "cxg_label"
dc.var_names = dc.var_names.str.strip()
dc.var_names_make_unique()
sc.settings.figdir = outdir
sc.pl.dotplot(dc, markers, groupby="leiden", standard_scale="var", save="_DC_recluster.pdf")
sc.pl.umap(dc, color=["leiden", "LILRA4", "CD1C"], save="_DC_recluster.pdf")
dc.write_h5ad(f"{outdir}/DC_reclustered.h5ad")
pdc_new = dc[dc.obs["leiden"].isin(["2", "5"])].copy()
print(pdc_new.n_obs, "pDCs")
ncells = pdc_new.obs["individual"].astype(str).value_counts()
print(ncells.size, "donors")
for k in [1, 3, 5, 10]:
    print(f">= {k} cells: {(ncells >= k).sum()} donors")

pdc_new.write_h5ad(f"{outdir}/pDC_from_updated.h5ad")

### extract specific position
pos=105391005-105391005   # hg19 position
grep -P "^14\t14:${pos}\t" genotypes/plink_merged_980_donors.bim   # confirm it exists, alleles T/C

#plink --bfile genotypes/plink_merged_980_donors --snp 14:${pos} \
#      --recode A --out rs2841277
#head -3 rs2841277.raw

import statsmodels.formula.api as smf
import anndata as ad
import pandas as pd
import scanpy as sc
import numpy as np
import scipy.sparse as sp

raw = pd.read_csv(f"{outdir}/rs2841277.raw", sep=r"\s+")
geno = raw.set_index("IID").iloc[:, -1].rename("T_count")
cov = pd.read_csv(f"{outdir}/covariates/updated_OneK1K_980_donors_PFs/DC_mx_pf50.txt",
                  sep=r"\s+", index_col="sampleid")
cov = cov[["sex", "age"] + [f"pc{i}" for i in range(1, 7)]]

# pDC pseudobulk: log2 CPM per donor
pdc_new = ad.read_h5ad(f"{outdir}/pDC_from_updated.h5ad")
donors = pdc_new.obs["individual"].astype(str).values
D = pd.get_dummies(donors)
pb = sp.csr_matrix(D.values.T.astype(np.float32)) @ pdc_new.layers["counts"]
pb = pd.DataFrame(pb.toarray(), index=D.columns, columns=pdc_new.var_names)
ncells = pd.Series(donors).value_counts()
cpm = np.log2(pb.div(pb.sum(1), axis=0) * 1e6 + 1)

genes = [g for g in ["PLD4"] if g in cpm.columns]
df = cpm[genes].join(geno, how="inner").join(cov, how="inner")
df["ncells"] = ncells.reindex(df.index).values
for g in genes:
    for k in [1, 3]:
        d = df[df.ncells >= k]
        m = smf.wls(f"{g} ~ T_count + sex + age + pc1+pc2+pc3+pc4+pc5+pc6",
                    d, weights=d.ncells).fit()
        print(f"{g:10s} >={k} cells  n={len(d):3d}  beta(T)={m.params['T_count']:+.3f}  p={m.pvalues['T_count']:.2e}")
        
        
##
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import rankdata, norm
def int_transform(x):
    r = rankdata(x); return norm.ppf((r - 0.5) / len(r))
def test_all(OneK1K, label, gene="PLD4", donor_col="individual", min_cells=3):
    dn = adata_sub.obs[donor_col].astype(str).values
    Dm = pd.get_dummies(dn)
    X = adata_sub.layers["counts"] if "counts" in adata_sub.layers else adata_sub.X
    j = list(adata_sub.var_names).index(gene)
    col = X[:, j].toarray().ravel(); lib = np.asarray(X.sum(1)).ravel()
    p = pd.DataFrame((sp.csr_matrix(Dm.values.T.astype(np.float32)) @ X).toarray(),
                     index=Dm.columns, columns=adata_sub.var_names)
    t = np.log2(p[[gene]].div(p.sum(1), axis=0) * 1e6 + 1).join(geno, how="inner").join(cov, how="inner")
    t["ncells"] = pd.Series(dn).value_counts().reindex(t.index).values
    t = t[t.ncells >= min_cells].copy()
    t["expr_int"] = int_transform(t[gene])
    covs = "sex + age + pc1+pc2+pc3+pc4+pc5+pc6"
    m1 = smf.wls(f"{gene} ~ T_count + {covs}", t, weights=t.ncells).fit()
    m2 = smf.ols(f"expr_int ~ T_count + {covs}", t).fit()
    return dict(cell_type=label, n_donors=len(t), n_cells=adata_sub.n_obs,
                mean_cpm=1e6 * col.sum() / lib.sum(), frac_expr=np.mean(col > 0),
                beta_log2=m1.params["T_count"], p_log2=m1.pvalues["T_count"],
                beta_int=m2.params["T_count"], p_int=m2.pvalues["T_count"])

rows = []
for ct in new.obs["cell_label"].cat.categories:
    a = new[(new.obs["cell_label"] == ct).values].to_memory()
    a.var_names = a.var["GeneSymbol"].astype(str).str.strip().values
    a.var_names_make_unique()
    rows.append(test_all(a, ct)); print(rows[-1]["cell_type"], "done")
    del a
rows.append(test_all(pdc_new, "pDC (reclustered)"))
rows.append(test_all(cdc, "cDC2 (reclustered)"))
res = pd.DataFrame(rows).sort_values("p_int")
pd.set_option("display.width", 200)
print(res.round(4).to_string(index=False))
res.to_csv(f"{outdir}/PLD4_rs2841277_all_celltypes.csv", index=False)


######## export PLD4 data from pDCs

out = df[["PLD4", "T_count", "ncells", "sex", "age"]].copy()
out.insert(0, "donor_id", out.index)
out["genotype"] = out["T_count"].map({0: "CC", 1: "CT", 2: "TT"})
out.to_csv(f"{outdir}/PLD4_rs2841277_pDC_boxplot_data.csv", index=False)

stats = []
for k in [1, 3]:
    d = df[df.ncells >= k]
    m = smf.wls("PLD4 ~ T_count + sex + age + pc1+pc2+pc3+pc4+pc5+pc6", d, weights=d.ncells).fit()
    stats.append(dict(min_cells=k, n_donors=len(d), beta_per_T=m.params["T_count"],
                      se=m.bse["T_count"], p=m.pvalues["T_count"]))
pd.DataFrame(stats).to_csv(f"{outdir}/PLD4_rs2841277_pDC_stats.csv", index=False)
