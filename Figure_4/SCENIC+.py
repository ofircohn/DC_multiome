# Figure 4 - pycisTopic
import warnings
warnings.simplefilter(action='ignore')
import pandas as pd
import pycisTopic
import matplotlib.pyplot as plt
pycisTopic.__version__
import os 
import scanpy as sc
import pyranges as pr
import pickle
import os
from pycisTopic.pseudobulk_peak_calling import *
from pycisTopic.pseudobulk_peak_calling import peak_calling
from pycisTopic.pseudobulk_peak_calling import export_pseudobulk
from pycisTopic.iterative_peak_calling import *
from pycisTopic.qc import *
import pybiomart as pbm

adata = sc.read_h5ad(os.path.join('DC.RNA.adata.h5ad'))

# fragments files 
fragments_dict = {'DC1': 'atac_fragments.tsv.gz',
                 'DC2': 'atac_fragments.tsv.gz',
                 'DC3': 'atac_fragments.tsv.gz'}

# create pseudobulk by major DC groups
cell_data = adata.obs
cell_data['group'] = cell_data['group'].astype(str) # set data type of the celltype column to str, otherwise the export_pseudobulk function will complain.
cell_data['group'].value_counts()
cell_data['barcode'] = [x.split('_')[1] for x in cell_data.index.tolist()]
target_url='http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes'
chromsizes=pd.read_csv(target_url, sep='\t', header=None)
chromsizes.columns=['Chromosome', 'End']
chromsizes['Start']=[0]*chromsizes.shape[0]
chromsizes=chromsizes.loc[:,['Chromosome', 'Start', 'End']]
chromsizes['Chromosome'] = [chromsizes['Chromosome'][x].replace('v', '.') for x in range(len(chromsizes['Chromosome']))]
chromsizes['Chromosome'] = [chromsizes['Chromosome'][x].split('_')[1] if len(chromsizes['Chromosome'][x].split('_')) > 1 else chromsizes['Chromosome'][x] for x in range(len(chromsizes['Chromosome']))]
chromsizes=pr.PyRanges(chromsizes)

bw_paths, bed_paths = export_pseudobulk(input_data = cell_data,
                 variable = 'group', # cell type annotation
                 sample_id_col = 'replicate',
                 chromsizes = chromsizes,
                 bed_path = outDir + 'consensus_peak_calling/pseudobulk_bed_files/',
                 bigwig_path = outDir + 'consensus_peak_calling/pseudobulk_bw_files/',
                 path_to_fragments = fragments_dict,
                 n_cpu = 1,
                 temp_dir = "/SCENIC/outs/",
                 normalize_bigwig = True,split_pattern = '__')

# Save
with open(outDir + 'consensus_peak_calling/pseudobulk_bed_files/bed_paths.pkl', 'wb') as f:
  pickle.dump(bed_paths, f)
with open(outDir + 'consensus_peak_calling/pseudobulk_bed_files/bw_paths.pkl', 'wb') as f:
  pickle.dump(bw_paths, f)

# peak calling
macs_path="/bin/macs2"
macs_outdir = outDir + 'consensus_peak_calling/MACS/'
if not os.path.exists(macs_outdir):
  os.mkdir(macs_outdir)
narrow_peaks_dict = peak_calling(macs_path,
                                 bed_paths,
                                 macs_outdir,
                                 genome_size=1.4e9, ## modified this to zebrafish genome size
                                 n_cpu=1,
                                 input_format='BEDPE',
                                 shift=73,
                                 ext_size=146,
                                 keep_dup = 'all',
                                 q_value = 0.05,
                                 _temp_dir = tmpDir + 'ray_spill')
                                 
with open(outDir + 'consensus_peak_calling/MACS/narrow_peaks_dict.pkl', 'wb') as f:
  pickle.dump(narrow_peaks_dict, f)
  
## Peak annotation
peak_half_width=250
path_to_blacklist="hg38-blacklist.v2.bed"
consensus_peaks=get_consensus_peaks(narrow_peaks_dict, peak_half_width, chromsizes=chromsizes, path_to_blacklist=path_to_blacklist)
consensus_peaks
# Write to bed
consensus_peaks.to_bed(path= outDir + 'consensus_peak_calling/consensus_regions.bed', keep=True, compression='infer', chain=False)

# QC
outDir = 'outs/'
pycistopic tss gene_annotation_list | grep Human
pycistopic tss get_tss \
    --output tss.bed \
    --name "hsapiens_gene_ensembl" \
    --to-chrom-source ucsc \
    --ucsc hg38
regions_bed_filename = os.path.join(outDir, "consensus_peak_calling/consensus_regions.bed")
tss_bed_filename = os.path.join(outDir, "qc", "tss.bed")

# Create text file with all pycistopic qc command lines.
pycistopic qc \
    --fragments atac_fragments.tsv.gz \
    --regions consensus_regions.bed \
    --tss tss.bed \
    --output outs/qc/DC1

pycistopic qc \
    --fragments atac_fragments.tsv.gz \
    --regions consensus_regions.bed \
    --tss tss.bed \
    --output outs/qc/DC2
    
pycistopic qc \
    --fragments atac_fragments.tsv.gz \
    --regions consensus_regions.bed \
    --tss tss.bed \
    --output outs/qc/DC3
                       
                       
                       
# creating a cisTopic object
path_to_regions = os.path.join("consensus_regions.bed")
path_to_blacklist="hg38-blacklist.v2.bed"
pycistopic_qc_output_dir = "outs/qc"
from pycisTopic.cistopic_class import create_cistopic_object_from_fragments
from pycisTopic.cistopic_class import *
import polars as pl

cistopic_obj_list = []
for sample_id in fragments_dict:
    sample_metrics = pl.read_parquet(
        os.path.join(pycistopic_qc_output_dir, f'{sample_id}.fragments_stats_per_cb.parquet')
    ).to_pandas().set_index("CB").loc[ sample_id_to_barcodes_passing_filters[sample_id] ]
    cistopic_obj = create_cistopic_object_from_fragments(
        path_to_fragments = fragments_dict[sample_id],
        path_to_regions = path_to_regions,
        path_to_blacklist = path_to_blacklist,
        metrics = sample_metrics,
        valid_bc = sample_id_to_barcodes_passing_filters[sample_id],
        n_cpu = 1,
        project = sample_id,
        split_pattern = '-'
    )
    cistopic_obj_list.append(cistopic_obj)
cistopic_obj = merge(cistopic_obj_list)
import pickle
pickle.dump(
    cistopic_obj,
    open(os.path.join("outs/cistopic_obj.pkl"), "wb")
)
print(cistopic_obj)          
              
              
# add metadata
adata = sc.read_h5ad(os.path.join('DC.RNA.adata.h5ad'))
cell_data = adata.obs
cell_data['group'] = cell_data['group'].astype(str) # set data type of the celltype column to str, otherwise the export_pseudobulk function will complain.
cell_data['barcode'] = [x.split('_')[1] for x in cell_data.index.tolist()]
cell_data.index = cell_data['barcode'] + '-' + cell_data['replicate'].astype("object")
cell_data.index = cell_data['barcode'] + '___' + cell_data['replicate'].astype("object")
cell_data.head()
cistopic_obj.add_cell_data(cell_data, split_pattern='-')
cistopic_obj.add_cell_data(cell_data)
                   
# load models
from pycisTopic.lda_models import *
import pickle
from pycisTopic.cistopic_class import *
from pycisTopic.lda_models import run_cgs_models
from pycisTopic.lda_models import run_cgs_models_mallet
outDir = 'outs/'

#wget https://github.com/mimno/Mallet/releases/download/v202108/Mallet-202108-bin.tar.gz
infile = open(outDir + 'cistopic_obj.pkl', 'rb')
cistopic_obj = pickle.load(infile)
print(cistopic_obj.cell_data.head())

mallet_path="bin/mallet"
os.environ['MALLET_MEMORY'] = '150G'
models=run_cgs_models_mallet(
    cistopic_obj,
    n_topics=[2, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50],
    n_cpu=10,
    n_iter=500,
    random_state=555,
    alpha=50,
    alpha_by_topic=True,
    eta=0.1,
    eta_by_topic=False,
    tmp_path = "/outs",
    save_path="/outs",
    mallet_path=mallet_path)
    
models = []
for n_topics in [2, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50]:
    # infile = open(outDir+ 'models/mallet_nt' + str(n_topics) + '.pkl', 'rb')
    infile = open('/outs/' + 'Topic' + str(n_topics) + '.pkl', 'rb')
    models.append(pickle.load(infile)) 
    infile.close()
pickle.dump(
    models,
    open(os.path.join('/outs/models.pkl'), "wb")
)
from pycisTopic.lda_models import *
outDir = 'outs/'

if not os.path.exists(os.path.join(outDir,"models")):
    os.makedirs(os.path.join(outDir,"models"))
model=evaluate_models(models,
                     select_model=None,
                     return_model=True,
                     metrics=['Arun_2010','Cao_Juan_2009', 'Minmo_2011', 'loglikelihood'],
                     plot_metrics=False,
                     save= outDir + 'models/model_selection.pdf')
model=evaluate_models(models,
                     select_model=50,
                     return_model=True,
                     metrics=['Arun_2010','Cao_Juan_2009', 'Minmo_2011', 'loglikelihood'],
                     plot_metrics=False,
                     save= outDir + 'models/model_selection_50.pdf')
cistopic_obj.add_LDA_model(model)
with open(outDir + 'cisTopicObject.pkl', 'wb') as f:
  pickle.dump(cistopic_obj, f)


# topic binarize
from pycisTopic.topic_binarization import binarize_topics
if not os.path.exists(os.path.join(outDir, 'topic_binarization')):
    os.makedirs(os.path.join(outDir, 'topic_binarization'))
from pycisTopic.topic_binarization import *
region_bin_topics_otsu = binarize_topics(cistopic_obj, method='otsu', plot=True, num_columns=5, save= outDir + 'topic_binarization/otsu.pdf')
region_bin_topics_top3k = binarize_topics(cistopic_obj, method='ntop', ntop = 3000)

# Save
with open(outDir + 'cisTopicObject.pkl', 'wb') as f:
  pickle.dump(cistopic_obj, f)

# DARs per cell type.
import pickle
import os
from pycisTopic.diff_features import *
import numpy as np

imputed_acc_obj = impute_accessibility(
    cistopic_obj,
    selected_cells=None,
    selected_regions=None,
    scale_factor=10**6
)
normalized_imputed_acc_obj = normalize_scores(imputed_acc_obj, scale_factor=10**4)
variable_regions = find_highly_variable_features(normalized_imputed_acc_obj, plot = False)
len(variable_regions)
markers_dict= find_diff_features(
    cistopic_obj,
    imputed_acc_obj,
    variable='group',
    var_features=None,
    contrasts=None,
    adjpval_thr=0.05,
    log2fc_thr=np.log2(1.5),
    n_cpu=10,
    split_pattern = '___'
)
print(markers_dict)

# Save
if not os.path.exists(os.path.join(outDir, 'candidate_enhancers')):
    os.makedirs(os.path.join(outDir, 'candidate_enhancers'))
outDir = '/outs/'
pickle.dump(region_bin_topics_otsu, open(os.path.join(outDir, 'candidate_enhancers/region_bin_topics_otsu.pkl'), 'wb'))
pickle.dump(region_bin_topics_top3k, open(os.path.join(outDir, 'candidate_enhancers/region_bin_topics_top3k.pkl'), 'wb'))
pickle.dump(markers_dict, open(os.path.join(outDir, 'candidate_enhancers/markers_dict.pkl'), 'wb'))
region_bin_topics_otsu = pickle.load(open(os.path.join(outDir, 'candidate_enhancers/region_bin_topics_otsu.pkl'), 'rb'))
region_bin_topics_top3k = pickle.load(open(os.path.join(outDir, 'candidate_enhancers/region_bin_topics_top3k.pkl'), 'rb'))
markers_dict = pickle.load(open(os.path.join(outDir, 'candidate_enhancers/markers_dict.pkl'), 'rb'))

# Save region sets
os.makedirs(os.path.join(outDir, "region_sets"), exist_ok = True)
os.makedirs(os.path.join(outDir, "region_sets", "Topics_otsu"), exist_ok = True)
os.makedirs(os.path.join(outDir, "region_sets", "Topics_top_3k"), exist_ok = True)
os.makedirs(os.path.join(outDir, "region_sets", "DARs_cell_type"), exist_ok = True)
from pycisTopic.utils import region_names_to_coordinates
for topic in region_bin_topics_otsu:
    region_names_to_coordinates(
        region_bin_topics_otsu[topic].index
    ).sort_values(
        ["Chromosome", "Start", "End"]
    ).to_csv(
        os.path.join(outDir, "region_sets", "Topics_otsu", f"{topic}.bed"),
        sep = "\t",
        header = False, index = False
    )
for topic in region_bin_topics_top3k:
    region_names_to_coordinates(
        region_bin_topics_top3k[topic].index
    ).sort_values(
        ["Chromosome", "Start", "End"]
    ).to_csv(
        os.path.join(outDir, "region_sets", "Topics_top_3k", f"{topic}.bed"),
        sep = "\t",
        header = False, index = False
    )
for group in markers_dict:
    region_names_to_coordinates(
        markers_dict[group].index
    ).sort_values(
        ["Chromosome", "Start", "End"]
    ).to_csv(
        os.path.join(outDir, "region_sets", "DARs_cell_type", f"{group}.bed"),
        sep = "\t",
        header = False, index = False
    )
    
# Prepare fasta from consensus regions
export PATH=$PATH:$/outs/genome/bedtools2/bin
BEDTOOLS=/outs/genome/bedtools2/bin/bedtools
# genome.fa = reference file
fasta_filename=/GRCh38-2020-A-build/GRCh38/fasta/genome.fa
consensus_regions=/outs/consensus_peak_calling/consensus_regions.bed
# bed file with peaks in scATAC file
bedtools getfasta \
    -fi $fasta_filename \
    -bed $consensus_regions \
    -fo /outs/genome/consensus_regions_hg38.fa
${create_cistarget_databases_dir}create_fasta_with_padded_bg_from_bed.sh $fasta_filename data/raw/macs-peaks/grouped-peaks.bed  500 yes

#bedtools getfasta -fi ${GENOME_FASTA} -bed ${REGION_BED}
#cat /broad/sankaranlab/cohn/data/GRN/SCENIC/outs/genome/hg38.fa | head

# REGION_BED="/outs/consensus_peak_calling/consensus_regions.bed"
# GENOME_FASTA="/outs/genome/consensus_regions.fa"
# CHROMSIZES=/outs/genome/hg38.chrom.sizes
# DATABASE_PREFIX="1kb_bg_with_mask"
# SCRIPT_DIR="create_cisTarget_databases"
# create_cistarget_databases_dir="create_cisTarget_databases/"
# ${create_cistarget_databases_dir}create_fasta_with_padded_bg_from_bed.sh \
#   $fasta_filename \
#   $consensus_regions \
#   hg38.1kb_bg_padding.fa \ 
#   1000 \
#   yes


# cisTarget motif database output prefix.
# consensdir=/outs/genome
# CBDIR=/aertslab_motif_colleciton/v10nr_clust_public/singletons
# outdir=/outs
# tag='motif'
# motif_list=/outs/genome/motifs.txt
# CLUSTER_BUSTER_PATH=/SCENIC
# export CLUSTER_BUSTER_PATH=SCENIC/cbust

# ${create_cistarget_databases_dir}/create_cistarget_motif_databases.py \
#          -f $consensdir/consensus_regions_hg38.fa \
#          -M $CBDIR \
#          -m $motif_list \
#          -o $outdir/$tag \
#          --cbust $CLUSTER_BUSTER_PATH \
#          -t 30
#          
# OUT_DIR=""${PWD}""
# CBDIR="/aertslab_motif_colleciton/v10nr_clust_public/singletons"
# FASTA_FILE="/outs/genome/consensus_regions.fa"
# MOTIF_LIST="/outs/genome/motifs.txt"
# SCRIPT_DIR="create_cisTarget_databases"
# OUT_DIR="/outs"
# DATABASE_PREFIX="motif"
# "${SCRIPT_DIR}/create_cistarget_motif_databases.py" \
#     -f ${FASTA_FILE} \
#     -M ${CBDIR} \
#     --cbust CLUSTER_BUSTER_PATH
#     -m ${MOTIF_LIST} \
#     -o ${OUT_DIR}/${DATABASE_PREFIX} \
#     -t 10


# SCENIC+ pipline
# fix adata barcodes
# convert scRNA-seq barcodes to scATAC-seq
import scanpy as sc
import pyranges as pr
import os
import pycisTopic
import pickle
adata = sc.read_h5ad(os.path.join('DC.RNA.adata.h5ad'))
adata.obs['barcode'] = [x.split('_')[1] for x in adata.obs.index.tolist()]
# Ensure that 'barcode' and 'replicate' are of string type before concatenation
adata.obs['barcode'] = adata.obs['barcode'].astype(str)
adata.obs['replicate'] = adata.obs['replicate'].astype(str)
adata.obs_names = adata.obs['barcode'] + '-' + adata.obs['replicate'] + '___' + adata.obs['replicate']
output_path = "DC_data_fixed.h5ad"
adata.write(output_path)
cistopic_obj = pickle.load(
            open('/outs/cisTopicObject.pkl', 'rb'))
print(cistopic_obj.cell_data)
print(adata.obs_names)

# run cistopic 
pycistarget cistarget  \
    --bed_fname /outs/region_sets/Topics_top_3k/Topic50.bed \
    --cistarget_db_fname /outs/targets/motif_con.regions_vs_motifs.rankings.feather \
    --path_to_motif_annotations /SCENIC/aertslab_motif_colleciton/v10nr_clust_public/snapshots/motifs-v10-nr.hgnc-m0.00001-o0.0.tbl \
    --output_mode tsv \
    --output_folder /outs/pycistarget/ \
    --species homo_sapiens 

## Regulon specificity calculation
from scenicplus.RSS import *
from scenicplus.eregulon_enrichment import  binarize_AUC
binarize_AUC(scplus_obj,
             auc_key='eRegulon_AUC',
             out_key='eRegulon_AUC_thresholds',
             signature_keys=['Gene_based', 'Region_based'],
             n_cpu=1)

# Snakemake pipline ; described here: https://scenicplus.readthedocs.io/en/latest/human_cerebellum.html


