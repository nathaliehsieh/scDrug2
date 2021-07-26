import numpy as np
import pandas as pd
import scanpy as sc

results_file = '../write/GSE156625_HCCscanpyobj.h5ad'

adata = sc.read_10x_mtx(
        '../data/GSE156625_HCC/',
        var_names='gene_symbols',
        cache=True)
adata.var_names_make_unique()


## Preprocessing
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)

adata.var['mt'] = adata.var_names.str.startswith('MT-')  # annotate the group of mitochondrial genes as 'mt'
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

adata = adata[adata.obs.pct_counts_mt < 30, :]

sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

sc.pp.highly_variable_genes(adata)
# sc.pl.highly_variable_genes(adata)

adata.raw = adata
adata = adata[:, adata.var.highly_variable]
sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'])

sc.pp.scale(adata)


## Principal component analysis
sc.tl.pca(adata, svd_solver='arpack')
# sc.pl.pca_variance_ratio(adata, log=True)

adata.write(results_file)


## Batch Correction with Harmony
section_df = pd.read_csv('../data/section_info.csv')
patient_id = [int(barcode.split('-')[-1]) for barcode in adata.obs.index]
adata.obs['PatientID'] = [section_df['PatientID'][pid-1] for pid in patient_id]
adata.obs['NormalvsTumor'] = [section_df['NormalvsTumor'][pid-1] for pid in patient_id]
sc.external.pp.harmony_integrate(adata, 'PatientID', adjusted_basis='X_pca')

## Computing the neighborhood graph
sc.pp.neighbors(adata, n_pcs=20)

## Embedding the neighborhood graph
sc.tl.umap(adata)

## Clustering the neighborhood graph
sc.tl.louvain(adata, resolution=0.6)
# sc.pl.umap(adata, color=['louvain'], use_raw=False)
# sc.pl.umap(adata, color=['PatientID'], use_raw=False)

adata.write(results_file)

## Finding differentially expressed genes
method = "t-test"
sc.tl.rank_genes_groups(adata, 'louvain', method=method)
# sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False)


## scMatch
# Export csv used by scMatch
result = adata.uns['rank_genes_groups']
groups = result['names'].dtype.names
mat = np.zeros((len(adata.raw.var.index), len(groups)), dtype=float)
for group in groups:
    mat[: , int(group)] = adata.raw.X[adata.obs['louvain']==group].mean(axis=0)
dat = pd.DataFrame(mat, index = adata.raw.var.index, columns = groups)
dat.to_csv("../data/harmony_scanpy_cluster_" + method + "_GSE156625_HCC_scMatch.csv")

# sh scMatch_cluster.sh

# Cell annotation result
scMatch_cluster_df = pd.read_csv('../data/harmony_scanpy_cluster_' + method + '_GSE156625_HCC_scMatch/annotation_result_keep_all_genes/human_Spearman_top_ann.csv')
scMatch_cluster_names = [group + "_" + scMatch_cluster_df.loc[scMatch_cluster_df['cell']==int(group)]\
                          ['cell type'].tolist()[0] for group in groups]
# adata.rename_categories('louvain', scMatch_cluster_names)
# sc.pl.umap(adata, color=['louvain'], legend_loc='on data', use_raw=False)


## Specify Tumor cluster
total_tumor_n = sum(adata.obs['NormalvsTumor'])
total_normal_n = len(adata.obs) - total_tumor_n
cluster_tumor_percentage = np.zeros(len(groups), dtype = float)
cluster_normal_percentage = np.zeros(len(groups), dtype = float)
for group in groups:
    cells = adata.obs['louvain'] == group
    cluster_tumor_n = sum(adata.obs['NormalvsTumor'][cells])
    cluster_tumor_percentage[int(group)] = cluster_tumor_n / total_tumor_n
    cluster_normal_percentage[int(group)] = (sum(cells) - cluster_tumor_n) / total_normal_n

# import matplotlib.pyplot as plt
# plt.scatter(cluster_normal_percentage, cluster_tumor_percentage) 
# plt.show()

tumor_cluster = np.where(cluster_tumor_percentage/cluster_normal_percentage > 2)[0].astype(str)
tumor_adata = adata.raw.to_adata()[adata.obs['louvain'].isin(tumor_cluster)]

sc.pp.highly_variable_genes(tumor_adata)
tumor_adata = tumor_adata[:, tumor_adata.var.highly_variable]
sc.pp.regress_out(tumor_adata, ['total_counts', 'pct_counts_mt'])

sc.pp.scale(tumor_adata)

# Principal component analysis
sc.tl.pca(tumor_adata, svd_solver='arpack')
# sc.pl.pca_variance_ratio(tumor_adata, log=True)

# Batch Correction with Harmony
section_df = pd.read_csv('../data/section_info.csv')
patient_id = [int(barcode.split('-')[-1]) for barcode in tumor_adata.obs.index]
tumor_adata.obs['PatientID'] = [section_df['PatientID'][pid-1] for pid in patient_id]
tumor_adata.obs['NormalvsTumor'] = [section_df['NormalvsTumor'][pid-1] for pid in patient_id]
# sc.external.pp.harmony_integrate(tumor_adata, 'PatientID')
sc.external.pp.harmony_integrate(tumor_adata, 'PatientID', adjusted_basis='X_pca')

# Computing the neighborhood graph
sc.pp.neighbors(tumor_adata, n_pcs=20)

# Embedding the neighborhood graph
sc.tl.umap(tumor_adata)


## SCCAF
from SCCAF import *

sc.tl.louvain(tumor_adata, resolution=2, key_added = 'L2_Round0')
sc.pl.umap(tumor_adata, color=['L2_Round0'], use_raw=False)

# SCCAF recover
SCCAF_optimize_all(ad=tumor_adata, plot=False, min_acc=0.9,  prefix = 'L2', use='pca')
sc.pl.umap(tumor_adata, color=['L2_result'], use_raw=False)

# Calculate silhouette score
from sklearn.metrics import silhouette_samples, silhouette_score
import matplotlib.pyplot as plt

X = tumor_adata.obsm['X_pca'][:, :20]
SCCAF_silhouette_avg = silhouette_score(X, tumor_adata.obs['L2_result'])
# SCCAF_sample_silhouette_values = silhouette_samples(tumor_adata.obsm['X_pca'], tumor_adata.obs['L2_result'])
# SCCAF_cluster_silhouette_values = [np.mean(SCCAF_sample_silhouette_values[np.where(tumor_adata.obs['louvain']==group)])
#                                    for group in tumor_adata.obs['louvain'].unique()]

## Louvain with silhouette score
R = np.arange(0.1, 1.6, 0.1)
silhouette_avg = []
sample_silhouette_values = []
for res in R: 
    sc.tl.louvain(tumor_adata, resolution=res)
    # sc.pl.umap(tumor_adata, color=['louvain'], title='louvain, res=' + str(np.round(res, 1)))
    
    ## Calculate silhouette score
    silhouette_avg.append(silhouette_score(X, tumor_adata.obs['louvain']))
    # sample_silhouette_values.append(silhouette_samples(X, tumor_adata.obs['louvain']))
    # cluster_silhouette_values = [np.mean(sample_silhouette_values[0][np.where(tumor_adata.obs['louvain']==group)])
    #                               for group in tumor_adata.obs['louvain'].unique()]
best_res = R[np.argmax(silhouette_avg)]


## TODO: Determine sub-clustering result

## Finding differentially expressed genes
method = "t-test"
sc.tl.rank_genes_groups(tumor_adata, 'louvain', method=method)
# sc.pl.rank_genes_groups(tumor_adata, n_genes=25, sharey=False)


## CaDRRes-SC
result = tumor_adata.uns['rank_genes_groups']
cluster = '0'
dat = pd.DataFrame({'cluster' + cluster + '_log_fc' : result['logfoldchanges'][cluster]}, index = result['names'][cluster])
dat.to_csv("../data/GSE156625_HCC_tumor_cluster' + cluster + '_log_fc.csv")

# python3 CaDRRes-SC_predicting_cluster_drug_response.py
