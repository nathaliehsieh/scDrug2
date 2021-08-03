import argparse, os, sys
import numpy as np
import pandas as pd
import scanpy as sc

def runGSEAPY(adata, group_by='louvain', gene_sets=['GO_Biological_Process_2021'], organism='Human', cutoff=0.05, logfc_threshold=2):
    import gseapy as gp

    df_list = []
    cluster_list = []
    for celltype in set(adata.obs[group_by]):
        indlist_logfc = adata.uns['rank_genes_groups']['logfoldchanges'][celltype] >= logfc_threshold
        indlist_adjp = adata.uns['rank_genes_groups']['pvals_adj'][celltype] <= 1e-2
        indlist_p = adata.uns['rank_genes_groups']['pvals'][celltype] <= 1e-2
        #indlist_pts = adata.uns['rank_genes_groups']['pts'][celltype] >= 0.1
        
        indlist = indlist_logfc * indlist_adjp * indlist_p 

        ind = [x for x in range(0, len(indlist)) if indlist[x] ]
        degs = adata.uns['rank_genes_groups']['names'][celltype][ind].tolist()
        
        if not degs:
            continue

        enr = gp.enrichr(gene_list=degs,
                gene_sets=gene_sets,
                organism=organism, 
                description=celltype,
                no_plot=True
                )
        df_list.append(enr.res2d)
        cluster_list.append(celltype)
    
    columns = ['Cluster', 'Gene_set', 'Term', 'Overlap', 'P-value', 'Adjusted P-value', 'Genes']

    df = pd.DataFrame(columns = columns)
    for cluster_ind, df_ in zip(cluster_list, df_list):
        df_ = df_[df_['Adjusted P-value'] <= cutoff]
        df_ = df_.assign(Cluster = cluster_ind)
        if(df_.shape[0] > 0):
            df = pd.concat([df, df_[columns]], sort=False)
        else:
            print('No pathway with an adjusted P-value less than the cutoff (={}) for cluster {}'.format(cutoff, cluster_ind))
    
    return df


## Parse command-line arguments
# process arguments
parser = argparse.ArgumentParser(description="scRNA-seq data analysis")

parser.add_argument("-i", "--input", required=True, help="path to input 10x directory or CSV file")
parser.add_argument("-f", "--format", default='10x', help="input format, 10x (default) | csv | h5ad (Anndata object for subclustering with --clusters CLUSTERS)")
parser.add_argument("-o", "--output", default='./', help="path to output directory, default='./'")
parser.add_argument("-r", "--resolution", type=float, default=0.6, help="resolution for clustering, default=0.6")
# parser.add_argument("--auto-resolution", action="store_true", help="set True to automatically determine resolution for clustering, default=False")
parser.add_argument("-m", "--metadata", default=None, help="path to metadata CSV file for batch correction (index as input in first column)")
parser.add_argument("-b", "--batch", default=None, help="column in metadata (or adata.obs) for batch correction, e.g. 'PatientID'")
parser.add_argument("-c", "--clusters", default=None, help="perform single cell analysis only on specified clusters, e.g. '1,3,8,9'")
parser.add_argument("-a", "--annotation", action="store_true", help="perform cell type annotation")
parser.add_argument("-g", "--gsea", action="store_true", help="perform gene set enrichment analsis (GSEA)")

args = parser.parse_args()

# check format, input and clusters
if not os.path.exists(args.input):
    sys.exit("The input path does not exist.")
if args.format == 'csv':
    if args.input[-4:] != '.csv':
        sys.exit("The input file is not a CSV file.")
elif args.format == '10x':
    if not os.path.exists(os.path.join(args.input, 'matrix.mtx')) and not os.path.exists(os.path.join(args.input, 'matrix.mtx.gz')):
        sys.exit("Cannot find 'matrix.mtx' file in the input directory.")
    if not os.path.exists(os.path.join(args.input, 'genes.tsv')) and not os.path.exists(os.path.join(args.input, 'genes.tsv.gz')):
        if not os.path.exists(os.path.join(args.input, 'features.tsv')) and not os.path.exists(os.path.join(args.input, 'features.tsv.gz')):
            sys.exit("Cannot find 'genes.tsv' or 'features.tsv' file in the input directory.")
    if not os.path.exists(os.path.join(args.input, 'barcodes.tsv')) and not os.path.exists(os.path.join(args.input, 'barcodes.tsv.gz')):
        sys.exit("Cannot find 'barcodes.tsv' file in the input directory.")
elif args.format == 'h5ad':
    if args.input[-5:] != '.h5ad':
        sys.exit("The input file is not a h5ad file.")
    if args.clusters is None:
        sys.exit("Need to speficy clusters to be analyzed with a h5ad file.")
else:
     sys.exit("The format can only be '10x' or 'csv'.")

# check output
if not os.path.isdir(args.output):
    sys.exit("The output directory does not exist.")

# check metadata
if not args.metadata is None:
    if not os.path.exists(args.metadata):
        sys.exit("The metadata file does not exist.")
    if args.metadata[-4:] != '.csv':
        sys.exit("The metadata file is not a CSV file.")

# check batch in metadata
if not args.batch is None:
    if args.metadata is None and  args.format != 'h5ad':
        sys.exit("Please provide the metadata file for batch correction with --metadata METADATA.")
    if not args.metadata is None:
        metadata_df = pd.read_csv(args.metadata, index_col=0)
        if not args.batch in metadata_df.columns:
            sys.exit("The batch column is not in the metadata file.")


## Preprocessing
results_file = os.path.join(args.output, 'scanpyobj.h5ad')

if args.format == 'h5ad':
    adata = sc.read(args.input)
    clusters = [x.strip() for x in args.clusters.split(',')]
    adata = adata.raw.to_adata()[adata.obs['louvain'].isin(clusters)]
else:
    if args.format == 'csv':
        adata = sc.read_csv(args.input)
    elif args.format == '10x':
        adata = sc.read_10x_mtx(args.input, var_names='gene_symbols', cache=True)
    adata.var_names_make_unique()

    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=3)
    
    adata.var['mt'] = adata.var_names.str.startswith('MT-')
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
    
    adata = adata[adata.obs.pct_counts_mt < 30, :]

# GEP preperation
adata_GEP = adata.copy()
sc.pp.normalize_total(adata_GEP, target_sum=1e6)
mat = adata_GEP.X.transpose()
GEP_df = pd.DataFrame(mat.toarray(), index=adata_GEP.var.index)

sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

sc.pp.highly_variable_genes(adata)
adata.raw = adata
adata = adata[:, adata.var.highly_variable]
sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'])

sc.pp.scale(adata)


## Principal component analysis
sc.tl.pca(adata, svd_solver='arpack')

# adata.write(results_file)


## Batch Correction with Harmony
if not args.batch is None:
    if not args.metadata is None:
        adata.obs[args.batch] = metadata_df.loc[adata.obs.index][args.batch]
    elif args.format == 'h5ad' and not args.batch in adata.obs.columns:
        sys.exit("The batch column is not in the Anndata object.")
    
    sc.external.pp.harmony_integrate(adata, args.batch, adjusted_basis='X_pca')


## Clustering
sc.pp.neighbors(adata, n_pcs=20)
sc.tl.umap(adata)
sc.tl.louvain(adata, resolution=args.resolution)

# adata.write(results_file)


groups = sorted(adata.obs['louvain'].unique(), key=int)
if args.annotation:
    ## scMatch
    # Export csv used by scMatch
    mat = np.zeros((len(adata.raw.var.index), len(groups)), dtype=float)
    for group in groups:
        mat[: , int(group)] = adata.raw.X[adata.obs['louvain']==group].mean(axis=0)
    dat = pd.DataFrame(mat, index = adata.raw.var.index, columns = groups)
    dat.to_csv(os.path.join(args.output, 'cluster_mean_exp.csv'))
    
    os.system('python /opt/scMatch/scMatch.py --refDS /opt/scMatch/refDB/FANTOM5 \
              --dFormat csv --testDS ' + os.path.join(args.output, 'cluster_mean_exp.csv'))
    
    # Cell annotation result
    scMatch_cluster_df = pd.read_csv(os.path.join(args.output, 'cluster_mean_exp') + '/annotation_result_keep_all_genes/human_Spearman_top_ann.csv')
    scMatch_cluster_names = [group + " " + scMatch_cluster_df.loc[scMatch_cluster_df['cell']==int(group)]\
                              ['cell type'].tolist()[0] for group in groups]
    adata.obs['cell_type'] = adata.obs['louvain'].cat.rename_categories(scMatch_cluster_names)

sc.settings.autosave = True
sc.settings.figdir = args.output
sc.pl.umap(adata, color=['louvain'], use_raw=False, show=False, save='_cluster.png')
if not args.batch is None:
    sc.pl.umap(adata, color=[args.batch], use_raw=False, show=False, save='_batch.png')
if args.annotation:
    sc.pl.umap(adata, color=['cell_type'], use_raw=False, show=False, save='_cell_type.png')

# GEP format
GEP_df.columns = adata.obs['louvain'].tolist()
# GEP_df = GEP_df.loc[adata.var.index[adata.var.highly_variable==True]]
GEP_df.dropna(axis=1, inplace=True)
GEP_df.to_csv(os.path.join(args.output, 'GEP.txt'), sep='\t')


## Finding differentially expressed genes
method = "t-test"
sc.tl.rank_genes_groups(adata, 'louvain', method=method, pts=True)

adata.write(results_file)

# cluster DEGs
result = adata.uns['rank_genes_groups']
dat = pd.DataFrame({group + '_' + key: result[key][group] for group in groups for key in ['names', 'logfoldchanges','scores','pvals']})
dat.to_csv(os.path.join(args.output, 'cluster_DEGs.csv'))

# perform GSEA
if args.gsea:
    df_gsea = runGSEAPY(adata)
    df_gsea.to_csv(os.path.join(args.output, 'GSEA_results.csv'))


