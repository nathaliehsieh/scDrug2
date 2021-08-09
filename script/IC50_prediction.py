import argparse, sys, os, pickle
from collections import Counter
import importlib
from ipywidgets import widgets
import warnings
warnings.filterwarnings('ignore')
import pandas as pd
import numpy as np
import scanpy as sc


## Parse command-line arguments
# process arguments
parser = argparse.ArgumentParser(description="IC50 prediction")

parser.add_argument("-i", "--input", required=True, help="path to input Anndata object (h5ad file)")
parser.add_argument("-o", "--output", default='./', help="path to output directory, default='./'")
parser.add_argument("-c", "--clusters", default='All', type=str, help="perform IC50 prediction on specified clusters, e.g. '1,3,8,9', default='All'")
parser.add_argument("-m", "--method", default='average', help="method used to produce clusterwise IC50 prediction, options: average | logfc, default: average")
parser.add_argument("--hvg", action="store_true", help="only use highly variable genes to predict IC50")

args = parser.parse_args()

# check input
if not os.path.exists(args.input):
    sys.exit("The input path does not exist.")
if args.input[-5:] != '.h5ad':
    sys.exit("The input file is not a h5ad file.")

# check output
if not os.path.isdir(args.output):
    sys.exit("The output directory does not exist.")

# check method
if not args.method in ['average', 'logfc']:
    sys.exit("The method can only be 'average' or 'logfc'.")

scriptpath = '../CaDRReS-Sc-master'
sys.path.append(os.path.abspath(scriptpath))

from cadrres_sc import pp, model, evaluation, utility

## Read pre-trained model
model_dir = '../CaDRReS-Sc-model/'
obj_function = widgets.Dropdown(options=['cadrres-wo-sample-bias', 'cadrres-wo-sample-bias-weight'], description='Objetice function')
model_spec_name = obj_function.value
model_file = model_dir + '{}_param_dict_all_genes.pickle'.format(model_spec_name)
cadrres_model = model.load_model(model_file)

## Read test data
gene_exp_df = pd.read_csv(scriptpath + '/data/GDSC/GDSC_exp.tsv', sep='\t', index_col=0)
gene_exp_df = gene_exp_df.groupby(gene_exp_df.index).mean()

# Calculate fold-change
cell_line_log2_mean_fc_exp_df, cell_line_mean_exp_df = pp.gexp.normalize_log2_mean_fc(gene_exp_df)

## Load cluster-specific gene expression profile
adata = sc.read(args.input)
if args.clusters == 'All':
    clusters = sorted(adata.obs['louvain'].unique(), key=int)
else:
    clusters = [x.strip() for x in args.clusters.split(',')]

## Read essential genes list
ess_gene_list = adata.var.highly_variable.index if args.hvg else gene_exp_df.index.dropna().tolist()

if args.method == 'average':
    adata_exp_mean = pd.Series(adata.raw.X.mean(axis=0).tolist()[0], index=adata.raw.var.index)
    clusters_pred_df = pd.DataFrame()
    for cluster in clusters:
        cluster_adata = adata[adata.obs['louvain']==cluster]
        cluster_exp_df = pd.DataFrame(cluster_adata.raw.X.transpose().toarray(), index=cluster_adata.raw.var.index, columns=cluster_adata.obs.index)
        
        ## Calculate fold-change
        cluster_log2_mean_fc_exp_df = cluster_exp_df.sub(adata_exp_mean, axis=0)
        
        ## Calculate kernel feature
        test_kernel_df = pp.gexp.calculate_kernel_feature(cluster_log2_mean_fc_exp_df, cell_line_log2_mean_fc_exp_df, ess_gene_list)
        
        ## Drug response prediction
        print('Predicting drug response for cluster {} using CaDRReS: {}'.format(cluster, model_spec_name))
        pred_df, P_test_df= model.predict_from_model(cadrres_model, test_kernel_df, model_spec_name)
        print('done!')
        
        pred_mean_df = pred_df.mean()
        drug_df = pd.read_csv(scriptpath + '/preprocessed_data/GDSC/drug_stat.csv', index_col=0)
        pred_mean_df.index = [drug_df.loc[int(drug_id)]['Drug Name'] for drug_id in pred_mean_df.index]
        pred_mean_df = pred_mean_df.groupby(pred_mean_df.index).mean()
        # pred_mean_df = pred_mean_df.sort_values()
        pred_mean_df = pred_mean_df.to_frame().transpose()
        pred_mean_df.index = [cluster]
        clusters_pred_df = clusters_pred_df.append(pred_mean_df)
    clusters_pred_df.to_csv(os.path.join(args.output, 'IC50_prediction.csv'))

elif args.method == 'logfc':
    cluster_norm_exp_df = pd.DataFrame(columns=clusters, index=adata.raw.var.index)
    for cluster in clusters:
        log_fc_df = pd.DataFrame({cluster: adata.uns['rank_genes_groups']['logfoldchanges'][cluster]}, 
                                  index=adata.uns['rank_genes_groups']['names'][cluster])
        cluster_norm_exp_df[cluster] = log_fc_df.loc[cluster_norm_exp_df.index]
    
    ## Calculate kernel feature
    test_kernel_df = pp.gexp.calculate_kernel_feature(cluster_norm_exp_df, cell_line_log2_mean_fc_exp_df, ess_gene_list)
    
    ## Drug response prediction
    print('Predicting drug response for cluster {} using CaDRReS: {}'.format(cluster, model_spec_name))
    pred_df, P_test_df= model.predict_from_model(cadrres_model, test_kernel_df, model_spec_name)
    print('done!')
    
    drug_df = pd.read_csv(scriptpath + '/preprocessed_data/GDSC/drug_stat.csv', index_col=0)
    pred_df.columns = [drug_df.loc[int(drug_id)]['Drug Name'] for drug_id in pred_df.columns]
    pred_df = pred_df.groupby(pred_df.columns, axis=1).mean()
    pred_df.to_csv(os.path.join(args.output, 'IC50_prediction.csv'))

