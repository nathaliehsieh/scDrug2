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
parser = argparse.ArgumentParser(description='Drug response prediction')

parser.add_argument('-i', '--input', required=True, help='path to input Anndata object (h5ad file)')
parser.add_argument('-o', '--output', default='./', help='path to output directory, default='./'')
parser.add_argument('-c', '--clusters', default='All', type=str, help='perform IC50 prediction on specified clusters, e.g. '1,3,8,9', default='All'')

args = parser.parse_args()

# check input
if not os.path.exists(args.input):
    sys.exit('The input path does not exist.')
if args.input[-5:] != '.h5ad':
    sys.exit('The input file is not a h5ad file.')

# check output
if not os.path.isdir(args.output):
    sys.exit('The output directory does not exist.')

scriptpath = '/opt/CaDRReS-Sc'
sys.path.append(os.path.abspath(scriptpath))

from cadrres_sc import pp, model, evaluation, utility

## Read drug statistics
drug_info_df = pd.read_csv(scriptpath + '/preprocessed_data/GDSC/drug_stat.csv', index_col=0)
drug_info_df.index = drug_info_df.index.astype(str)

### IC50 prediction
## Read pre-trained model
model_dir = '/scDrug/CaDRReS-Sc-model/'
obj_function = widgets.Dropdown(options=['cadrres-wo-sample-bias', 'cadrres-wo-sample-bias-weight'], description='Objetice function')
model_spec_name = obj_function.value
model_file = model_dir + '{}_param_dict_all_genes.pickle'.format(model_spec_name)
cadrres_model = model.load_model(model_file)

## Read test data
gene_exp_df = pd.read_csv(scriptpath + '/data/GDSC/GDSC_exp.tsv', sep='\t', index_col=0)
gene_exp_df = gene_exp_df.groupby(gene_exp_df.index).mean()

## Load cluster-specific gene expression profile
adata = sc.read(args.input)
if args.clusters == 'All':
    clusters = sorted(adata.obs['louvain'].unique(), key=int)
else:
    clusters = [x.strip() for x in args.clusters.split(',')]

cluster_norm_exp_df = pd.DataFrame(columns=clusters, index=adata.raw.var.index)
for cluster in clusters:
    cluster_norm_exp_df[cluster] =  adata.raw.X[adata.obs['louvain']==cluster].mean(axis=0).T \
                                    if np.sum(adata.raw.X[adata.obs['louvain']==cluster]) else 0.0

## Read essential genes list
ess_gene_list = gene_exp_df.index.dropna().tolist()

## Calculate fold-change
cell_line_log2_mean_fc_exp_df, cell_line_mean_exp_df = pp.gexp.normalize_log2_mean_fc(gene_exp_df)
    
adata_exp_mean = pd.Series(adata.raw.X.mean(axis=0).tolist()[0], index=adata.raw.var.index)
cluster_norm_exp_df = cluster_norm_exp_df.sub(adata_exp_mean, axis=0)

## Calculate kernel feature
test_kernel_df = pp.gexp.calculate_kernel_feature(cluster_norm_exp_df, cell_line_log2_mean_fc_exp_df, ess_gene_list)

## Drug response prediction
print('Predicting drug response for using CaDRReS: {}'.format(model_spec_name))
pred_ic50_df, P_test_df= model.predict_from_model(cadrres_model, test_kernel_df, model_spec_name)
print('done!')

### Drug kill prediction
ref_type = 'log2_median_ic50'
masked_drugs = ['293','1062','193','255','119','166','147','1038','202','37','1133','136','35','86','34','170','1069','156','71','207','88','185','180','1053','1066','165','52','63','186','1023','172','17','1058','59','163','94','1042','127','89','106','1129','6','1067','199','64','1029','111','1072','192','1009','104','1039','1043','110','91']
drug_list = [ x for x in pred_ic50_df.columns if x not in masked_drugs]
drug_info_df = drug_info_df.loc[drug_list]

## Predict cell death percentage at the ref_type dosage
pred_delta_df = pd.DataFrame(pred_ic50_df.values - drug_info_df[ref_type].values, columns=drug_list)
pred_cv_df = 100 / (1 + (np.power(2, -pred_delta_df)))
pred_kill_df = 100 - pred_cv_df

drug_df = pd.DataFrame({'Drug ID': drug_list, 
                        'Drug Name': [drug_info_df.loc[drug_id]['Drug Name'] for drug_id in drug_list]})
pred_ic50_df.columns = pd.MultiIndex.from_frame(drug_df)
pred_ic50_df.round(3).to_csv(os.path.join(args.output, 'IC50_prediction.csv'))
pred_kill_df.columns = pd.MultiIndex.from_frame(drug_df)
pred_kill_df.round(3).to_csv(os.path.join(args.output, 'drug_kill_prediction.csv'))


