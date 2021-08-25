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
parser = argparse.ArgumentParser(description="Drug response prediction")

parser.add_argument("-i", "--input", required=True, help="path to input Anndata object (h5ad file)")
parser.add_argument("-o", "--output", default='./', help="path to output directory, default='./'")
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


scriptpath = '/opt/CaDRReS-Sc'
sys.path.append(os.path.abspath(scriptpath))

from cadrres_sc import pp, model, evaluation, utility

### IC50 prediction
## Read pre-trained model
model_dir = '/single-cell-analysis/CaDRReS-Sc-model/'
obj_function = widgets.Dropdown(options=['cadrres-wo-sample-bias', 'cadrres-wo-sample-bias-weight'], description='Objetice function')
model_spec_name = obj_function.value
model_file = model_dir + '{}_param_dict_all_genes.pickle'.format(model_spec_name)
cadrres_model = model.load_model(model_file)

## Read test data
gene_exp_df = pd.read_csv(scriptpath + '/data/GDSC/GDSC_exp.tsv', sep='\t', index_col=0)
gene_exp_df = gene_exp_df.groupby(gene_exp_df.index).mean()

## Load cluster-specific gene expression profile
adata = sc.read(args.input)
adata_exp_df = pd.DataFrame(adata.raw.X.transpose().toarray(), index=adata.raw.var.index, columns=adata.obs.index)

## Read essential genes list
ess_gene_list = adata.var.highly_variable.index if args.hvg else gene_exp_df.index.dropna().tolist()

## Calculate fold-change
cell_line_log2_mean_fc_exp_df, cell_line_mean_exp_df = pp.gexp.normalize_log2_mean_fc(gene_exp_df)
adata_log2_mean_fc_exp_df, adata_mean_exp_df = pp.gexp.normalize_log2_mean_fc(adata_exp_df)

## Calculate kernel feature
test_kernel_df = pp.gexp.calculate_kernel_feature(adata_log2_mean_fc_exp_df, cell_line_log2_mean_fc_exp_df, ess_gene_list)

## Drug response prediction
print('Predicting drug response for using CaDRReS: {}'.format(model_spec_name))
cell_pred_df, P_test_df= model.predict_from_model(cadrres_model, test_kernel_df, model_spec_name)
print('done!')

pred_ic50_df = cell_pred_df.mean().to_frame().T


### Drug kill prediction
## Read drug statistics
drug_info_df = pd.read_csv(scriptpath + '/preprocessed_data/GDSC/drug_stat.csv', index_col=0)
drug_info_df.index = drug_info_df.index.astype(str)
ref_type = 'log2_median_ic50'
drug_list = pred_ic50_df.columns

## Predict cell death percentage at the ref_type dosage
pred_delta_df = pd.DataFrame(pred_ic50_df.values - drug_info_df[ref_type].values, columns=drug_list)
pred_cv_df = 100 / (1 + (np.power(2, -pred_delta_df)))
pred_kill_df = 100 - pred_cv_df

pred_df = pd.DataFrame({'Drug ID': pred_ic50_df.columns,
                        'Drug Name': [drug_info_df.loc[drug_id]['Drug Name'] for drug_id in pred_ic50_df.columns],
                        'IC50 Prediction': cell_pred_df.mean().tolist(), 
                        'Drug Kill Prediction': pred_kill_df.loc[0].tolist()})

pred_df.round(3).to_csv(os.path.join(args.output, 'drug_response_prediction.csv'), index=False)


