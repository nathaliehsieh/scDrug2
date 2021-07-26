import sys, os, pickle
from collections import Counter
import importlib
from ipywidgets import widgets
import warnings
warnings.filterwarnings('ignore')
import pandas as pd
import numpy as np

scriptpath = '../CaDRReS-Sc-master'
sys.path.append(os.path.abspath(scriptpath))
pd.set_option('precision', 2)

from cadrres_sc import pp, model, evaluation, utility

## Read pre-trained model
model_dir = '../CaDRReS-Sc-model/'
obj_function = widgets.Dropdown(options=['cadrres-wo-sample-bias', 'cadrres-wo-sample-bias-weight'], description='Objetice function')

## Load the pre-trained model based on your selection
model_name = obj_function.value
model_suffix = 'all_genes'
model_file = model_dir + '{}_param_dict_{}.pickle'.format(model_name, model_suffix)
cadrres_model = model.load_model(model_file)

## Read test data
gene_exp_df = pd.read_csv('../CaDRReS-Sc-master/data/GDSC/GDSC_exp.tsv', sep='\t', index_col=0)
gene_exp_df = gene_exp_df.groupby(gene_exp_df.index).mean()
print("Dataframe shape:", gene_exp_df.shape, "\n")

## Calculate fold-change
cell_line_log2_mean_fc_exp_df, cell_line_mean_exp_df = pp.gexp.normalize_log2_mean_fc(gene_exp_df)

## Load cluster-specific gene expression profile
output_dir = '../write/'

# log fc data
cluster_norm_exp_fname = '../data/GSE156625_HCC_cluster_14_log_fc.csv'
cluster_norm_exp_df = pd.read_csv(cluster_norm_exp_fname, index_col=0)

## Read essential genes list
# ess_gene_list = utility.get_gene_list('../../data/essential_genes.txt')
ess_gene_list = gene_exp_df.index.dropna().tolist()
selected_gene_list = [g for g in ess_gene_list if g in cluster_norm_exp_df.index]
len(selected_gene_list)

## Calculate kernel feature
test_kernel_df = pp.gexp.calculate_kernel_feature(cluster_norm_exp_df, cell_line_log2_mean_fc_exp_df, selected_gene_list)
print("Dataframe shape:", test_kernel_df.shape, "\n")

## Drug response prediction
print('Predicting drug response using CaDRReS: {}'.format(model_name))
pred_df, P_test_df= model.predict_from_model(cadrres_model, test_kernel_df)
print('done!')

# drug ID to drug name
# drug_df = pd.read_csv('../CaDRReS-Sc-master/preprocessed_data/GDSC/drug_stat.csv', index_col=0)
# drug_name = [drug_df.loc[int(drug_id)]['Drug Name'] for drug_id in pred_df.columns]
# pred_df.columns = drug_name

# cluster vs drugs
print('Saving ' + model_dir + '{}_test_pred.csv'.format(model_name))
pred_df.to_csv(output_dir + '{}_cluster_14_pred_by_log_fc.csv'.format(model_name))

## Drug discovery
drug_df = pd.read_csv('../CaDRReS-Sc-master/preprocessed_data/GDSC/drug_stat.csv', index_col=0)

## IC50 value
drugs_by_value = pred_df.loc['14_l'].sort_values()[:10]
drug_name = [drug_df.loc[int(drug_id)]['Drug Name'] for drug_id in drugs_by_value.index]
drugs_by_value.index = drug_name

## IC50 z-score
# Read data
cell_line_obs_df = pd.read_csv('../CaDRReS-Sc-master/data/GDSC/gdsc_all_abs_ic50_bayesian_sigmoid_only9dosages.csv', index_col=0)
cell_line_obs_df.index = cell_line_obs_df.index.astype(str)
# drug_name = [drug_df.loc[int(drug_id)]['Drug Name'] for drug_id in cell_line_obs_df.columns]
# cell_line_obs_df.columns = drug_name

pred_df_z = pd.DataFrame()
for drug in pred_df.columns:
    pred_df_z[drug] = (pred_df[drug] - np.mean(cell_line_obs_df[drug])) / np.std(cell_line_obs_df[drug])

drugs_by_z = pred_df_z.loc['14_l'].sort_values()[:10]
drug_name = [drug_df.loc[int(drug_id)]['Drug Name'] for drug_id in drugs_by_z.index]
drugs_by_z.index = drug_name

