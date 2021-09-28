import argparse, sys, os
from ipywidgets import widgets
import warnings
warnings.filterwarnings('ignore')
import pandas as pd
import numpy as np
import scanpy as sc

import seaborn as sns
from matplotlib import pyplot as plt


## Parse command-line arguments
# process arguments
parser = argparse.ArgumentParser(description='Drug response prediction')

parser.add_argument('-i', '--input', required=True, help='path to input Anndata object (h5ad file)')
parser.add_argument('-o', '--output', default='./', help='path to output directory, default=\'./\'')
parser.add_argument('-c', '--clusters', default='All', type=str, help='perform IC50 prediction on specified clusters, e.g. \'1,3,8,9\', default=\'All\'')
parser.add_argument('-p', '--platform', default='GDSC', type=str, help='the sensitivity screening is from GDSC ic50/PRISM auc, e.g. GDSC, PRISM')

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

class Drug_Response:
    def __init__(self):
        self.load_model()
        self.drug_info()
        self.bulk_exp()
        self.sc_exp()
        self.kernel_feature_preparartion()
        self.sensitivity_prediction()
        if args.platform == 'GDSC':
            self.masked_drugs = list(pd.read_csv('../data/masked_drugs.csv')['GDSC'].dropna().astype('int64').astype('str'))
            self.cell_death_proportion()
        else:
            self.masked_drugs = list(pd.read_csv('masked_drugs.csv')['PRISM'])
        self.output_result()
        self.figure_output()

    def load_model(self):
        ### IC50/AUC prediction
        ## Read pre-trained model
        model_dir = '/scDrug/CaDRReS-Sc-model/'
        obj_function = widgets.Dropdown(options=['cadrres-wo-sample-bias', 'cadrres-wo-sample-bias-weight'], description='Objetice function')
        self.model_spec_name = obj_function.value
        if args.platform == 'GDSC':
            model_file = model_dir + '{}_param_dict_all_genes.pickle'.format(self.model_spec_name)
        elif args.platform == 'PRISM':
            model_file = model_dir + '{}_param_dict_prism.pickle'.format(self.model_spec_name)
        else:
            sys.exit('Wrong platform name.')
        self.cadrres_model = model.load_model(model_file)
    
    def drug_info(self):
        ## Read drug information
        if args.platform == 'GDSC':
            self.drug_info_df = pd.read_csv(scriptpath + '/preprocessed_data/GDSC/drug_stat.csv', index_col=0)
            self.drug_info_df.index = self.drug_info_df.index.astype(str)
        else:
            self.drug_info_df = pd.read_csv(scriptpath + '/preprocessed_data/PRISM/PRISM_drug_info.csv', index_col='broad_id')
        
    def bulk_exp(self):
        ## Read test data
        if args.platform == 'GDSC':
            self.gene_exp_df = pd.read_csv(scriptpath + '/data/GDSC/GDSC_exp.tsv', sep='\t', index_col=0)
            self.gene_exp_df = self.gene_exp_df.groupby(self.gene_exp_df.index).mean()
        else:
            self.gene_exp_df = pd.read_csv(scriptpath + '/data/CCLE/CCLE_expression.csv', low_memory=False, index_col=0).T
            self.gene_exp_df.index = [gene.split(sep=' (')[0] for gene in self.gene_exp_df.index]

    def sc_exp(self):
        ## Load cluster-specific gene expression profile
        self.adata = sc.read(args.input)
        if args.clusters == 'All':
            clusters = sorted(self.adata.obs['louvain'].unique(), key=int)
        else:
            clusters = [x.strip() for x in args.clusters.split(',')]

        self.cluster_norm_exp_df = pd.DataFrame(columns=clusters, index=self.adata.raw.var.index)
        for cluster in clusters:
            self.cluster_norm_exp_df[cluster] =  self.adata.raw.X[self.adata.obs['louvain']==cluster].mean(axis=0).T \
                                                 if np.sum(self.adata.raw.X[self.adata.obs['louvain']==cluster]) else 0.0

    def kernel_feature_preparartion(self):
        ## Read essential genes list
        if args.platform == 'GDSC':
            ess_gene_list = self.gene_exp_df.index.dropna().tolist()
        else:
            ess_gene_list = utility.get_gene_list(scriptpath + '/preprocessed_data/PRISM/feature_genes.txt')

        ## Calculate fold-change
        cell_line_log2_mean_fc_exp_df, cell_line_mean_exp_df = pp.gexp.normalize_log2_mean_fc(self.gene_exp_df)
            
        self.adata_exp_mean = pd.Series(self.adata.raw.X.mean(axis=0).tolist()[0], index=self.adata.raw.var.index)
        cluster_norm_exp_df = self.cluster_norm_exp_df.sub(self.adata_exp_mean, axis=0)

        ## Calculate kernel feature
        self.test_kernel_df = pp.gexp.calculate_kernel_feature(cluster_norm_exp_df, cell_line_log2_mean_fc_exp_df, ess_gene_list)
    
    def sensitivity_prediction(self):
        ## Drug response prediction
        if args.platform == 'GDSC':
            print('Predicting drug response for using CaDRReS(GDSC): {}'.format(self.model_spec_name))
            self.pred_ic50_df, P_test_df= model.predict_from_model(self.cadrres_model, self.test_kernel_df, self.model_spec_name)
            print('done!')
        else:
            print('Predicting drug response for using CaDRReS(PRISM): {}'.format(self.model_spec_name))
            self.pred_auc_df, P_test_df= model.predict_from_model(self.cadrres_model, self.test_kernel_df, self.model_spec_name)
            print('done!')

    def cell_death_proportion(self):
        ### Drug kill prediction
        ref_type = 'log2_median_ic50'
        self.drug_list = [x for x in self.pred_ic50_df.columns if not x in self.masked_drugs]
        self.drug_info_df = self.drug_info_df.loc[self.drug_list]
        self.pred_ic50_df = self.pred_ic50_df.loc[:,self.drug_list]

        ## Predict cell death percentage at the ref_type dosage
        pred_delta_df = pd.DataFrame(self.pred_ic50_df.values - self.drug_info_df[ref_type].values, columns=self.pred_ic50_df.columns)
        pred_cv_df = 100 / (1 + (np.power(2, -pred_delta_df)))
        self.pred_kill_df = 100 - pred_cv_df
    
    def output_result(self):
        if args.platform == 'GDSC':
            drug_df = pd.DataFrame({'Drug ID': self.drug_list, 
                                    'Drug Name': [self.drug_info_df.loc[drug_id]['Drug Name'] for drug_id in self.drug_list]})
            self.pred_ic50_df.columns = pd.MultiIndex.from_frame(drug_df)
            self.pred_ic50_df.round(3).to_csv(os.path.join(args.output, 'IC50_prediction.csv'))
            self.pred_kill_df.columns = pd.MultiIndex.from_frame(drug_df)
            self.pred_kill_df.round(3).to_csv(os.path.join(args.output, 'drug_kill_prediction.csv'))
        else:
            drug_list = list(self.pred_auc_df.columns)
            drug_list  = [d for d in drug_list if d not in self.masked_drugs]
            drug_df = pd.DataFrame({'Drug ID':drug_list,
                                    'Drug Name':[self.drug_info_df.loc[d, 'name'] for d in drug_list]})
            scaling = pd.Series([240]*1448, index=self.pred_auc_df.columns)
            translation = pd.Series([-120]*1448, index=self.pred_auc_df.columns)
            self.pred_auc_df = (self.pred_auc_df-translation)/scaling
            self.pred_auc_df.columns = pd.MultiIndex.from_frame(drug_df)
            self.pred_auc_df.round(3).to_csv(os.path.join(args.output, 'AUC_prediction.csv'))
    
    def draw_plot(self, df, name='', figsize=()):
        if args.platform == 'GDSC':
            fig, ax = plt.subplots(figsize=figsize) 
            sns.heatmap(df.iloc[:,:-1], cmap='Blues', \
                        linewidths=0.5, linecolor='lightgrey', cbar=True, cbar_kws={'shrink': .2, 'label': name}, ax=ax)
            ax.set(xlabel='Cluster', ylabel='Drug')
            for _, spine in ax.spines.items():
                spine.set_visible(True)
                spine.set_color('lightgrey') 
            plt.savefig(os.path.join(args.output, '{}.png'.format(name)), bbox_inches='tight', dpi=200)
            plt.close()

        else:
            fig, ax = plt.subplots(figsize=figsize) 
            sns.heatmap(df.iloc[:,:-1], cmap='Blues', \
                        linewidths=0.5, linecolor='lightgrey', cbar=True, cbar_kws={'shrink': .2, 'label': name}, ax=ax, vmin=0, vmax=1)
            ax.set(xlabel='Cluster', ylabel='Drug')
            for _, spine in ax.spines.items():
                spine.set_visible(True)
                spine.set_color('lightgrey') 
            plt.savefig(os.path.join(args.output, '{}.png'.format(name)), bbox_inches='tight', dpi=200)
            plt.close()

    def figure_output(self):
        print('Ploting...')
        ## GDSC figures
        if args.platform == 'GDSC':
            tmp_pred_ic50_df = self.pred_ic50_df.iloc[1:,:].T
            tmp_pred_ic50_df = tmp_pred_ic50_df.assign(sum=tmp_pred_ic50_df.sum(axis=1)).sort_values(by='sum', ascending=True)
            self.draw_plot(tmp_pred_ic50_df, name='predicted IC50', figsize=(12,40))
            tmp_pred_kill_df = self.pred_kill_df.iloc[1:,:].T
            tmp_pred_kill_df = tmp_pred_kill_df.loc[(tmp_pred_kill_df>=50).all(axis=1)]
            tmp_pred_kill_df = tmp_pred_kill_df.assign(sum=tmp_pred_kill_df.sum(axis=1)).sort_values(by='sum', ascending=False)
            self.draw_plot(tmp_pred_kill_df, name='predicted cell death', figsize=(12,8))

        ## PRISM figures
        else:
            tmp_pred_auc_df = self.pred_auc_df.iloc[1:,:].T
            tmp_pred_auc_df = tmp_pred_auc_df.assign(sum=tmp_pred_auc_df.sum(axis=1)).sort_values(by='sum', ascending=True)
            self.draw_plot(tmp_pred_auc_df, name='predicted AUC', figsize=(12,60))  
        print('done!') 


job = Drug_Response()
