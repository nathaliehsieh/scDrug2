import argparse, sys, os, pickle
from collections import Counter
import importlib
from ipywidgets import widgets
import warnings
warnings.filterwarnings('ignore')
import pandas as pd
import numpy as np
import scanpy as sc

# Parse command-line arguments
# process arguments
parser = argparse.ArgumentParser(description="Drug kill prediction")

parser.add_argument("-i", "--input", required=True, help="path to cluster-specific IC50 prediction (CSV file)")
parser.add_argument("-p", "--proportion", required=True, help="path to cluster proportion CSV file")
parser.add_argument("-o", "--output", default='./', help="path to output directory, default='./'")

args = parser.parse_args()

# check input
if not os.path.exists(args.input):
    sys.exit("The IC50 prediction file does not exist.")
if args.input[-4:] != '.csv':
    sys.exit("The IC50 prediction file is not a CSV file.")

# check proportion
if not os.path.exists(args.proportion):
    sys.exit("The proportion file does not exist.")
if args.input[-4:] != '.csv':
    sys.exit("The proportion file is not a CSV file.")

# check output
if not os.path.isdir(args.output):
    sys.exit("The output directory does not exist.")


scriptpath = '/opt/CaDRReS-Sc'
sys.path.append(os.path.abspath(scriptpath))

from cadrres_sc import pp, model, evaluation, utility


## Predicting overall drug response and cell death percentage
# for each patient, if cell cluster is less than 5%, then we don't consider that cluster 
freq_cutoff = 0.05
# estimate cell death percentage based on log2 of the median IC50 observed in HNSC cell lines (GDSC)
ref_type = 'log2_median_ic50_hn'

## Read drug statistics
drug_info_df = pd.read_csv(scriptpath + '/preprocessed_data/GDSC/drug_stat.csv', index_col=0)
drug_info_df.index = drug_info_df.index.astype(str)
drug_id_name_dict = dict(zip(drug_info_df.index, drug_info_df['Drug Name']))

## Load cluster-specific drug response prediction
cadrres_cluster_df = pd.read_csv(args.input, header=[0,1], index_col=0)

#load prediction for a certain set of drugs
drug_list = drug_info_df.index
cluster_list = cadrres_cluster_df.index.astype(str)
drug_info_df = drug_info_df.loc[drug_list]
cadrres_cluster_df = cadrres_cluster_df[drug_list]

## Load cluster proportion information
freq_df = pd.read_csv(args.proportion, header=0, index_col=0)
patient_list = freq_df.index

## Predict cell death percentage at the ref_type dosage
pred_delta_df = pd.DataFrame(cadrres_cluster_df.values - drug_info_df[ref_type].values, columns=drug_list, index=cluster_list)
pred_cv_df = 100 / (1 + (np.power(2, -pred_delta_df)))
pred_kill_df = 100 - pred_cv_df

rows = []
print('List of cluster in each patient')
for p in patient_list:
    c_list = freq_df.loc[p][freq_df.loc[p] >= freq_cutoff].index.values
    freqs = freq_df.loc[p][freq_df.loc[p] >= freq_cutoff].values

    print(p, c_list, freqs)

    p_pred_delta_weighted = np.matmul(pred_delta_df.loc[c_list].values.T, freqs)
    p_pred_delta_mat = pred_delta_df.loc[c_list].values
    
    p_pred_kill_weighted = np.matmul(pred_kill_df.loc[c_list].values.T, freqs)
    p_pred_kill_mat = pred_kill_df.loc[c_list].values

    for d_i, d_id in enumerate(drug_list):
        rows += [[p, d_id, drug_id_name_dict[d_id]] + ['|'.join(c_list)] + ['|'.join(["{:.14}".format(f) for f in freqs])] + 
                  ['|'.join(["{:.14}".format(f) for f in p_pred_delta_mat[:, d_i]])] + 
                  ["{:.14}".format(p_pred_delta_weighted[d_i])] +
                  ['|'.join(["{:.14}".format(f) for f in p_pred_kill_mat[:, d_i]])] + 
                  ["{:.14}".format(p_pred_kill_weighted[d_i])]
                ]

single_drug_pred_df = pd.DataFrame(rows, columns=['patient', 'drug_id', 'drug_name', 'cluster', 'cluster_percentage', 'cluster_delta', 'delta', 'cluster_cell_death', 'cell_death'])
single_drug_pred_df.to_csv(os.path.join(args.output, 'drug_kill_prediction.csv'), index=False)
