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
parser = argparse.ArgumentParser(description="Combinatorial drug kill prediction")

parser.add_argument("-i", "--input", required=True, help="path to cluster-specific IC50 prediction (CSV file)")
parser.add_argument("-d", "--drugs", required=True, help="path to drugs list file")
parser.add_argument("-o", "--output", default='./', help="path to output directory, default='./'")

args = parser.parse_args()

# check input
if not os.path.exists(args.input):
    sys.exit("The IC50 prediction file does not exist.")
if args.input[-4:] != '.csv':
    sys.exit("The IC50 prediction file is not a CSV file.")

# check drug
if not os.path.exists(args.drugs):
    sys.exit("The drugs list file does not exist.")

# check output
if not os.path.isdir(args.output):
    sys.exit("The output directory does not exist.")


scriptpath = '/opt/CaDRReS-Sc'
sys.path.append(os.path.abspath(scriptpath))

from cadrres_sc import pp, model, evaluation, utility


## Predict patient respose to combinatorial drugs
with open(args.drugs, 'r') as f:
    tested_drug_list = f.read().splitlines()

single_drug_pred_df = pd.read_csv(args.input)
single_drug_pred_df['drug_id'] = single_drug_pred_df['drug_id'].astype(str)
patient_list = list(set(single_drug_pred_df['patient']))

## Setup all drug combinations by patient
drug_combi_list = []
n_drugs = len(tested_drug_list)
for p in patient_list:
    for x in range(0, n_drugs-1):
        for y in range(x+1, n_drugs):
            drug_x = tested_drug_list[x]
            drug_y = tested_drug_list[y]

            drug_combi_list += [[p, drug_x, drug_y]]

drug_combi_df = pd.DataFrame(drug_combi_list, columns=['patient', 'A', 'B'])


merge_df = pd.merge(drug_combi_df, single_drug_pred_df, how='left', left_on=['patient', 'A'], right_on=['patient', 'drug_id'])
drug_combi_pred_df = pd.merge(merge_df, single_drug_pred_df[['patient', 'drug_id', 'drug_name', 'cluster_delta', 'delta', 'cluster_cell_death', 'cell_death']], how='left', left_on=['patient', 'B'], right_on=['patient', 'drug_id'], suffixes=['_A', '_B'])

rows = []
for _, data in drug_combi_pred_df.iterrows():
    
    cluster_p = np.array([float(p) for p in data['cluster_percentage'].split('|')])
    
    cluster_kill_A = np.array([float(k) for k in data['cluster_cell_death_A'].split('|')])
    cluster_kill_B = np.array([float(k) for k in data['cluster_cell_death_B'].split('|')])
    
    kill_A = float(data['cell_death_A'])
    kill_B = float(data['cell_death_B'])
    
    cluster_kill_C = cluster_kill_A + cluster_kill_B - np.multiply(cluster_kill_A/100, cluster_kill_B/100)*100
    kill_C = np.sum(cluster_p * cluster_kill_C)
    
    best_kill = np.max([kill_A, kill_B])
    improve = kill_C - best_kill
    improve_p = (kill_C - best_kill) / best_kill
    
    ##### specificity (entropy) #####
    
    temp_A = np.sum(cluster_p[cluster_kill_A > cluster_kill_B])
    temp_B = np.sum(cluster_p[cluster_kill_A <= cluster_kill_B])
    if temp_A == 0 or temp_B == 0:
        entropy = 0
    else:
        entropy = -(temp_A * np.log2(temp_A) + temp_B * np.log2(temp_B))
    
    sum_kill_dif = np.sum(np.abs(cluster_kill_A - cluster_kill_B))
    
    ##### save output #####
    
    rows += [['|'.join(["{:.14}".format(k) for k in cluster_kill_C])] + [kill_C, improve, improve_p, entropy, sum_kill_dif]]

drug_combi_pred_df = pd.concat([drug_combi_pred_df, pd.DataFrame(rows, columns=['cluster_cell_death_combi', 'cell_death_combi', 'improve', 'improve_percentage', 'kill_entropy', 'sum_kill_dif'])], axis=1)

## Final drug combination predictions for patients
drug_combi_pred_df = drug_combi_pred_df[['patient', 'drug_id_A', 'drug_name_A', 'drug_id_B', 'drug_name_B', 'cluster', 'cluster_percentage', 'cluster_cell_death_A', 'cluster_cell_death_B', 'cluster_cell_death_combi', 'cell_death_A', 'cell_death_B', 'cell_death_combi', 'improve', 'improve_percentage', 'kill_entropy', 'sum_kill_dif']]
drug_combi_pred_df.to_csv(os.path.join(args.output, 'combinatorial_drug_kill_prediction.csv'), index=False)