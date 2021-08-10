from cmapPy.pandasGEXpress.parse import parse
import argparse
import pandas as pd
import numpy as np
import os
import sys


def downloadFromGEO(filename, url):
    print('downloading {} from the GEO website...'.format(filename))
    from urllib import request
    import shutil
    with request.urlopen(url) as response, open(filename, 'wb') as out_file:
        shutil.copyfileobj(response, out_file)
    if('gctx' in filename):
        import gzip
        with gzip.open(filename, 'rb') as f_in:
            with open(filename.rsplit('.',1)[0], 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description="Generate bulk profiles of the specified cell type from the LINCS L1000 database.")
    parser.add_argument("-o", "--outdir", default='./LINCS/', help="Path to output directory, default='./'")
    parser.add_argument("-c", "--celltype", default=None, help='Cell Line name. Options:  A375 (malignant melanoma),  A549 (non-small cell lung carcinoma),  HCC515 (non-small cell lung adenocarcinoma),  HEPG2 (hepatocellular carcinoma), HT29 (colorectal adenocarcinoma),  MCF7 (breast adenocarcinoma),  PC3 (prostate adenocarcinoma),  YAPC (Pancreatic carcinoma)')
    parser.add_argument("--inst", type=str, default='GSE70138_Broad_LINCS_inst_info_2017-03-06.txt.gz', help="Path to inst_info file (.txt.gz)")
    parser.add_argument("--gctx", type=str, default='GSE70138_Broad_LINCS_Level3_INF_mlr12k_n345976x12328_2017-03-06.gctx', help="Path to LINCS L1000 level 3 GEPs (.gctx)")
    parser.add_argument("--gene", type=str, default='GSE70138_Broad_LINCS_gene_info_2017-03-06.txt.gz', help="Path to gene_info file (.txt.gz)")
    parser.add_argument("--auto", type=bool, default=False, help="Whether to automatically determine the reference cell type, default=False.")
    parser.add_argument("--gep", type=str, help="Path to GEP.txt for cell type determination.")

    args = parser.parse_args()

    cell_types = ['A375','A549','HEPG2','HT29','MCF7','PC3','YAPC']
    do_detect = args.auto

    if not args.celltype and not do_detect:
         sys.exit("Please provide a specific cell type name or set --auto=TRUE. See --help for more information.")

    if args.celltype and (not args.celltype in cell_types):
        sys.exit('Unacceptable cell type: {}\nSee \'--help\' for acceptable cell line names. '.format(args.celltype))
    
    if do_detect and not os.path.isfile(args.gep):
        sys.exit("The single-cell gene expression profile ('GEP.txt') for cell type determination cannot be found. Check the file or set --auto to False to manually select the cell type.")

    if not os.path.isdir(args.outdir):
        os.mkdir(args.outdir)

    if not (os.path.isfile(args.inst) and args.inst.endswith('.gz')):
        if args.inst == 'GSE70138_Broad_LINCS_inst_info_2017-03-06.txt.gz':
            downloadFromGEO('GSE70138_Broad_LINCS_inst_info_2017-03-06.txt.gz', 'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE70nnn/GSE70138/suppl/GSE70138%5FBroad%5FLINCS%5Finst%5Finfo%5F2017%2D03%2D06%2Etxt%2Egz')
        else:
            sys.exit("The inst_info file does not exist or is not .txt. gz file.")
    if not (os.path.isfile(args.gctx) and args.gctx.endswith('.gctx')):
        if args.gctx == 'GSE70138_Broad_LINCS_Level3_INF_mlr12k_n345976x12328_2017-03-06.gctx':
            downloadFromGEO('GSE70138_Broad_LINCS_Level3_INF_mlr12k_n345976x12328_2017-03-06.gctx.gz', 'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE70nnn/GSE70138/suppl/GSE70138%5FBroad%5FLINCS%5FLevel3%5FINF%5Fmlr12k%5Fn345976x12328%5F2017%2D03%2D06%2Egctx%2Egz')
        else:
            sys.exit("The gctx file does not exist or is not .gctx file.")
    if not (os.path.isfile(args.gene) and args.gene.endswith('.gz')):
        if args.gene == 'GSE70138_Broad_LINCS_gene_info_2017-03-06.txt.gz':
            downloadFromGEO('GSE70138_Broad_LINCS_gene_info_2017-03-06.txt.gz', 'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE70nnn/GSE70138/suppl/GSE70138%5FBroad%5FLINCS%5Fgene%5Finfo%5F2017%2D03%2D06%2Etxt%2Egz')
        else:
            sys.exit("The gene_info file does not exist or is not .txt.gz file.")

    cell = args.celltype

    if do_detect:
        from scipy.stats import pearsonr

        average_gep = pd.read_csv(args.gep, delimiter='\t', index_col=0).mean(axis=1)
        cellline_gep = pd.read_csv('./LINCS/CCLE_GEP.csv', index_col=0)
        # calculate Pearson's correlation between average scGEP and the bulk GEP for each cell line
        mutual_genes = [x for x in average_gep.index if x in cellline_gep.index]
        max_p = float('-inf')
        for c in cellline_gep.columns:
            if((p:=pearsonr(average_gep[mutual_genes].to_numpy(), cellline_gep.loc[mutual_genes, c].to_numpy())[0]) > max_p):
                max_p = p
                cell = c
        cell = cell.split('_',1)[0]
        print('selected cell type = {}, with Pearsonr = {:.2f}'.format(cell, max_p))


    # check if the corresponding GEP already exists
    if(os.path.isfile('{}/LINCS_L1000_GEP_{}.txt'.format(args.outdir, cell))):
        print('GEP already exists. Location: {}/LINCS_L1000_GEP_{}.txt'.format(args.outdir, cell))
        exit(0)
    
    # read inst_info and gene_info
    inst_info = pd.read_csv(args.inst, sep='\t', compression='gzip')
    sig_info = pd.read_csv(args.gene, sep='\t', usecols=['pr_gene_id','pr_gene_symbol'], compression='gzip')

    # select instance ids for a specific cell type
    inst_ids = inst_info['inst_id'][inst_info['cell_id'] == cell]
    # read gctx
    gctoo = parse(args.gctx, cid=inst_ids)
    gctoo.data_df.index = gctoo.data_df.index.astype(int)
    # covert rowids to gene names
    named_df = pd.merge(gctoo.data_df, sig_info, left_index=True, right_on=['pr_gene_id'], validate='1:1')
    max_df = named_df.groupby('pr_gene_symbol').max().dropna().drop(labels='pr_gene_id',axis=1)
    # reverse to non-log for CIBERSORTx
    exp_df = 2**max_df
    exp_df.to_csv('{}/LINCS_L1000_GEP_{}.txt'.format(args.outdir, cell), sep='\t')