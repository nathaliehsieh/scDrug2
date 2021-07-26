#!/bin/bash
python ../scMatch-master/scMatch.py --refDS ../scMatch-master/refDB/FANTOM5 --dFormat csv --testDS ../data/harmony_scanpy_cluster_t-test_GSE156625_HCC_scMatch.csv --coreNum 2
# python ../scMatch-master/toTerms.py --splF ../data/harmony_scanpy_cluster_t-test_GSE156625_HCC_scMatch/annotation_result_keep_all_genes --refDS ../scMatch-master/refDB/FANTOM5 --coreNum 2
# python ../scMatch-master/visAnnos.py --dFormat csv --testDS ../data/harmony_scanpy_cluster_t-test_GSE156625_HCC_scMatch.csv --annoFile ../data/harmony_scanpy_cluster_t-test_GSE156625_HCC_scMatch/annotation_result_keep_all_genes/human_Spearman_top_ann.csv --visMethod u
# python ../scMatch-master/visAnnos.py --dFormat csv --testDS ../data/harmony_scanpy_cluster_t-test_GSE156625_HCC_scMatch.csv --annoFile ../data/harmony_scanpy_cluster_t-test_GSE156625_HCC_scMatch/annotation_result_keep_all_genes/human_Spearman_Avg_top_ann.csv --visMethod u
