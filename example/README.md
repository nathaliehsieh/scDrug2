## Example Usage for The Single-Cell Analysis Pipeline

```
cd example
unzip data/10x_mtx.zip
unzip data/metadata.csv.zip
mkdir write
```

```
mkdir write/clustering

python3 ../script/single_cell_analysis.py --input data/10x_mtx --output write/clustering \
--metadata data/metadata.csv --batch PatientID --resolution 0.6 \
--annotation --gsea --GEP False
```

```
mkdir write/subclustering

python3 ../script/single_cell_analysis.py --input write/clustering/scanpyobj.h5ad --output write/subclustering \
--format h5ad  --clusters '1,5,9' --batch PatientID --resolution 0.8
```

```
mkdir write/drug_response_prediction

python3 ../script/drug_response_prediction.py --input write/subclustering/scanpyobj.h5ad --output write/drug_response_prediction
```

```
mkdir write/CIBERSORTx_fractions

python3 ../script/CIBERSORTx_fractions.py --input write/subclustering/GEP.txt --output write/CIBERSORTx_fractions \
--username USERNAME --token TOKEN --celltype HEPG2
```

```
mkdir write/treatment_selection

python3 ../script/treatment_selection.py --input data/CIBERSORTx_Results.txt --output write/treatment_selection \
--celltype HEPG2 --metadata data/GSE70138_Broad_LINCS_inst_info_2017-03-06.txt
```

```
mkdir write/draw_effect

python3 ../script/draw_effect.py --input write/treatment_selection --output write/draw_effect --drugs "palbociclib,NVP-BEZ235,selumetinib"
```
