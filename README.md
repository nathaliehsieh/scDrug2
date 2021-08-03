# Single-Cell Analysis Pipeline

(Introduction)

## Download and Installation

```
git clone https://github.com/b06902075/single-cell-analysis.git ./single-cell-analysis
docker build -t single-cell-analysis ./single-cell-analysis
docker run -it --name single-cell-analysis -v /path/to/mount/in/server:/path/to/mount/in/docker --privileged IMAGE_ID 
/etc/init.d/docker start
docker pull cibersortx/fractions
```

## Usage

### Single-Cell Data Analysis

```
usage: single_cell_analysis.py [-h] -i INPUT [-f FORMAT] [-o OUTPUT] [-r RESOLUTION] [-m METADATA] [-b BATCH]
                               [-c CLUSTERS] [-a] [-g]

scRNA-seq data analysis

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        path to input 10x directory or CSV file
  -f FORMAT, --format FORMAT
                        input format, 10x (default) | csv | h5ad (Anndata object for subclustering with --clusters
                        CLUSTERS)
  -o OUTPUT, --output OUTPUT
                        path to output directory, default='./'
  -r RESOLUTION, --resolution RESOLUTION
                        resolution for clustering, default=0.6
  -m METADATA, --metadata METADATA
                        path to metadata CSV file for batch correction (index as input in first column)
  -b BATCH, --batch BATCH
                        column in metadata (or adata.obs) for batch correction, e.g. 'PatientID'
  -c CLUSTERS, --clusters CLUSTERS
                        perform single cell analysis only on specified clusters, e.g. '1,3,8,9'
  -a, --annotation      perform cell type annotation
  -g, --gsea            perform gene set enrichment analsis (GSEA)
```

### IC50 Prediction

```
usage: IC50_prediction.py [-h] -i INPUT [-o OUTPUT] [-c CLUSTERS]

IC50 prediction

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        path to input Anndata object (h5ad file)
  -o OUTPUT, --output OUTPUT
                        path to output directory, default='./'
  -c CLUSTERS, --clusters CLUSTERS
                        perform IC50 prediction on specified clusters, e.g. '1,3,8,9', default='All'
```

### Treatment Selection

```
usage: CIBERSORTx_fractions.py [-h] -i INPUT [-o OUTPUT] -u USERNAME -t TOKEN -r REFSAMPLE -m MIXTURE

impute cell fractions

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        path to input directory
  -o OUTPUT, --output OUTPUT
                        path to output directory, default='./'
  -u USERNAME, --username USERNAME
                        email address registered on CIBERSORTx website
  -t TOKEN, --token TOKEN
                        token obtained from CIBERSORTx website
  -r REFSAMPLE, --refsample REFSAMPLE
                        reference sample file from single cell RNA sequencing data
  -m MIXTURE, --mixture MIXTURE
                        mixture matrix required for running CIBERSORTx

```

```
usage: treatment_selection.py [-h] -i INPUT [-o OUTDIR] [-t THRESHOLD] [-c CON_THRESHOLD] --celltype CELLTYPE
                              [--metadata METADATA]

Select treatment combination from the LINCS L1000 database.

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        CIBERSORTx output file.
  -o OUTDIR, --outdir OUTDIR
                        path to output directory, default='./'
  -t THRESHOLD, --threshold THRESHOLD
                        Sensitivity threshold. Range: [-1,0), default:-0.9
  -c CON_THRESHOLD, --con_threshold CON_THRESHOLD
                        Consistency threshold. Range: [-1,0), default:-0.75
  --celltype CELLTYPE   Same as the cell type for decomposition. Options: A375|A549|ASC|BT20|CD34|HA1E|HCC515|HELA|HEP
                        G2|HME1|HS578T|HT29|HUES3|HUVEC|JURKAT|LNCAP|MCF10A|MCF7|MDAMB231|MNEU|NEU|NPC|PC3|SKBR3|SKL|Y
                        APC
  --metadata METADATA   the L1000 instance info file, e.g., 'GSE70138_Broad_LINCS_inst_info_2017-03-06.txt'
```
