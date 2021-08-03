# Single-Cell Analysis Pipeline

(Introduction)

## Download and Installation

1.  Clone the repository to local directory, ex: `./single-cell-analysis`.

    ```
    git clone https://github.com/b06902075/single-cell-analysis.git ./single-cell-analysis
    ```


2.  Build the docker image tagged `single-cell-analysis`.

    ```
    docker build -t single-cell-analysis ./single-cell-analysis
    ```


3.  Run the docker container named `single-cell-analysis` with `/docker/path` mounted to `/server/path` to access files within the docker container.
    
    ```
    docker run -it --name single-cell-analysis -v /server/path:/docker/path --privileged IMAGE_ID
    ```
    
    Note: Get `IMAGE_ID` with command `docker images`.
    
4.  In the docker container `single-cell-analysis`, pull the docker image `cibersortx/fractions` used in treatment selection.

    ```
    /etc/init.d/docker start
    docker pull cibersortx/fractions
    ```
    
    Note: Get `CONTAINER_ID` with command `docker ps -a` and start the container with `docker start -i CONTAINER_ID`.


## Usage

The single-cell analysis pipeline goes through fundamental data analysis, IC50 prediction and treatment selection on scRNA-seq data.

- **Single-Cell Data Analysis** performs clustering, cell type annotation and Gene Set Enrichment Analysis (GSEA) on the transcriptomic data. 

- **IC50 Prediction** estimates the half maximal inhibitory concentration (IC50) of given cell clusters.

- **Treatment Selection** lists treatment combinations of given cell clusters.


### Single-Cell Data Analysis

**Single-Cell Data Analysis** takes the scRNA-seq data in a 10x-Genomics-formatted mtx directory or a CSV file as input, and outputs a Scanpy Anndata object `scanpyobj.h5ad`, the UMAP `umap_cluster.png` and differentially expressed genes (DEGs) `cluster_DEGs.csv` of the clustering result, and a gene expression profile (GEP) file `GEP.txt`.

Optionally, **Single-Cell Data Analysis** carries out batch correction, cell type annotation and Gene Set Enrichment Analysis (GSEA), and provides additional UMAPs showing batch effects and cell type (`umap_batch.png` and `umap_cell_type.png`), and the GSEA result `GSEA_results.csv`.

Furthermore, **Single-Cell Data Analysis** can take previously produced Anndata as input and apply sub-clustering on specified clusters.


- Run `python single_cell_analysis.py -h` to show the help messages as follow for **Single-Cell Data Analysis**.

```
usage: single_cell_analysis.py [-h] -i INPUT [-f FORMAT] [-o OUTPUT] [-r RESOLUTION] [-m METADATA] [-b BATCH]
                               [-c CLUSTERS] [-a] [-g]

scRNA-seq data analysis

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        path to input 10x directory or CSV file
  -f FORMAT, --format FORMAT
                        input format, 10x (default) | csv | h5ad (Anndata object for subclustering with --clusters CLUSTERS)
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
  -g, --gsea            perform gene set enrichment analysis (GSEA)
```

- Apply **Single-Cell Data Analysis** with batch correction, cell type annotation and GSEA.

```
python single_cell_analysis.py --input INPUT --metadata METADATA --batch BATCH --annotation --gsea
```

- **Single-Cell Data Analysis** for sub-clustering with batch correction.

```
python single_cell_analysis.py --input scanpyobj.h5ad --batch BATCH
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
  --celltype CELLTYPE   Same as the cell type for decomposition. Options: A375|A549|ASC|BT20|CD34|HA1E|HCC515|HELA|HEPG2|
                        HME1|HS578T|HT29|HUES3|HUVEC|JURKAT|LNCAP|MCF10A|MCF7|MDAMB231|MNEU|NEU|NPC|PC3|SKBR3|SKL|YAPC
  --metadata METADATA   the L1000 instance info file, e.g., 'GSE70138_Broad_LINCS_inst_info_2017-03-06.txt'
```
