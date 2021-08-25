# Single-Cell Analysis Pipeline

The Single-Cell Analysis Pipeline constructed a workflow for comprehensive analysis on single-cell RNA sequencing (scRNA-seq) data. It provided a powerful tool with various functions, from fundamental data analysis to drug response prediction, and treatment suggestions.

The Single-Cell Analysis Pipeline went through three parts on raw scRNA-seq data investigation: **Single-Cell Data Analysis**, **Drug Response Prediction**, and **Treatment Selection**.

- **Single-Cell Data Analysis** performed data preprocessing, clustering, cell type annotation and Gene Set Enrichment Analysis (GSEA). 

- **Drug Response Prediction** estimated the half maximal inhibitory concentration (IC50) of cell clusters, and reported the cell death percentages as drug kill efficacy.

- **Treatment Selection** lists treatment combinations of given cell clusters.


## Download and Installation

1.  Clone the repository to local directory, ex: `./single-cell-analysis`.

    ```
    git clone https://github.com/b06902075/single-cell-analysis.git ./single-cell-analysis
    ```


2.  Build the Docker image tagged `single-cell-analysis`.

    ```
    docker build -t single-cell-analysis ./single-cell-analysis
    ```


3.  Run the Docker container named `single-cell-analysis` with `/docker/path` mounted to `/server/path` to access files within the Docker container.
    
    ```
    docker run -it --name single-cell-analysis -v /server/path:/docker/path --privileged IMAGE_ID
    ```
    
    Note: Get `IMAGE_ID` with command `docker images`.
    
4.  In the Docker container `single-cell-analysis`, pull the Docker image `cibersortx/fractions` used in treatment selection.

    ```
    /etc/init.d/docker start
    docker pull cibersortx/fractions
    ```
    
    Note: Get `CONTAINER_ID` with command `docker ps -a` and start the container with `docker start -i CONTAINER_ID`.


## Usage

### Single-Cell Data Analysis

**Single-Cell Data Analysis** took the scRNA-seq data in a 10x-Genomics-formatted mtx directory or a CSV file as input, performed fundamental data analysis, and output a Scanpy Anndata object `scanpyobj.h5ad`, a UMAP `umap_cluster.png` and differentially expressed genes (DEGs) `cluster_DEGs.csv` of the clustering result, and a gene expression profile (GEP) file `GEP.txt`.

Optionally, **Single-Cell Data Analysis** carried out batch correction, cell type annotation and Gene Set Enrichment Analysis (GSEA), and provided additional UMAPs showing batch effects and cell types (`umap_batch.png` and `umap_cell_type.png`), and the GSEA result `GSEA_results.csv`. For cell type annotation, we used [scMatch: a single-cell gene expression profile annotation tool using reference datasets](https://github.com/asrhou/scMatch).

Furthermore, **Single-Cell Data Analysis** could take previously produced Anndata as input and applied sub-clustering on specified clusters.


- Run `python3 single_cell_analysis.py -h` to show the help messages as follow for **Single-Cell Data Analysis**.

```
usage: single_cell_analysis.py [-h] -i INPUT [-f FORMAT] [-o OUTPUT] [-r RESOLUTION] [--auto-resolution] [-m METADATA] [-b BATCH] [-c CLUSTERS] 
                               [--GEP GEP] [--annotation] [--gsea] [--cpus CPUS]

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
  --auto-resolution     automatically determine resolution for clustering
  -m METADATA, --metadata METADATA
                        path to metadata CSV file for batch correction (index as input in first column)
  -b BATCH, --batch BATCH
                        column in metadata (or adata.obs) for batch correction, e.g. 'PatientID'
  -c CLUSTERS, --clusters CLUSTERS
                        perform single cell analysis only on specified clusters, e.g. '1,3,8,9'
  --GEP GEP             whether to generate Gene Expression Profile file, default=True
  --annotation          perform cell type annotation
  --gsea                perform gene set enrichment analysis (GSEA)
  --cpus CPUS           number of CPU used for auto-resolution and annotation, default=1
```

- Apply **Single-Cell Data Analysis** with batch correction, clustering resolution 1.0, cell type annotation and GSEA.

```
python3 single_cell_analysis.py --input INPUT --metadata METADATA --batch BATCH --resolution 1.0 --annotation --gsea
```

- **Single-Cell Data Analysis** for sub-clustering with batch correction and automatically determined clustering resolution run under 2 cpus.

```
python3 single_cell_analysis.py --input scanpyobj.h5ad --batch BATCH --auto-resolution --cpus 2
```


### Drug Response Prediction

**Drug Response Prediction** examined  `scanpyobj.h5ad` generated in **Single-Cell Data Analysis**, reported IC50 and cell death percentages to drugs in GDSC database via [CaDRReS-Sc](https://github.com/CSB5/CaDRReS-SC) (a recommender system framework for *in silico* drug response prediction), and output the prediction result `drug_response_prediction.csv`.

- Run `python3 drug_response_prediction.py -h` to show the help messages as follow for **Drug Response Prediction**.

```
usage: drug_response_prediction.py [-h] -i INPUT [-o OUTPUT] [--hvg]

Drug response prediction

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        path to input Anndata object (h5ad file)
  -o OUTPUT, --output OUTPUT
                        path to output directory, default='./'
  --hvg                 only use highly variable genes to predict IC50
```

- Predict IC50 and drug kill response on input Anndata with **Drug Response Prediction**. For efficiency, run the following command with `--hvg`.

```
python3 drug_response_prediction.py --input scanpyobj.h5ad
```


### Treatment Selection

In **Treatment Selection**, we first **imputed cell fractions** of bulk GEPs from the LINCS L1000 database with single-cell GEP `GEP.txt` created in **Single-Cell Data Analysis** via Docker version of [CIBERSORTx Cell Fractions](https://cibersortx.stanford.edu), which enumerated the proportions of distinct cell subpopulations in tissue expression profiles. Then, we **selected treatment combinations** from the LINCS L1000 database with the CIBERSORTx result.

#### Impute Cell Fractions

**Impute Cell Fractions** took the input directory containing the reference sample file `GEP.txt` as input to run CIBERSORTx Cell Fractions with bulk GEP of user specified or automatically determined cell type, and output CIBERSORTx result files to the output directory, including `CIBERSORTx_Adjusted.txt`. The cell type for bulk GEP involved A375 (malignant melanoma),  A549 (non-small cell lung carcinoma),  HCC515 (non-small cell lung adenocarcinoma),  HEPG2 (hepatocellular carcinoma), HT29 (colorectal adenocarcinoma),  MCF7 (breast adenocarcinoma),  PC3 (prostate adenocarcinoma),  YAPC (Pancreatic carcinoma).

- Run `python3 CIBERSORTx_fractions.py -h` to show the help messages as follow for **Impute Cell Fractions**.

```
usage: CIBERSORTx_fractions.py [-h] -i INPUT [-o OUTPUT] -u USERNAME -t TOKEN -r REFSAMPLE
                               [--celltype CELLTYPE]

impute the fractions of previous identified cell subsets under each bulk sample in the LINCS L1000
database.

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
                        absolute path to the single-cell sample file
  --celltype CELLTYPE   choose a cell line from the options. If no name is provided, we will automatically
                        determine the cell type. Options: A375 (malignant melanoma), A549 (non-small cell
                        lung carcinoma), HCC515 (non-small cell lung adenocarcinoma), HEPG2 (hepatocellular
                        carcinoma), HT29 (colorectal adenocarcinoma), MCF7 (breast adenocarcinoma), PC3
                        (prostate adenocarcinoma), YAPC (Pancreatic carcinoma)
```

-  **Impute Cell Fractions** via CIBERSORTx Cell Fractions with single-cell GEP `GEP.txt` and LINCS L1000 bulk GEP of automatically determined cell type.

```
python3 CIBERSORTx_fractions.py --input INPUT --username USERNAME --token TOKEN --refsample GEP.txt
```

Note: To obtain `USERNAME` and `TOKEN`, register and request for access to CIBERSORTx Docker on [CIBERSORTx](https://cibersortx.stanford.edu) website.

#### Select Treatment Combinations

**Select Treatment Combinations** took the CIBERSORTx result `CIBERSORTx_Adjusted.txt` and the L1000 instance info file as input, selects treatment combinations for given cell type from the LINCS L1000 database, and outputs the treatment combinations list `CIBERSORTx_Results_solution_list_*.csv`.

- Run `python3 treatment_selection.py -h` to show the help messages as follow for **Select Treatment Combinations**.

```
usage: treatment_selection.py [-h] -i INPUT [-o OUTDIR] [-t THRESHOLD] [-c CON_THRESHOLD] --celltype CELLTYPE [--metadata METADATA]

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
  --celltype CELLTYPE   Same as the cell type for decomposition. Options: A375 | A549 | HEPG2 | HT29 | MCF7 | PC3 | YAPC
  --metadata METADATA   the L1000 instance info file, e.g., 'GSE70138_Broad_LINCS_inst_info_2017-03-06.txt'
```

- **Select Treatment Combinations** with the L1000 metadata.

```
python3 treatment_selection.py --input CIBERSORTx_Adjusted.txt --celltype CELLTYPE --metadata METADATA
```
