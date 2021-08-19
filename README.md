# Single-Cell Analysis Pipeline

The Single-Cell Analysis Pipeline constructed a workflow for comprehensive analysis on single-cell RNA sequencing (scRNA-seq) data. It provided a powerful tool with various functions, from fundamental data analysis to drug response prediction, and treatment suggestions.

The Single-Cell Analysis Pipeline went through three parts on raw scRNA-seq data investigation: **Single-Cell Data Analysis**, **Drug Response Prediction**, and **Treatment Selection**.

- **Single-Cell Data Analysis** performed data preprocessing, clustering, cell type annotation and Gene Set Enrichment Analysis (GSEA). 

- **Drug Response Prediction** estimated the half maximal inhibitory concentration (IC50) of cell clusters, and reported the cell death percentages to drugs and drug combinations.

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

**Single-Cell Data Analysis** took the scRNA-seq data in a 10x-Genomics-formatted mtx directory or a CSV file as input, performed fundamental data analysis, and output a Scanpy Anndata object `scanpyobj.h5ad`, the UMAP `umap_cluster.png` and differentially expressed genes (DEGs) `cluster_DEGs.csv` of the clustering result, and a gene expression profile (GEP) file `GEP.txt`.

Optionally, **Single-Cell Data Analysis** carried out batch correction, cell type annotation and Gene Set Enrichment Analysis (GSEA), and provided additional UMAPs showing batch effects and cell type (`umap_batch.png` and `umap_cell_type.png`), and the GSEA result `GSEA_results.csv`. For cell type annotation, we used [scMatch: a single-cell gene expression profile annotation tool using reference datasets](https://github.com/asrhou/scMatch).

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

**Drug Response Prediction** examined  `scanpyobj.h5ad` generated in **Single-Cell Data Analysis**, and reported **IC50, single and combinatorial drug kill response** on specified clusters via [CaDRReS-Sc](https://github.com/CSB5/CaDRReS-SC) (a recommender system framework for *in silico* drug response prediction).

#### IC50 Prediction

**IC50 Prediction** took the Scanpy Anndata as input, and output cluster-specific `IC50_prediction.csv` as the prediction result.

- Run `python3 IC50_prediction.py -h` to show the help messages as follow for **IC50 Prediction**.

```
usage: IC50_prediction.py [-h] -i INPUT [-o OUTPUT] [-c CLUSTERS] [-m METHOD] [--hvg]

IC50 prediction

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        path to input Anndata object (h5ad file)
  -o OUTPUT, --output OUTPUT
                        path to output directory, default='./'
  -c CLUSTERS, --clusters CLUSTERS
                        perform IC50 prediction on specified clusters, e.g. '1,3,8,9', default='All'
  -m METHOD, --method METHOD
                        method used to produce clusterwise IC50 prediction, options: average | logfc, default: average
  --hvg                 only use highly variable genes to predict IC50
```

- Predict IC50 drug response on specified clusters with **IC50 Prediction**. For efficiency, run IC50 prediction with arguments `--metod logfc` or `--hvg`.

```
python3 IC50_prediction.py --input scanpyobj.h5ad --clusters CLUSTERS
```

#### Drug Kill Prediction

With cluster-specific `IC50_prediction.csv` and the cluster proportion information, **Drug Kill Prediction** computed cell death precentages `drug_kill_prediction.csv` based on the median IC50 observed in GDSC cell lines drug response.

- Run `python3 drug_kill_prediction.py -h` to show the help messages as follow for **Drug Kill Prediction**.

```
usage: drug_kill_prediction.py [-h] -i INPUT -p PROPORTION [-o OUTPUT]

Drug kill prediction

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        path to cluster-specific IC50 prediction (CSV file)
  -p PROPORTION, --proportion PROPORTION
                        path to cluster proportion CSV file
  -o OUTPUT, --output OUTPUT
                        path to output directory, default='./'
```

- Predict drug kill response with **Drug Kill Prediction**.

```
python3 drug_kill_prediction.py --input IC50_prediction.csv --proportion PROPORTION
```

#### Combinatorial Drug Kill Prediction

With `drug_kill_prediction.csv` and the drug list of interest, **Combinatorial Drug Kill Prediction** output `combinatorial_drug_kill_prediction.csv` storing cell death percentages to drug pairs within combinations in the drug list.

- Run `python3 combinatorial_drug_kill_prediction.py -h` to show the help messages as follow for **Combinatorial Drug Kill Prediction**.

```
usage: combinatorial_drug_kill_prediction.py [-h] -i INPUT -d DRUGS [-o OUTPUT]

Combinatorial drug kill prediction

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        path to cluster-specific IC50 prediction (CSV file)
  -d DRUGS, --drugs DRUGS
                        path to drugs list file
  -o OUTPUT, --output OUTPUT
                        path to output directory, default='./'
```

- Predict combinatorial drug kill response with **Drug Kill Prediction**.

```
python3 combinatorial_drug_kill_prediction.py --input drug_kill_prediction.csv --drugs DRUGS
```

### Treatment Selection

In **Treatment Selection**, we first **generated bulk GEPs** from the LINCS L1000 database and **imputed cell fractions** of `GEP.txt` created in **Single-Cell Data Analysis** via Docker version of [CIBERSORTx Cell Fractions](https://cibersortx.stanford.edu), which enumerated the proportions of distinct cell subpopulations in tissue expression profiles. Then, we **selected treatment combinations** from the LINCS L1000 database with the CIBERSORTx result.

#### Generate L1000 GEP

From `GEP.txt` created in **Single-Cell Data Analysis**, **Generate L1000 GEP** generated bulk GEP `LINCS_L1000_GEP_*.txt` from the LINCS L1000 database of a user-specified or automatically determined cell type of eight reference cell lines, including A375 (malignant melanoma), A549 (non-small cell lung carcinoma), HCC515 (non-small cell lung adenocarcinoma), HEPG2 (hepatocellular carcinoma), HT29 (colorectal adenocarcinoma), MCF7 (breast adenocarcinoma), PC3 (prostate adenocarcinoma), and YAPC (Pancreatic carcinoma).

- Run `python3 generate_L1000_GEP.py -h` to show the help messages as follow for **Generate L1000 GEP**.

```
usage: generate_L1000_GEP.py [-h] [-o OUTDIR] [-c CELLTYPE] [--inst INST] [--gctx GCTX] [--gene GENE] [--auto AUTO] [--gep GEP]

Generate bulk profiles of the specified cell type from the LINCS L1000 database.

optional arguments:
  -h, --help            show this help message and exit
  -o OUTDIR, --outdir OUTDIR
                        Path to output directory, default='./'
  -c CELLTYPE, --celltype CELLTYPE
                        Cell Line name. Options: A375 (malignant melanoma), A549 (non-small cell lung carcinoma), 
                        HCC515 (non-small cell lung adenocarcinoma), HEPG2 (hepatocellular carcinoma), HT29 (colorectal adenocarcinoma), 
                        MCF7 (breast adenocarcinoma), PC3 (prostate adenocarcinoma), YAPC (Pancreatic carcinoma)
  --inst INST           Path to inst_info file (.txt.gz)
  --gctx GCTX           Path to LINCS L1000 level 3 GEPs (.gctx)
  --gene GENE           Path to gene_info file (.txt.gz)
  --auto AUTO           Whether to automatically determine the reference cell type, default=False.
  --gep GEP             Path to GEP.txt for cell type determination.
```

- Generate LINCS L1000 bulk GEP from automatically determined cell type.

```
python3 generate_L1000_GEP.py --auto TRUE --gep GEP.txt
```

- Generate LINCS L1000 bulk GEP from sepcified cell type with information files.

```
python3 generate_L1000_GEP.py --celltype CELLTYPE --gctx GCTX --gene GENE --inst INST --gep GEP.txt
```

#### Impute Cell Fractions

With the single-cell and bulk GEPs, `GEP.txt` created in **Single-Cell Data Analysis** and `LINCS_L1000_GEP_*.txt` from the previous step, **Impute Cell Fractions** took a reference sample file from scRNA-seq data and a mixture matrix, both within the input directory, as input, to run CIBERSORTx Cell Fractions, and output CIBERSORTx result files to the output directory, including `CIBERSORTx_Results.txt`.

- Run `python3 CIBERSORTx_fractions.py -h` to show the help messages as follow for **Impute Cell Fractions**.

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

- **Impute Cell Fractions** with CIBERSORTx Cell Fractions.

```
python3 CIBERSORTx_fractions.py --input INPUT --username USERNAME --token TOKEN --refsample GEP.txt --mixture MIXTURE
```

Note: To obtain `USERNAME` and `TOKEN`, register and request for access to CIBERSORTx Docker on [CIBERSORTx](https://cibersortx.stanford.edu) website.

#### Select Treatment Combinations

**Select Treatment Combinations** took the CIBERSORTx result `CIBERSORTx_Results.txt` and the L1000 instance info file as input, selects treatment combinations for given cell type from the LINCS L1000 database, and outputs the treatment combinations list `CIBERSORTx_Results_solution_list_*.csv`.

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
python3 treatment_selection.py --input CIBERSORTx_Results.txt --celltype CELLTYPE --metadata METADATA
```
